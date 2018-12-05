import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from datetime import datetime, timedelta
import matplotlib.dates as dts
import pandas as pd
import os
import locale
import warnings

def visc(temp):
    """ Calculate viscosity of air """

    nyy_ref=18.203e-6
    S=110.4
    temp_ref=293.15
    nyy=nyy_ref*((temp_ref+S)/(temp+S))*((temp/temp_ref)**(3./2.))
    return nyy

def rlambda(t,press):
    """ Calculate mean-free path """

    kaasuv=8.3143
    dm=3.7e-10
    avoc=6.022e23
    return kaasuv*t/(np.sqrt(2.)*avoc*press*np.pi*dm*dm)

def cunn(Dp,t,press):
    """ Calculate Cunningham correction """

    lambd = rlambda(t,press)
    return 1. + 2.*lambd/Dp * (1.165 + 0.483 \
           * np.exp(-0.997*Dp/lambd))

def diffuus(dpp,temp,press):
    """ Calculate diffusion coefficient """

    K=1.38e-23
    return (K*temp*cunn(dpp,temp,press))/(3.*np.pi*visc(temp)*dpp)

def tubeloss(dpp,pflow,plength,temp,press):
    """ Laminar diffusion losses in circular straight conduit """

    rmuu = np.pi*diffuus(dpp,temp,press)*plength/pflow;
    pene = np.zeros(rmuu.shape)
    for i in range(0,len(dpp)):
        if rmuu[i] < 0.02:
            pene[i]=1. - 2.56*rmuu[i]**(2./3.) \
                    + 1.2*rmuu[i]+0.177*rmuu[i]**(4./3.)
        else:
            pene[i]=0.819*np.exp(-3.657*rmuu[i])+0.097 \
                    *np.exp(-22.3*rmuu[i])+0.032*np.exp(-57.*rmuu[i])
    return pene

def x2dlogx(x):
    """ Calculate log-differences """

    logx = np.log10(x)
    logx_mid = (logx[1:] + logx[:-1])/2.0 
    logx_mid_first_value = logx[0] + (logx[0] - logx_mid[0])
    logx_mid_last_value  = logx[-1] + (logx_mid[-1] - logx[-1])
    logx_mid = np.insert(logx_mid,0,logx_mid_first_value)
    logx_mid = np.append(logx_mid,logx_mid_last_value)
    dlogx = np.abs(np.diff(logx_mid))
    return dlogx

def datetime2datenum(dt):
    """ Convert from python datetime to matlab datenum """

    ord = dt.toordinal()
    mdn = dt + timedelta(days = 366)
    frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds \
           / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac

def datenum2datetime(matlab_datenum):
    """ Convert from matlab datenum to python datetime """

    return datetime.fromordinal(int(matlab_datenum)) \
    + timedelta(days=matlab_datenum%1) \
    - timedelta(days = 366)


def plot_sumfile(handle,v,clim=(10,100000)):
    """ Plot UHEL's sum-formatted aerosol number-size distribution """
    
    plt.rcParams['image.cmap'] = 'viridis'
    plt.rcParams.update({'font.size': 10,
                         'axes.titlesize': 10,
                         'figure.titlesize': 10,
                         'legend.fontsize': 10})
    
    time = v[1:,0]
    dp = v[0,2:]
    data = v[1:,2:]
    mesh_dp,mesh_time = np.meshgrid(dp,time)
    pcolorplot = handle.pcolormesh(mesh_time,mesh_dp,data,
                                   norm=colors.LogNorm(),
    linewidth=0,rasterized=True)
    handle.set_yscale('log')
    pcolorplot.set_clim(clim)
    pcolorplot.set_edgecolor('face')
    handle.autoscale(tight='true')

    handle.grid('on',which='both',linestyle='--',color='w',lw=0.5)
    handle.xaxis.set_major_locator(dts.MinuteLocator(interval=120))
    handle.xaxis.set_major_formatter(dts.DateFormatter('%H'))
    plt.setp(handle.get_xticklabels(), rotation=60)
    handle.set_ylabel('Dp, [m]')
    handle.set_xlabel('Time, [h]')

    cbar = plt.colorbar(pcolorplot,ax=handle,
                        ticks=LogLocator(subs=range(10)))
    cbar.set_label('dN/dlogDp, [cm-3]')
    return pcolorplot
    


def nais_processor(load_path='../data',
                   save_path='../processed',
                   ideal_flow=27.0,
                   particle_files='block-particles.spectra',
                   ion_files='block-ions.spectra',
                   records_files='block.records',
                   delimiter='\t',
                   pipe_length=0.0):
    """ Function that processes data from the (N)AIS 

    The function reads the raw data files from (Neutral cluster) and 
    Air Ion Spectrometer (N)AIS, applies corrections (diffusion losses, 
    temperature/pressure changes, sample flow changes, Robert Wagner's 
    ion mode calibration) to the measured number concentrations and 
    saves the data as a University of Helsinki sum-formatted number-size
    distribution

    Note:
        Sum-format:

            n[1:,0] = Time (days) 
            n[0,2:] = Geometric mean diameter of size-channel (m)
            n[1:,2:] = Normalized number concentration (cm-3)

    Args:
        path (str): path to data files
        ideal_flow (float): total sample flow / 2.0
        particle_files (str): string that identifies the particle files
        ion_files (str): string that identifies the ion files
        recrods_files (str): string that identifies the records files
        delimiter (str): delimiter in the files 
        pipe_length (float): length of the inlet tube measured from the 
                             front of the instrument

    """
    # Ignore all warnings
    warnings.filterwarnings("ignore")
    
    # Change the temperature and pressure to fit ambient conditions if 
    # necessary
    temp = 293.15 # K  
    pres = 101325.0 # Pa 

    # Standard conditions
    temp_ref = 273.15 # K
    pres_ref = 101325.0 # Pa

    # Load mobilities and diameters
    dp_ion_lowres = np.loadtxt('dp_ion_lowres.dat') # nm
    dp_ion_hires = np.loadtxt('dp_ion_hires.dat') # nm
    dp_par_lowres = np.loadtxt('dp_par_lowres.dat') # nm
    dp_par_hires = np.loadtxt('dp_par_hires.dat') # nm
    mob_ion_lowres = np.loadtxt('mob_ion_lowres.dat') # cm2 s-1 V-1
    mob_ion_hires = np.loadtxt('mob_ion_hires.dat') # cm2 s-1 V-1    

    # Find the filenames in the destination directory and sort them
    filenames = os.listdir(load_path)
    filenames.sort()

    # List particle/ion/records filenames
    particle_filenames = [x for x in filenames if particle_files in x] 
    ion_filenames = [x for x in filenames if ion_files in x]
    records_filenames = [x for x in filenames if records_files in x]

    # Determine what quantities were measured
    ions_measured = 1 if len(ion_filenames)>0 else 0
    records_measured = 1 if len(records_filenames)>0 else 0
    particles_measured = 1 if len(particle_filenames)>0 else 0

    number_of_data_days = len(ion_filenames)

    if ions_measured==0:
        raise Exception("Need at least ion data")

    # Process the last day found in the data
    for i in range(0,number_of_data_days):
        # Read ion data
        ions = pd.read_table(load_path+'/'+ion_filenames[i],
                             sep=delimiter,
                             comment='#',
                             engine='python')

        # Determine inverter type
        if len(ions.columns)>400:
            dp_ion = dp_ion_hires
            mob_ion = mob_ion_hires
            dp_par = dp_par_hires
        else:
            dp_ion = dp_ion_lores                                   
            mob_ion = mob_ion_lores                                 
            dp_par = dp_par_lores

        # Calculate log-differences
        dlogdp_ion = x2dlogx(dp_ion)
        dlogmob_ion = x2dlogx(mob_ion)
        dlogdp_par = x2dlogx(dp_par)

        # Extract ion number-size distribution
        ions = ions.set_index(ions.columns[0])
        ions.index = ions.index.map(pd.to_datetime)
        neg_ions = ions.iloc[:,3:3+len(dp_ion)].interpolate().values
        pos_ions = ions.iloc[:,3+len(dp_ion)+len(dp_ion):3+\
                            len(dp_ion)+len(dp_ion)+len(dp_ion)]\
                            .interpolate().values
        neg_ions = neg_ions * dlogmob_ion / dlogdp_ion
        pos_ions = pos_ions * dlogmob_ion / dlogdp_ion

        if records_measured:
            # Read the diagnostic data
            records = pd.read_table(load_path+'/'+records_filenames[i],
                                    sep=delimiter,
                                    comment='#',
                                    engine='python')
            records = records.set_index('opmode')

            # Extract the diagnostic data during the ion measurements
            ion_records = records.loc['ions']\
                                          .set_index(records.columns[0])
            ion_records.index = ion_records.index.map(pd.to_datetime)
            ion_records = ion_records.reindex(index=ions.index)

            # Find the names of diagnostic variables that will be used
            # in the corrections
            if i==0:
                print str(list(records)).replace(',','\n')[1:-1]
                print 'FIND THE NAME OF THE BELOW QUANTITIES FROM THE\n\
                       LIST AND PRESS ENTER. LEAVE EMPTY IF THE\n\
                       VARIABLE CANNOT BE FOUND'
                temp_name = raw_input ("TEMPERATURE: ")
                pres_name = raw_input ("AIR PRESSURE: ")
                flow_name_pos = raw_input ("SAMPLE FLOW (-): ")
                flow_name_neg = raw_input ("SAMPLE FLOW (+): ")

#            # Alternatively you can list the variable names here
#            temp_name='temperature.mean'
#            pres_name='baro.mean'
#            flow_name_pos='pos_sampleflow.mean'
#            flow_name_neg='neg_sampleflow.mean'

            temperatures = 1 if temp_name!='' else 0
            pressures = 1 if pres_name!='' else 0
            flows_neg = 1 if flow_name_neg!='' else 0
            flows_pos = 1 if flow_name_pos!='' else 0

        # Extract temperatures, pressures and sample flows
        if temperatures & pressures:
            temp_ions = (273.15 + ion_records[temp_name]\
                        .interpolate()).values.flatten()
            pres_ions = (100.0 * ion_records[pres_name]\
                        .interpolate()).values.flatten()
        else:
            temp_ions = temp*np.ones(neg_ions.shape[0])
            pres_ions = pres*np.ones(neg_ions.shape[0])
        
        if (flows_neg & flows_pos):
            flow_ions_neg = (ion_records[flow_name_neg].interpolate())\
                             .values.flatten()
            flow_ions_pos = (ion_records[flow_name_pos].interpolate())\
                             .values.flatten()
        else:
            flow_ions_neg = ideal_flow*np.ones(neg_ions.shape[0])
            flow_ions_pos = ideal_flow*np.ones(pos_ions.shape[0])

        # flow correction
        flow_corr_ions_neg = (flow_ions_neg/ideal_flow)
        flow_corr_ions_pos = (flow_ions_pos/ideal_flow)
        neg_ions = (flow_corr_ions_neg*neg_ions.T).T
        pos_ions = (flow_corr_ions_pos*pos_ions.T).T

        # Pressure and temperature correction
        stp_corr_ions = (pres_ref*temp_ions)/(temp_ref*pres_ions)
        neg_ions = (stp_corr_ions*neg_ions.T).T
        pos_ions = (stp_corr_ions*pos_ions.T).T

        # Diffusion loss correction
        throughput_ions = np.zeros(neg_ions.shape)
        for j in range(0,throughput_ions.shape[0]):
            throughput_ions[j,:] = tubeloss(dp_ion*1e-9,
                           (flow_ions_neg[j]+flow_ions_pos[j])\
                           *1e-6,pipe_length,temp_ions[j],pres_ions[j])
        neg_ions = throughput_ions*neg_ions
        pos_ions = throughput_ions*pos_ions

        # Robert Wagner's calibration (only ions)
        roberts_corr = 0.713*dp_ion**0.120
        neg_ions = roberts_corr*neg_ions
        pos_ions = roberts_corr*pos_ions

        # Integrate total ion number concentration
        total_neg_ions = np.nansum(neg_ions*dlogdp_ion,axis=1)
        total_pos_ions = np.nansum(pos_ions*dlogdp_ion,axis=1)

        # Define time vector for ions
        time_ions = np.array([datetime2datenum(x) for x in ions.index])
        
        # Construct the sum matrices
        negions = np.concatenate((np.insert(dp_ion*1e-9,(0,0),0)\
            [np.newaxis],np.concatenate((time_ions[np.newaxis].T,\
            total_neg_ions[np.newaxis].T,neg_ions),axis=1)))
        posions = np.concatenate((np.insert(dp_ion*1e-9,(0,0),0)\
            [np.newaxis],np.concatenate((time_ions[np.newaxis].T,\
            total_pos_ions[np.newaxis].T,pos_ions),axis=1)))
        
        # Save the sum matrices
        np.savetxt(save_path+'/'+ion_filenames[i]\
                   +'_processed_neg_ions.sum',negions)
        np.savetxt(save_path+'/'+ion_filenames[i]\
                   +'_processed_pos_ions.sum',posions)

        # Apply same processing to particles if they exist
        if particles_measured:
            particles = pd.read_table(load_path+'/'+particle_filenames[i],
                                 sep=delimiter,
                                 comment='#',
                                 engine='python')

            particles = particles.set_index(particles.columns[0])
            particles.index = particles.index.map(pd.to_datetime)
            neg_particles = particles.iloc[:,3:3+len(dp_par)]\
                              .interpolate().values
            pos_particles = particles.iloc[:,3+len(dp_par)+len(dp_par)\
                               :3+len(dp_par)+len(dp_par)+len(dp_par)]\
                               .interpolate().values
            neg_particles = neg_particles
            pos_particles = pos_particles

            if records_measured:
                particle_records = records.loc['particles']\
                                         .set_index(records.columns[0])
                particle_records.index = particle_records\
                                         .index.map(pd.to_datetime)
                particle_records = particle_records\
                                        .reindex(index=particles.index)

            if temperatures & pressures:
                temp_particles = (273.15 + particle_records[temp_name]\
                            .interpolate()).values.flatten()
                pres_particles = (100.0 * particle_records[pres_name]\
                            .interpolate()).values.flatten()
            else:
                temp_particles = temp*np.ones(neg_particles.shape[0])
                pres_particles = pres*np.ones(neg_particles.shape[0])
            
            if (flows_neg & flows_pos):
                flow_particles_neg = (particle_records[flow_name_neg]\
                                     .interpolate()).values.flatten()
                flow_particles_pos = (particle_records[flow_name_pos]\
                                     .interpolate()).values.flatten() 
            else:
                flow_particles_neg = ideal_flow*np.ones(neg_particles\
                                                     .shape[0])
                flow_particles_pos = ideal_flow*np.ones(pos_particles\
                                                     .shape[0])

            flow_corr_particles_neg = (flow_particles_neg/ideal_flow)
            flow_corr_particles_pos = (flow_particles_pos/ideal_flow)
            neg_particles = (flow_corr_particles_neg*neg_particles.T).T
            pos_particles = (flow_corr_particles_pos*pos_particles.T).T

            stp_corr_particles = (pres_ref*temp_particles)\
                                 / (temp_ref*pres_particles)
            neg_particles = (stp_corr_particles*neg_particles.T).T
            pos_particles = (stp_corr_particles*pos_particles.T).T

            throughput_particles = np.zeros(neg_particles.shape)
            for j in range(0,throughput_particles.shape[0]):
                throughput_particles[j,:] = tubeloss(dp_par*1e-9,\
                   (flow_particles_neg[j]+flow_particles_pos[j])*1e-6,\
                    pipe_length,temp_particles[j],pres_particles[j])
            neg_particles = throughput_particles*neg_particles
            pos_particles = throughput_particles*pos_particles

            total_neg_particles = np.nansum(neg_particles*dlogdp_par,
                                            axis=1)
            total_pos_particles = np.nansum(pos_particles*dlogdp_par,
                                            axis=1)

            time_particles = np.array([datetime2datenum(x)\
                                             for x in particles.index])
            
            negparticles = np.concatenate((np.insert(dp_par*1e-9,(0,0)\
            ,0)[np.newaxis],np.concatenate((time_particles[np.newaxis]\
            .T,total_neg_particles[np.newaxis].T,neg_particles),axis=1\
            )))
            posparticles = np.concatenate((np.insert(dp_par*1e-9,(0,0)\
            ,0)[np.newaxis],np.concatenate((time_particles[np.newaxis]\
            .T,total_pos_particles[np.newaxis].T,pos_particles),axis=1\
            )))
            
            np.savetxt(save_path+'/'+particle_filenames[i]
                          +'_processed_neg_particles.sum',negparticles)
            np.savetxt(save_path+'/'+particle_filenames[i]
                          +'_processed_pos_particles.sum',posparticles)

        print "Files ready"

def nais_plotter(load_path='../processed',
                 save_path='../figures'):
   """ Function to plot the processed NAIS data """
  
   plt.style.use('dark_background')

   # list all processed data and sort it
   filenames = os.listdir(load_path)
   filenames.sort()

   # Extract the list of different kinds of data files
   pos_particle_filenames = [x for x in filenames if 'pos_particles' \
                            in x]
   neg_particle_filenames = [x for x in filenames if 'neg_particles' \
                            in x]
   pos_ion_filenames = [x for x in filenames if 'pos_ions' in x]
   neg_ion_filenames = [x for x in filenames if 'neg_ions' in x]

   num_files = len(pos_ion_filenames)

   particles_measured = 1 if len(neg_particle_filenames)>0 else 0

   # Take the second lastest files and make a plot
   for i in range(0,num_files):
       if particles_measured:
           locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

           pos_ion_data = np.loadtxt(load_path+'/'+pos_ion_filenames[i])
           neg_ion_data = np.loadtxt(load_path+'/'+neg_ion_filenames[i])
           pos_particle_data = np.loadtxt(load_path+'/'\
                                          +pos_particle_filenames[i])
           neg_particle_data = np.loadtxt(load_path+'/'\
                                          +neg_particle_filenames[i])

           date = datenum2datetime(pos_ion_data[1,0])
           date_str =  date.strftime("%b")\
                       +'-'+date.strftime("%d")\
                       +'-'+date.strftime('%Y')

           fig,ax = plt.subplots(2,2,figsize=(9,6))
           ax = ax.flatten()

           plot_sumfile(ax[0],pos_ion_data,clim=(10,10000))
           plot_sumfile(ax[1],neg_ion_data,clim=(10,10000))
           plot_sumfile(ax[2],pos_particle_data,clim=(10,100000))
           plot_sumfile(ax[3],neg_particle_data,clim=(10,100000))
                       
           ax[0].set_title('Positive ions')
           ax[1].set_title('Negative ions')
           ax[2].set_title('Particles (positive polarity)')
           ax[3].set_title('Particles (negative polarity)')

           ax[2].set_ylim(bottom=2.0e-9)
           ax[3].set_ylim(bottom=2.0e-9)

           plt.tight_layout()
           
           fig.suptitle(date_str, y=1)
           
           plt.savefig(save_path+'/'+'NAIS_'+date_str,
                       bbox_inches='tight',dpi=300)
           plt.close()
           print "Plot ready"
