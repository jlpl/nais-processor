import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from datetime import datetime, timedelta
from scipy.optimize import minimize
import matplotlib.dates as dts
import re
import easygui
import pandas as pd
import os
import errno

# J. Lampilahti 28.11.2017

################################################################### FUNCTIONS

def set_plotting_defaults(**kwargs):
    """
    Enables LaTeX
    Choose fontsize and default colormap
    e.g. set_plotting_defaults(fontsize=15,cmap='jet')
    """
    plt.rcParams['image.cmap'] = kwargs['cmap']
    plt.rcParams.update({'text.usetex': True,
                     'font.size': kwargs['fontsize'],
                     'axes.titlesize': kwargs['fontsize'],
                     'legend.fontsize': kwargs['fontsize'],
                     'font.family':'Times',
                     'text.latex.unicode' : True,
                     'text.latex.preamble' : [r'\usepackage{mathptmx}']})
    return

def plot_sumfile2(handle,time,dp,data):
    """
    Plot number size distribution given in sum-matrix format
    """
    mesh_dp,mesh_time = np.meshgrid(dp,time)
    pcolorplot = handle.pcolormesh(mesh_time,mesh_dp,data,norm=colors.LogNorm(),linewidth=0)
    handle.set_yscale('log')
    pcolorplot.set_edgecolor('face')
    handle.autoscale(tight='true')
    cbar = plt.colorbar(pcolorplot,ax=handle,ticks = LogLocator(subs=range(10)))
    return (pcolorplot,cbar)

def tubeloss(Dp,pipeflow,pipediameter,pipelength,t,press):
    """
    Laminar diffusion losses in circular straight conduit
    Translated from code by T. Petaja
    """
    def visc(temp):
       nyy_ref=18.203e-6
       S=110.4
       temp_ref=293.15
       nyy=nyy_ref*((temp_ref+S)/(temp+S))*((temp/temp_ref)**(3./2.))
       return nyy
       
    def cunn(Dp,t,press):
        lambd = rlambda(t,press)
        return 1. + 2.*lambd/Dp * (1.165 + 0.483 * np.exp(-0.997*Dp/lambd))
        
    def rlambda(t,press):
        kaasuv=8.3143
        dm=3.7e-10
        avoc=6.022e23
        return kaasuv*t/(np.sqrt(2.)*avoc*press*np.pi*dm*dm)

    def diffuus(dpp,temp,press):
        K=1.38e-23
        return (K*temp*cunn(dpp,temp,press))/(3.*3.14*visc(temp)*dpp)
        
    def reynolds(pipeflow,pipediameter,t,press):
        density=1.29*(273.15/t)*(press/101325.)
        pipearea=np.pi/4.*pipediameter**2
        velocity=pipeflow/pipearea
        return density*velocity*pipediameter/visc(t)
        
    def ltubefl(dpp,plength,pflow,temp,press):
        rmuu=np.pi*diffuus(dpp,temp,press)*plength/pflow;
        res = np.zeros(rmuu.shape)
        for i in range(0,len(dpp)):
            if rmuu[i] < 0.02:
                res[i]=1.-2.56*rmuu[i]**(2./3.)+1.2*rmuu[i]+0.177*rmuu[i]**(4./3.)
            else:
                res[i]=0.819*np.exp(-3.657*rmuu[i])+0.097*np.exp(-22.3*rmuu[i])+0.032*np.exp(-57.*rmuu[i])
        return res

    difflam = ltubefl(Dp,pipelength,pipeflow,t,press)

    return difflam

def datetime2datenum(dt):
    """
    Convert from python datetime to matlab datenum
    """
    ord = dt.toordinal()
    mdn = dt + timedelta(days = 366)
    frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac

def mkdir_p(path):
    """
     Analogous to UNIX mkdir -p
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def Cc(Dp):
    """ 
    Makela et al. (1996)
    """
    return 1.+2.*64.5/Dp*(1.246+0.420*np.exp(-0.87*Dp/(2.*64.5))) 

def diam_to_mob(Dp): # [Dp] = nm
    """
    Electrical mobility diameter -> Electrical mobility
    """
    e = 1.60217662e-19 # Coulomb
    return (e*Cc(Dp))/(3.*np.pi*1.83245e-5*Dp*1e-9) # m2 s-1 V-1

def mob_to_diam(Zp):
    """
    Electrical mobility -> electrical mobility diameter
    """
    def minimize_this(Dp,Zp): # [Dp] = nm, [Zp] = m2 s-1 V-1
        return np.abs(diam_to_mob(Dp)-Zp)
    Dp0 = 1. # initial guess in nm
    return minimize(minimize_this, Dp0, args=(Zp,), tol=1e-15).x[0] # m

###################################################################### END FUNCTIONS

set_plotting_defaults(fontsize=11,cmap='jet')

title = "NAIS processer"

# Enter measurement information (GUI)
field_names = ["Path to files",
               "Start date (yyyymmdd)",
               "End date (yyyymmdd)",
               "Inverter resolution (lores/hires)",
               "Sampleflow (lpm)",
               "STP correction (yes/no)",
               "Tubeloss correction (yes/no)",
               "Inlet length (m)",
               "Inlet diameter (m)",
               "Delimiter in files (\\t/,/etc.)",
               "Date format in filenames (e.g. %Y%m%d or %Y-%m-%d)",
               "Records/diagnostics filename ending",
               "Ion filename ending",
               "Particle filename ending",
               "Temperature/pressure recorded? (yes/no)",
               "Deviation from UTC in hours (negative or positive integer)"]

# These are the default values, change if needed
field_values = ['/path/to/files/',
                '20170524',
                '20170716',
                'lores',
                '54.0',
                'no',
                'yes',
                '0.4',
                '0.035',
                '\t',
                '%Y%m%d',
                '-block.records',
                '-block-ions.spectra',
                '-block-particles.spectra',
                'yes',
                '0']

field_values = easygui.multenterbox("Enter information", title, field_names, values = field_values)

# Assign values to variables
path = field_values[0]
date1 = field_values[1]
date2 = field_values[2]
inv = field_values[3]
stp_corr = field_values[5]
tubeloss_corr = field_values[6]
delimiter = field_values[9]
datefmt = field_values[10]
recordfmt = field_values[11]
ionfmt = field_values[12]
particlefmt = field_values[13]
tp_data = field_values[14]
sf = float(field_values[4])   # lpm convert to m3/s
inlet_length = float(field_values[7])
inlet_diameter = float(field_values[8])
time_zone = float(field_values[15])

# Create list of dates
datetime1 = datetime.strptime(date1,'%Y%m%d')
datetime2 = datetime.strptime(date2,'%Y%m%d')
datetimes = pd.date_range(start=datetime1,end=datetime2).tolist()

# Reference temperature and pressure
temp_ref=293.15     # K
pres_ref=101325.0   # Pa

# Temperature and pressure conditions at the instrument.
# These will be used in corrections if the instrument does not record temp/press
# Disable tubeloss/STP corrections if you want to correct later with own data.
temp_inst=293.15     # K
pres_inst=101325.0   # Pa

# Create Figure and Data directories if they do not exist
mkdir_p(path+'figures')
mkdir_p(path+'sum')

# Low resolution mobilities (ions) in cm2 s-1 V-1
mob_lores=1e-4*np.array([3.16227766016838, 2.37137370566166, 1.77827941003892, 1.33352143216332, 1, 0.749894209332456, 0.562341325190349, 0.421696503428582, 0.316227766016838, 0.237137370566166, 0.177827941003892, 0.133352143216332, 0.1, 0.0749894209332456, 0.0562341325190349, 0.0421696503428583, 0.0316227766016838, 0.0237137370566166, 0.0177827941003892, 0.0133352143216333, 0.01, 0.00749894209332457, 0.0056234132519035, 0.00421696503428583, 0.00316227766016838, 0.00237137370566166, 0.00177827941003892, 0.00133352143216333])

# Low resolution diameters (particles) in nm
dps_lores=2.*np.array([0.374947104666227, 0.432982161680032, 0.499999999999999, 0.577390992344728, 0.666760716081661, 0.769963263029745, 0.88913970501946, 1.02676251322857, 1.18568685283083, 1.36920981713218, 1.58113883008419, 1.82587063627419, 2.10848251714291, 2.43483762582931, 2.81170662595174, 3.24690815788105, 3.74947104666228, 4.32982161680032, 5, 5.77390992344729, 6.66760716081662, 7.69963263029745, 8.89139705019461, 10.2676251322857, 11.8568685283083, 13.6920981713218, 15.8113883008419, 18.2587063627419, 21.0848251714291])

# High resolution mobilities (ions)
mob_hires=1e-4*np.array([3.162277660168379, 2.942727176209282, 2.738419634264362, 2.5482967479793475, 2.371373705661656, 2.206734069084591, 2.053525026457147, 1.9109529749704417, 1.7782794100389239, 1.6548170999431828, 1.5399265260594934, 1.4330125702369643, 1.3335214321633255, 1.2409377607517211, 1.1547819846894598, 1.0746078283213192, 1.0000000000000016, 0.9305720409297006, 0.8659643233600669, 0.8058421877614834, 0.7498942093324574, 0.6978305848598679, 0.6493816315762129, 0.6042963902381343, 0.5623413251903505, 0.5232991146814961, 0.48696752516586445, 0.45315836376008306, 0.4216965034285835, 0.39241897584845475, 0.36517412725483883, 0.339820832894257, 0.31622776601683894, 0.29427271762092916, 0.2738419634264371, 0.25482967479793556, 0.23713737056616638, 0.2206734069084598, 0.2053525026457154, 0.1910952974970448, 0.177827941003893, 0.16548170999431883, 0.15399265260594985, 0.1433012570236969, 0.13335214321633299, 0.12409377607517251, 0.11547819846894634, 0.10746078283213226, 0.10000000000000048, 0.09305720409297036, 0.08659643233600697, 0.08058421877614859, 0.07498942093324598, 0.069783058485987, 0.06493816315762148, 0.06042963902381361, 0.056234132519035224, 0.05232991146814976, 0.04869675251658659, 0.04531583637600844, 0.042169650342858474, 0.03924189758484559, 0.03651741272548399, 0.033982083289425806, 0.03162277660168399, 0.029427271762093, 0.027384196342643784, 0.025482967479793627, 0.023713737056616703, 0.022067340690846038, 0.020535250264571595, 0.01910952974970453, 0.017782794100389347, 0.016548170999431927, 0.015399265260595025, 0.014330125702369726, 0.013335214321633334, 0.012409377607517284, 0.011547819846894663, 0.01074607828321325, 0.010000000000000071, 0.009305720409297056, 0.008659643233600717, 0.008058421877614878, 0.007498942093324614, 0.006978305848598716, 0.006493816315762163, 0.006042963902381375, 0.005623413251903534, 0.005232991146814988, 0.00486967525165867, 0.0045315836376008545, 0.004216965034285858, 0.003924189758484568, 0.003651741272548408, 0.0033982083289425882, 0.0031622776601684067, 0.002942727176209307, 0.002738419634264385, 0.002548296747979369, 0.002371373705661676, 0.0022067340690846097, 0.0020535250264571646, 0.001910952974970458, 0.0017782794100389392, 0.001654817099943197, 0.0015399265260595066, 0.0014330125702369763, 0.0013335214321633369])

# High resolution diameters (particles)
dps_hires=2*np.array([0.3749471046662255, 0.3886825151193854, 0.4029210938807383, 0.41768127347891043, 0.43298216168002995, 0.4488435662236543, 0.4652860204648466, 0.48233080995559663, 0.4999999999999969, 0.5183164642188458, 0.5373039141606555, 0.5569869299973979, 0.5773909923447257, 0.5985425152478614, 0.6204688803758561, 0.6431984724684835, 0.6667607160816582, 0.6911861136789459, 0.7165062851184774, 0.7427540085863835, 0.7699632630297418, 0.7981692721439668, 0.8274085499715863, 0.8577189481714349, 0.8891397050194567, 0.9217114962045504, 0.9554764874852153, 0.9904783892751643, 1.0267625132285678, 1.064375830898181, 1.1033670345422895, 1.1437866001591923, 1.185686852830822, 1.229122034460093, 1.2741483739896673, 1.3208241601930402, 1.3692098171321745, 1.4193679823793712, 1.4713635881046347, 1.5252639451335066, 1.5811388300841833, 1.6390605756967227, 1.6991041644712728, 1.7613473257365437, 1.8258706362741812, 1.8927576246293076, 1.9620948792422603, 2.033972160541516, 2.108482517142903, 2.1857224063055365, 2.2657918188004005, 2.3487944083532373, 2.4348376258293065, 2.5240328583337264, 2.616495573407464, 2.712345468505654, 2.811706625951736, 2.914707673568027, 3.021481951190654, 3.1321676832844174, 3.246908157881046, 3.3658519120724804, 3.489152924299321, 3.616970813683363, 3.749471046662268, 3.886825151193868, 4.029210938807398, 4.176812734789119, 4.329821616800315, 4.488435662236559, 4.652860204648483, 4.823308099555984, 4.999999999999988, 5.183164642188477, 5.373039141606575, 5.569869299973999, 5.773909923447278, 5.985425152478636, 6.204688803758584, 6.431984724684859, 6.667607160816606, 6.911861136789483, 7.165062851184799, 7.4275400858638605, 7.699632630297445, 7.981692721439697, 8.274085499715891, 8.57718948171438, 8.891397050194598, 9.217114962045537, 9.554764874852188, 9.904783892751679, 10.267625132285717, 10.643758308981848, 11.033670345422935, 11.437866001591965, 11.856868528308263, 12.291220344600974, 12.74148373989672, 13.20824160193045, 13.692098171321794, 14.193679823793762, 14.713635881046399, 15.252639451335117, 15.81138830084189, 16.390605756967286, 16.99104164471279, 17.613473257365502, 18.258706362741883, 18.92757624629315, 19.620948792422677, 20.339721605415235, 21.08482517142911])

# Calculate ion electrical mobility diameters from mobilities
# -Assuming standard conditions
# -On mountain tops and in airplanes one should use airborne type NAIS, for this conversion to remain valid.
dpmil_lores = []
for i in range(0,len(mob_lores)):
    dpmil_lores.append(mob_to_diam(mob_lores[i]))
dpmil_lores=np.array(dpmil_lores)

dpmil_hires = []
for i in range(0,len(mob_hires)):
    dpmil_hires.append(mob_to_diam(mob_hires[i]))
dpmil_hires=np.array(dpmil_hires)

# Inverter resolution
if inv=='hires':
    mob = mob_hires
    dpmil = dpmil_hires
    dps = dps_hires
if inv=='lores':
    mob = mob_lores
    dpmil = dpmil_lores
    dps = dps_lores

# Calculate dlogdpmil and dlogmob
ldpmil = np.log10(dpmil)
ldpmil_limits = np.array([(ldpmil[i] + ldpmil[i+1])/2.0 for i in range(0,len(ldpmil)-1)])
head = ldpmil[0] - (ldpmil_limits[0] - ldpmil[0])
tail = ldpmil[-1] + (ldpmil[-1] - ldpmil_limits[-1])
ldpmil_limits = np.insert(ldpmil_limits,0,head)
ldpmil_limits = np.append(ldpmil_limits,tail)
dlogdpmil = np.diff(ldpmil_limits)

lmob = np.log10(mob)
lmob_limits = np.array([(lmob[i] + lmob[i+1])/2.0 for i in range(0,len(lmob)-1)])
head = lmob[0] + (lmob[0] - lmob_limits[0])
tail = lmob[-1] - (lmob_limits[-1] - lmob[-1])
lmob_limits = np.insert(lmob_limits,0,head)
lmob_limits = np.append(lmob_limits,tail)
dlogmob = np.abs(np.diff(lmob_limits))

# Calculate dlogdps
ldps = np.log10(dps)
ldps_limits = np.array([(ldps[i] + ldps[i+1])/2.0 for i in range(0,len(ldps)-1)])
head = ldps[0] - (ldps_limits[0] - ldps[0])
tail = ldps[-1] + (ldps[-1] - ldps_limits[-1])
ldps_limits = np.insert(ldps_limits,0,head)
ldps_limits = np.append(ldps_limits,tail)
dlogdps = np.diff(ldps_limits)

# Find out the filenames in the destination directory.
filenames = os.listdir(path)

records_read = 0

# Now run through the date range given and check for existence of files
for i in range(0,len(datetimes)):

    #different date formats encountered in NAIS data files.
    date = datetime.strftime(datetimes[i],datefmt)

    #Date format used in the saved sum-files and figures
    mydatefmt = datetime.strftime(datetimes[i],'%Y%m%d')  

    # figure out filenames
    records_file = ''
    ions_file = ''
    particles_file = ''
    for x in filenames:
        if (date in x):
            # extract the different type of files
            if (date+recordfmt==x):
                records_file = x
            if (date+ionfmt==x):
                ions_file = x
            if (date+particlefmt==x):
                particles_file = x

    # Check if data exists
    if (records_file!='') and (ions_file!=''):

        # Read records (=diagnostic parameter file) data and index it with opmode
        records = pd.read_table(path+records_file,sep=delimiter,comment='#',engine='python')
        records = records.set_index('opmode')

        # Find the column names for time, temperature and pressure (other things?)
        # (their column number and name change from instrument to instrument so I rely on user input here)
        if records_read == 0:
            list_of_records = list(records)
            time_name = easygui.choicebox("Choose time from parameters (use begin time) and press OK", title, list_of_records)
            if tp_data == 'yes':
                temperature_name = easygui.choicebox("Choose temperature from parameters and press OK", title, list_of_records)
                pressure_name = easygui.choicebox("Choose pressure from parameters and press OK", title, list_of_records)
            records_read = 1

        # Separate the ion records data and index it with time
        # Time will be in UTC after this
        ion_records = records.loc['ions'].set_index(time_name)
        ion_records.index = ion_records.index.map(pd.to_datetime)

        # Read ion data and index it with time
        ions = pd.read_table(path+ions_file,sep=delimiter,comment='#',engine='python').set_index(time_name) 
        ions.index = ions.index.map(pd.to_datetime)

        # reindex the records to match the ions
        ion_records = ion_records.reindex(index=ions.index)

        # Get negative ions and positive ions -> transform to numpy arrays because matrix operations ensue
        neg_ions = ions.iloc[:,3:3+len(dpmil)].interpolate().values
        pos_ions = ions.iloc[:,3+len(dpmil)+len(dpmil):3+len(dpmil)+len(dpmil)+len(dpmil)].interpolate().values

        # Convert ion data from dndlogmob -> dndlogdp
        neg_ions = neg_ions * dlogmob / dlogdpmil
        pos_ions = pos_ions * dlogmob / dlogdpmil

        # do the same steps for particles if they exist
        if particles_file!='':
            particle_records=records.loc['particles'].set_index(time_name)
            particle_records.index=particle_records.index.map(pd.to_datetime)
            particles = pd.read_table(path+particles_file,sep=delimiter,comment='#',engine='python').set_index(time_name)
            particles.index = particles.index.map(pd.to_datetime)
            particle_records = particle_records.reindex(index=particles.index)                
            neg_particles = particles.iloc[:,3:3+len(dps)].interpolate().values
            pos_particles = particles.iloc[:,3+len(dps)+len(dps):3+len(dps)+len(dps)+len(dps)].interpolate().values

        # Get temperature and pressure (1-D arrays)
        if tp_data=='yes':
            temperature_ions = (273.15 + ion_records[temperature_name].interpolate()).values.flatten()
            pressure_ions = (100.0 * ion_records[pressure_name].interpolate()).values.flatten()
            if particles_file!='':
                temperature_particles = (273.15 + particle_records[temperature_name].interpolate()).values.flatten()
                pressure_particles = (100.0 * particle_records[pressure_name].interpolate()).values.flatten()
        # If temperature and pressure were not recorded, use user defiend mean conditions.
        else:
            temperature_ions=temp_inst*np.ones(ion_records.index.shape)
            pressure_ions=pres_inst*np.ones(ion_records.index.shape)
            if particles_file != '':
                temperature_particles=temp_inst*np.ones(particle_records.index.shape)
                pressure_particles=pres_inst*np.ones(particle_records.index.shape)

        # STP correction (each time row needs to be multiplied by a unique value)
        if stp_corr == 'yes':
            stp_corr_ions = (pres_ref*temperature_ions)/(temp_ref*pressure_ions)
            neg_ions = (stp_corr_ions*neg_ions.T).T
            pos_ions = (stp_corr_ions*pos_ions.T).T
            if particles_file!='':
                stp_corr_particles = (pres_ref*temperature_particles)/(temp_ref*pressure_particles)
                neg_particles = (stp_corr_particles*neg_particles.T).T
                pos_particles = (stp_corr_particles*pos_particles.T).T
        
        # tubeloss correction
        if tubeloss_corr=='yes':
            throughput_ions = np.zeros(neg_ions.shape)
            for i in range(0,throughput_ions.shape[0]):
                throughput_ions[i,:] = tubeloss(dpmil,sf*1.667e-5,inlet_diameter,inlet_length,temperature_ions[i],pressure_ions[i])
            neg_ions = throughput_ions*neg_ions
            pos_ions = throughput_ions*pos_ions
            if particles_file!='':
                throughput_particles = np.zeros(neg_particles.shape)
                for i in range(0,throughput_particles.shape[0]):
                    throughput_particles[i,:] = tubeloss(dps,sf*1.667e-5,inlet_diameter,inlet_length,temperature_particles[i],pressure_particles[i])
                neg_particles = throughput_particles*neg_particles
                pos_particles = throughput_particles*pos_particles
            
        # Calculate the total number concentration
        total_neg_ions = np.nansum(neg_ions*dlogdpmil,axis=1)
        total_pos_ions = np.nansum(pos_ions*dlogdpmil,axis=1)
        if particles_file!='': 
            total_neg_particles = np.nansum(neg_particles*dlogdps,axis=1)  
            total_pos_particles = np.nansum(pos_particles*dlogdps,axis=1)

        # Get the time vector for ions and particles and convert it to matlab datenum taking into account the time zone
        time_ions = np.array([datetime2datenum(x) for x in ions.index]) + time_zone/24.
        if particles_file!='':
            time_particles = np.array([datetime2datenum(x) for x in particles.index]) + time_zone/24.

        # Construct the sum-matrices and save them to data folder
        negions=np.concatenate((np.insert(dpmil*1e-9,(0,0),0)[np.newaxis],np.concatenate((time_ions[np.newaxis].T,total_neg_ions[np.newaxis].T,neg_ions),axis=1)))
        posions=np.concatenate((np.insert(dpmil*1e-9,(0,0),0)[np.newaxis],np.concatenate((time_ions[np.newaxis].T,total_pos_ions[np.newaxis].T,pos_ions),axis=1)))
        np.savetxt(path+'sum/negions_' + mydatefmt + '.sum',negions)
        np.savetxt(path+'sum/posions_' + mydatefmt + '.sum',posions)
       
        if particles_file!='':
            negparticles=np.concatenate((np.insert(dps*1e-9,(0,0),0)[np.newaxis],np.concatenate((time_particles[np.newaxis].T,total_neg_particles[np.newaxis].T,neg_particles),axis=1)))
            posparticles=np.concatenate((np.insert(dps*1e-9,(0,0),0)[np.newaxis],np.concatenate((time_particles[np.newaxis].T,total_pos_particles[np.newaxis].T,pos_particles),axis=1)))
            np.savetxt(path+'sum/negparticles_' + mydatefmt + '.sum',negparticles)
            np.savetxt(path+'sum/posparticles_' + mydatefmt + '.sum',posparticles)

        # Plot the number size distributions
        if particles_file!='':
            fig,ax = plt.subplots(2,2,figsize=(10,7))
            ax = ax.flatten()
        else:
            fig,ax = plt.subplots(2,1,figsize=(10,5))
            ax = ax.flatten()
  
        (sc_negion,cb_negion) = plot_sumfile2(ax[0],time_ions,dpmil,neg_ions)
        (sc_posion,cb_posion) = plot_sumfile2(ax[1],time_ions,dpmil,pos_ions)
        sc_negion.set_clim((1,10000))
        sc_posion.set_clim((1,10000))
  
        ax[0].set_title('Ions ($-$)')
        ax[1].set_title('Ions ($+$)')

        for i in range(0,2):
            ax[i].grid('on',which='both',linestyle='--',color='k')
            ax[i].xaxis.set_major_locator(dts.MinuteLocator(interval=60))
            ax[i].xaxis.set_major_formatter(dts.DateFormatter('%H'))
            plt.setp(ax[i].get_xticklabels(), rotation=60)
            ax[i].set_ylabel('$d_p$, [m]')
            ax[i].set_xlabel('Time')

        cb_negion.set_label('$dN/d\log d_p$, [cm$^{-3}$]')                      
        cb_posion.set_label('$dN/d\log d_p$, [cm$^{-3}$]')

        if particles_file!='':
            (sc_negpar,cb_negpar) = plot_sumfile2(ax[2],time_particles,dps,neg_particles)
            (sc_pospar,cb_pospar) = plot_sumfile2(ax[3],time_particles,dps,pos_particles)
            sc_negpar.set_clim((1,100000))
            sc_pospar.set_clim((1,100000))
      
            ax[2].set_title('Particles ($-$)')
            ax[3].set_title('Particles ($+$)')
      
            fig.suptitle(mydatefmt)
      
            for i in range(2,4):
                ax[i].grid('on',which='both',linestyle='--',color='k')
                ax[i].xaxis.set_major_locator(dts.MinuteLocator(interval=60))
                ax[i].xaxis.set_major_formatter(dts.DateFormatter('%H'))
                plt.setp(ax[i].get_xticklabels(), rotation=60)
                ax[i].set_ylabel('$d_p$, [m]')
                ax[i].set_xlabel('Time')
  
            cb_negpar.set_label('$dN/d\log d_p$, [cm$^{-3}$]')
            cb_pospar.set_label('$dN/d\log d_p$, [cm$^{-3}$]')
 
        fig.suptitle(mydatefmt,fontsize=11)

        plt.tight_layout()
  
        plt.savefig(path+'figures/'+'nais_'+mydatefmt+'.png',bbox_inches='tight',dpi=300)
        print mydatefmt+": Success!"
        plt.close

    else:
        print "No data for", mydatefmt
