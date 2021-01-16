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
import yaml
import re
import sys
from dateutil.parser import parse
from tinydb import TinyDB, Query
from tinydb.operations import add

# Fixed diameter and mobility bins
dp_ion = np.array([7.949610066873187275e-01,9.181737924552214603e-01,1.060513600503926179e+00,1.224959679823698799e+00,1.414958699738506631e+00,1.634499249798819331e+00,1.888198514085806856e+00,2.181403433339687226e+00,2.520308747865528165e+00,2.912095102815642989e+00,3.365090891236600878e+00,3.888962384293289887e+00,4.494937535166431353e+00,5.196070414640996837e+00,6.007554438162747701e+00,6.947095098447752193e+00,8.035355151375323857e+00,9.296489193192451594e+00,1.075878902024538242e+01,1.245546773082500103e+01,1.442561898219513949e+01,1.671539984850161886e+01,1.937950186998520152e+01,2.248299804137784363e+01,2.610368545677439300e+01,3.033508982931992648e+01,3.529036394466827886e+01,4.110740875515996606e+01])
dp_par = np.array([7.498942093324539870e-01,8.659643233600640144e-01,9.999999999999980016e-01,1.154781984689456031e+00,1.333521432163321974e+00,1.539926526059490097e+00,1.778279410038920094e+00,2.053525026457140079e+00,2.371373705661659947e+00,2.738419634264360081e+00,3.162277660168379967e+00,3.651741272548380213e+00,4.216965034285819591e+00,4.869675251658620141e+00,5.623413251903479626e+00,6.493816315762099833e+00,7.498942093324560076e+00,8.659643233600640144e+00,1.000000000000000000e+01,1.154781984689457985e+01,1.333521432163323972e+01,1.539926526059490008e+01,1.778279410038922137e+01,2.053525026457139901e+01,2.371373705661660125e+01,2.738419634264360170e+01,3.162277660168379967e+01,3.651741272548380124e+01,4.216965034285819769e+01])
mob_ion = np.array([3.162277660168379937e-04,2.371373705661659990e-04,1.778279410038920258e-04,1.333521432163320159e-04,1.000000000000000048e-04,7.498942093324559917e-05,5.623413251903490022e-05,4.216965034285820205e-05,3.162277660168380208e-05,2.371373705661660125e-05,1.778279410038919852e-05,1.333521432163319990e-05,1.000000000000000082e-05,7.498942093324561442e-06,5.623413251903490361e-06,4.216965034285830030e-06,3.162277660168380038e-06,2.371373705661659871e-06,1.778279410038920148e-06,1.333521432163330027e-06,1.000000000000000167e-06,7.498942093324570124e-07,5.623413251903499890e-07,4.216965034285829924e-07,3.162277660168379721e-07,2.371373705661660136e-07,1.778279410038920042e-07,1.333521432163329868e-07])*1e4

# Some other values calculated from the fixed bins
mob_ion_geomeans=np.array([2.73841963e-04, 2.05352503e-04, 1.53992653e-04, 1.15478198e-04,
       8.65964323e-05, 6.49381632e-05, 4.86967525e-05, 3.65174127e-05,
       2.73841963e-05, 2.05352503e-05, 1.53992653e-05, 1.15478198e-05,
       8.65964323e-06, 6.49381632e-06, 4.86967525e-06, 3.65174127e-06,
       2.73841963e-06, 2.05352503e-06, 1.53992653e-06, 1.15478198e-06,
       8.65964323e-07, 6.49381632e-07, 4.86967525e-07, 3.65174127e-07,
       2.73841963e-07, 2.05352503e-07, 1.53992653e-07])*1e4
dp_par_geomeans=np.array([ 0.80584219,  0.93057204,  1.07460783,  1.24093776,  1.43301257,
        1.6548171 ,  1.91095297,  2.20673407,  2.54829675,  2.94272718,
        3.39820833,  3.92418976,  4.53158364,  5.23299115,  6.0429639 ,
        6.97830585,  8.05842188,  9.30572041, 10.74607828, 12.40937761,
       14.3301257 , 16.548171  , 19.10952975, 22.06734069, 25.48296748,
       29.42727176, 33.98208329, 39.24189758])
dlogmob_ion=np.array([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
       0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
       0.125])
dlogdp_ion=np.array([0.06257907, 0.06258521, 0.06259845, 0.06261376, 0.06263147,
       0.06265194, 0.06267563, 0.06270305, 0.06273478, 0.06277153,
       0.06281409, 0.06286343, 0.06292064, 0.06298703, 0.06306411,
       0.06315368, 0.06325786, 0.06337916, 0.06352054, 0.06368553,
       0.06387836, 0.06410408, 0.06436873, 0.06467961, 0.06504553,
       0.06547715, 0.06598741, 0.06626396])
dlogdp_par=np.array([0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
       0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
       0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
       0.0625, 0.0625, 0.0625, 0.0625, 0.0625])

# All possible names and naming formats encountered in current NAIS data files
ion_filename_formats = [
'%Y-%m-%d.ions.nds',
'%Y%m%d-block-ions.spectra']

particle_filename_formats = [
'%Y-%m-%d.particles.nds',
'%Y%m%d-block-particles.spectra']

diagnostic_filename_formats = [
'%Y-%m-%d.log',
'%Y%m%d-block.records']

possible_sampleflow_names = [
'pos_sampleflow.mean',
'neg_sampleflow.mean',
'pos_sampleflow',
'neg_sampleflow',
'sampleflow',
'Flowaer']

possible_temperature_names = [
'temperature.mean',
'temperature',
]

possible_pressure_names = [
'baro.mean',
'baro']


def visc(temp):
    """ Calculate viscosity of air """

    nyy_ref=18.203e-6
    S=110.4
    temp_ref=293.15
    nyy=nyy_ref*((temp_ref+S)/(temp+S))*((temp/temp_ref)**(3./2.))
    return nyy

def rlambda(temp,press):
    """ Calculate mean-free path """

    kaasuv=8.3143
    dm=3.7e-10
    avoc=6.022e23
    return kaasuv*temp/(np.sqrt(2.)*avoc*press*np.pi*dm*dm)

def cunn(Dp,temp,press):
    """ Calculate Cunningham correction """

    lambd = rlambda(temp,press)
    return 1. + 2.*lambd/Dp * (1.165 + 0.483 * np.exp(-0.997*Dp/lambd))

def diffuus(dpp,temp,press):
    """ Calculate diffusion coefficient """

    K=1.38e-23
    return (K*temp*cunn(dpp,temp,press))/(3.*np.pi*visc(temp)*dpp)

def tubeloss(dpp,pflow,plength,temp,press):
    """ Laminar diffusion losses in circular straight conduit """

    DPP,TEMP = np.meshgrid(dpp,temp)
    DPP,PRESS = np.meshgrid(dpp,press)
    DPP,PFLOW = np.meshgrid(dpp,pflow)

    rmuu = np.pi*diffuus(DPP,TEMP,PRESS)*plength/PFLOW;
    pene = np.zeros(rmuu.shape)

    cond1=rmuu<0.02
    cond2=rmuu>=0.02
    pene[cond1]=1. - 2.56*rmuu[cond1]**(2./3.) + 1.2*rmuu[cond1]+0.177*rmuu[cond1]**(4./3.)
    pene[cond2]=1. - 2.56*rmuu[cond2]**(2./3.) + 1.2*rmuu[cond2]+0.177*rmuu[cond2]**(4./3.)

    return pene

#    for i in range(0,len(dpp)):
#        if rmuu[i] < 0.02:
#            pene[i]=1. - 2.56*rmuu[i]**(2./3.) + 1.2*rmuu[i]+0.177*rmuu[i]**(4./3.)
#        else:
#            pene[i]=0.819*np.exp(-3.657*rmuu[i])+0.097 \
#                    *np.exp(-22.3*rmuu[i])+0.032*np.exp(-57.*rmuu[i])
#    return pene

#def x2dlogx(x):
#    """ Calculate log-differences for monotonically
#    increasing or decreasing vector x """
#
#    logx = np.log10(x)
#    logx_mid = (logx[1:] + logx[:-1])/2.0 
#    logx_mid_first_value = logx[0] + (logx[0] - logx_mid[0])
#    logx_mid_last_value  = logx[-1] - (logx_mid[-1] - logx[-1])
#    logx_mid = np.insert(logx_mid,0,logx_mid_first_value)
#    logx_mid = np.append(logx_mid,logx_mid_last_value)
#    dlogx = np.abs(np.diff(logx_mid))
#    return dlogx

def datetime2datenum(dt):
    """ Convert from python datetime to matlab datenum """

    mdn = dt + timedelta(days = 366)
    frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0,tzinfo=dt.tzinfo)).seconds \
           / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac

#def datenum2datetime(matlab_datenum):
#    """ Convert from matlab datenum to python datetime """
#
#    return datetime.fromordinal(int(matlab_datenum)) \
#    + timedelta(days=matlab_datenum%1) \
#    - timedelta(days = 366)

def plot_sumfile(handle,v,clim=(10,100000)):
    """ Plot UHEL's sum-formatted aerosol number-size distribution """
    
    time = v[1:,0] # This is datenum
    dp = v[0,2:]
    data = v[1:,2:]
    data[data<=0]=1e-30 # No holes in plots
    mesh_dp, mesh_time = np.meshgrid(dp,time)
    pcolorplot = handle.pcolormesh(mesh_time,mesh_dp,data,
                                   norm=colors.LogNorm(),
                                   linewidth=0,rasterized=True,cmap='jet')
    handle.set_yscale('log')
    pcolorplot.set_clim(clim)
    pcolorplot.set_edgecolor('face')
    handle.autoscale(tight='true')
    #handle.set_xlim((np.floor(time[0]),np.floor(time[0])+1)) # makes 24-h axis

    handle.grid('on',which='both',linestyle='--',color='w',lw=0.5)
    handle.xaxis.set_major_locator(dts.HourLocator(interval=1))
    handle.xaxis.set_major_formatter(dts.DateFormatter('%H'))
    #handle.set_xticks(np.floor(time[0])+np.arange(0,25)/24.0)
    #handle.set_xticklabels(["%2.2i" % x for x in np.arange(0,25)])
    plt.setp(handle.get_xticklabels(), rotation=80)
    handle.set_ylabel('Dp, [m]')
    handle.set_xlabel('UTC'+'%+d'%v[0,0]+', [h]')
    cbar = plt.colorbar(pcolorplot, ax=handle, 
                        ticks=LogLocator(subs=range(10)))
    cbar.set_label('dN/dlogDp, [cm-3]')
    return pcolorplot

def find_delimiter(fn):
    with open(fn) as f:
        line = f.readline()
        while line:
            if line[0]=='#':
                line = f.readline()
                continue
            else:
                l = line
                break
    result = re.search('(.)opmode',l)
    delimiter = result.group(1)
    return delimiter


def average_mob(y,h):

    data = np.nan*np.ones((y.shape[0],len(mob_ion)))

    for i in range(0,len(mob_ion_geomeans)):
        if i==0:
            y_block = y[:,h>mob_ion_geomeans[i]]
        elif i==(len(mob_ion_geomeans)-1):
            y_block = y[:,h<=mob_ion_geomeans[i]]
        else:
            y_block = y[:,((h>mob_ion_geomeans[i]) & (h<=mob_ion_geomeans[i-1]))]
        data[:,i] = np.nanmean(y_block,axis=1)
        
    return data

def average_dp(y,h):

    data = np.nan*np.ones((y.shape[0],len(dp_par)))

    for i in range(0,len(dp_par_geomeans)):
        if i==0:
            y_block = y[:,h<dp_par_geomeans[i]]
        elif i==(len(dp_par_geomeans)-1):
            y_block = y[:,h>=dp_par_geomeans[i]]
        else:
            y_block = y[:,((h<dp_par_geomeans[i]) & (h>=dp_par_geomeans[i-1]))]
        data[:,i] = np.nanmean(y_block,axis=1)

    return data


def nais_processor(config_file):
    """ Function that processes data from the NAIS 
    
    The function reads the raw NAIS data files from the load_path, 
    applies corrections (diffusion losses in the inlet line, conversion 
    to standard conditions (optional) and R. Wagner's ion mode calibration) to the 
    measured number concentrations and saves the data as a University of 
    Helsinki sum-formatted number-size distribution to the save_path.

    A measurement setup specific configuration file is needed in the 
    processing (see README file for an example).

    The function tries to process raw files in the given time range
    and maintains a database of the processed files.

    The sum file format:
        [0,0]  = UTC offset in hours
        [1:,0] = Time (MATLAB datenum) 
        [0,2:] = Geometric mean diameter of size-channel (m)
        [1:,1] = Integrated total number concentration (cm-3)
        [1:,2:] = Normalized number concentrations, dN/dlogDp (cm-3)

    Function arguments:
        Name of the configuration file (str)

    Example:
        nais_processor('/home/user/data/config.yml')

    """
    warnings.filterwarnings("ignore")

    # Find out today
    today_dt = datetime.today()
    today = today_dt.strftime('%Y%m%d')

    # Check that the config file exists
    if os.path.isfile(config_file)==False:
        raise Exception('"%s" does not exist' % config_file)
    else:
        with open(config_file,'r') as stream:
            try:
                config = yaml.load(stream)
                load_path = config['data_folder']
                save_path = config['processed_folder']
                start_date = config['start_date']
                database = config['database_file']
                if 'end_date' in config:
                    end_date = config['end_date']
                else:
                    end_date = today
                if 'inlet_length'in config:
                    pipelength = config['inlet_length']
                else:
                    pipelength=0.0
                if 'sealevel_correction' in config:
                    sealevel_correction = config['sealevel_correction']
                else:
                    sealevel_correction = False
            except Exception as error_msg:
                raise Exception("bad configuration file")

    # Initialize the database
    try:
      db = TinyDB(database)
      check = Query()
    except Exception as error_msg:
        raise Exception(error_msg)
   
    # Test if the configuration information is valid
    try:
        float(pipelength)
    except:
        raise Exception('"%s" must be a number' % pipelength)
    if not os.path.exists(load_path):
        raise Exception('"%s" does not exist' % load_path)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    try:
       start_dt = pd.to_datetime(start_date)
       end_dt = pd.to_datetime(end_date)
    except:
       raise Exception('bad start_date or end_date')

    model = 'NAIS'

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    # Convert load and save paths to absolute paths
    load_path = os.path.abspath(load_path) + '/'
    save_path = os.path.abspath(save_path) + '/'

    last_date = end_date
    first_date = start_date

    # Make a list of datetimes
    list_of_datetimes = pd.date_range(start=first_date, end=last_date)

    # list existing files:
    list_of_existing_ion_files = [os.path.split(x['ions'])[1] for x in db.search(check.ions.exists())]
    list_of_existing_particle_files = [os.path.split(x['particles'])[1] for x in db.search(check.particles.exists())]
    list_of_existing_diagnostic_files = [os.path.split(x['diagnostics'])[1] for x in db.search(check.diagnostics.exists())]

    #List all possible filenames in the date range for particles, ions and diagnostic files
    list_of_ion_files=[]
    list_of_ion_file_dates=[]
    for y in ion_filename_formats:
        list_of_ion_files = list_of_ion_files + [x.strftime(y) for x in list_of_datetimes if x.strftime(y) not in list_of_existing_ion_files]
        list_of_ion_file_dates = list_of_ion_file_dates + [x.strftime('%Y%m%d') for x in list_of_datetimes if x.strftime(y) not in list_of_existing_ion_files] 

    list_of_particle_files = []
    list_of_particle_file_dates = []
    for y in particle_filename_formats:
        list_of_particle_files = list_of_particle_files + [x.strftime(y) for x in list_of_datetimes if x.strftime(y) not in list_of_existing_particle_files]
        list_of_particle_file_dates = list_of_particle_file_dates + [x.strftime('%Y%m%d') for x in list_of_datetimes if x.strftime(y) not in list_of_existing_particle_files]

    list_of_diagnostic_files=[]
    list_of_diagnostic_file_dates=[]
    for y in diagnostic_filename_formats:
        list_of_diagnostic_files = list_of_diagnostic_files + [x.strftime(y) for x in list_of_datetimes if x.strftime(y) not in list_of_existing_diagnostic_files]
        list_of_diagnostic_file_dates = list_of_diagnostic_file_dates + [x.strftime('%Y%m%d') for x in list_of_datetimes if x.strftime(y) not in list_of_existing_diagnostic_files]

    # Initialize entries to the database with the timestamps
    list_of_dates = [x.strftime('%Y%m%d') for x in list_of_datetimes]
    for x in list_of_dates:
        if bool(db.search(check.timestamp==x)):
            continue
        else:
            db.insert({'timestamp':x,'ion_error':'','particle_error':''})

    # Descend into the raw data folder
    for root, dirs, files in os.walk(load_path):
 
      # Find matching files
      ion_findex_all,ion_findex_actual = np.intersect1d(files,list_of_ion_files,return_indices=True)[1:] 
      particle_findex_all,particle_findex_actual = np.intersect1d(files,list_of_particle_files,return_indices=True)[1:]
      diagnostic_findex_all,diagnostic_findex_actual = np.intersect1d(files,list_of_diagnostic_files,return_indices=True)[1:] 

      # Put the files found into the database
      for i in range(0,len(ion_findex_all)):
        full_name = os.path.join(root, files[ion_findex_all[i]])
        ion_datestr = list_of_ion_file_dates[ion_findex_actual[i]]

        # If there is raw data file: don't replace it.
        if bool(db.search((check.timestamp==ion_datestr) & ~check.ions.exists())):
            db.update({'ions':full_name},check.timestamp==ion_datestr)

      for i in range(0,len(particle_findex_all)):
        full_name = os.path.join(root, files[particle_findex_all[i]])
        particle_datestr = list_of_particle_file_dates[particle_findex_actual[i]]
        db.update({'particles':full_name},check.timestamp==particle_datestr)

        if bool(db.search((check.timestamp==particle_datestr) & ~check.particles.exists())):
            db.update({'particles':full_name},check.timestamp==particle_datestr)

      for i in range(0,len(diagnostic_findex_all)):
        full_name = os.path.join(root, files[diagnostic_findex_all[i]])
        diagnostic_datestr = list_of_diagnostic_file_dates[diagnostic_findex_actual[i]]
        db.update({'diagnostics':full_name},check.timestamp==diagnostic_datestr)

        if bool(db.search((check.timestamp==diagnostic_datestr) & ~check.diagnostics.exists())):
            db.update({'diagnostics':full_name},check.timestamp==diagnostic_datestr)

    # From the database find the last days with processed data
    processed_ion_days = db.search( (check.processed_neg_ion_file.exists() &
                                     check.processed_pos_ion_file.exists()) )
    if len(processed_ion_days)!=0:
      last_ion_day=np.max([datetime.strptime(x['timestamp'],'%Y%m%d') for x in processed_ion_days]).strftime('%Y%m%d')
    else:
      last_ion_day='empty'
    
    processed_particle_days=db.search( (check.processed_neg_particle_file.exists() &
                                        check.processed_pos_particle_file.exists()) )
    if len(processed_particle_days)!=0:
      last_particle_day=np.max([datetime.strptime(x['timestamp'],'%Y%m%d') for x in processed_particle_days]).strftime('%Y%m%d')
    else:
      last_particle_day='empty'
 


    # Define standard conditions
    temp_ref = 273.15 # K
    pres_ref = 101325.0 # Pa
    
    # Try to process unprocessed files or there is some error in processing
    # only do the time range in the config file
    for x in iter(db.search( ((check.timestamp==last_ion_day) &
                             (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str)) |
                             (check.diagnostics.exists() & check.ions.exists() &
                             ~check.processed_neg_ion_file.exists() &
                             ~check.processed_pos_ion_file.exists() &
                             (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str))                             
                           )):

        try:

            # search for delimiter
            delimiter = find_delimiter(x['ions'])

            print('processing %s' % os.path.split(x['ions'])[1])

            # Read the ion data
            ions = pd.read_table(x['ions'],
                                 sep=delimiter,
                                 comment='#',
                                 error_bad_lines=False,
                                 engine='python',
                                 header=None)

            # Read the diagnostic data
            records = pd.read_table(x['diagnostics'],
                                    sep=delimiter,
                                    comment='#',
                                    engine='python',
                                    error_bad_lines = False,
                                    header=None)

            # Remove rows with too few fields
            ions = ions[ions.count(1)>len(ions.columns)/3]
            records = records[records.count(1)>len(records.columns)/3]

            # Remove duplicate rows based on first column (time)
            ions = ions.drop_duplicates(subset=0)
            records = records.drop_duplicates(subset=0)

            # Set the first row as the header and remove it from the actual data
            ions.columns = ions.iloc[0,:]
            ions = ions.drop(0)
            records.columns = records.iloc[0,:]
            records = records.drop(0)

            # Set the time row as index
            ions = ions.set_index(ions.columns[0])
            ions.index = [parse(y) for y in ions.index]

            # Get the list of column headers
            ion_columns = ions.columns

            # figure out inverter resolution
            inverter_reso = int((len(ion_columns)-2)/4)
     
            # get the mobilities after the particular inversion
            mob_ion_inv = np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0]) for y in ion_columns[2:2+inverter_reso]])

            # get the number densities
            neg_ions = ions.iloc[:,2:2+inverter_reso].astype(float).interpolate().values
            pos_ions = ions.iloc[:,2+2*inverter_reso:2+3*inverter_reso].astype(float).interpolate().values

            # average to the fixed ion mobs
            neg_ions = average_mob(neg_ions,mob_ion_inv)
            pos_ions = average_mob(pos_ions,mob_ion_inv)

            neg_ions = neg_ions * dlogmob_ion / dlogdp_ion
            pos_ions = pos_ions * dlogmob_ion / dlogdp_ion

            # Index records by the operation mode
            records = records.set_index('opmode')

            # Then extract the ion records
            ion_records = records.loc['ions'].set_index(records.columns[0])
            ion_records.index = [parse(y) for y in ion_records.index]            

            # Match ion records to ions according to time
            ion_records = ion_records.reindex(index=ions.index,method='nearest')

            # Get the list of diagnostic parameters
            diag_params = list(records)

            temperature_name = ''
            for temp_name in possible_temperature_names:
                if temp_name in diag_params:
                    temperature_name = temp_name
                    break

            # No temperature data, use standard conditions
            if temperature_name=='':
                temp_ions = temp_ref*np.ones(neg_ions.shape[0])
            # Read the temperature data
            else:
                temp_ions = 273.15 + ion_records[temperature_name].astype(float).interpolate().values.flatten()

            pressure_name = ''
            for pres_name in possible_pressure_names:
                if pres_name in diag_params:
                    pressure_name = pres_name
                    break

            # Extract pressure data from the diagnostic files
            if pressure_name=='':
                pres_ions = pres_ref*np.ones(neg_ions.shape[0])
            else:
                pres_ions = 100.0 * ion_records[pressure_name].astype(float).interpolate().values.flatten()

            # Correct the number concentrations to standard conditions (optional)
            if sealevel_correction:
                stp_corr_ions = (pres_ref*temp_ions)/(temp_ref*pres_ions)
                neg_ions = (stp_corr_ions*neg_ions.T).T
                pos_ions = (stp_corr_ions*pos_ions.T).T

            sampleflow_name = []
            for flow_name in possible_sampleflow_names:
                if flow_name in diag_params:
                    sampleflow_name.append(flow_name)

            # Extract sample flow rate data from the diagnostic files
            try:
                flow_ions = ion_records[sampleflow_name].astype(float).sum(axis=1).interpolate().values.flatten()
            except ValueError:
                flow_ions = ion_records[sampleflow_name].astype(float).interpolate().values.flatten()

            # Test if the sampleflow is in cm3/s (old models) or l/min and possibly convert to l/min
            if flow_ions[0]>100:
                flow_ions = (flow_ions/1000.0) * 60.0
            else:
                pass

            # Diffusion loss correction
            throughput_ions = tubeloss(dp_ion*1e-9,flow_ions*1.667e-5,pipelength,temp_ions,pres_ions)

            neg_ions = neg_ions / throughput_ions
            pos_ions = pos_ions / throughput_ions

            # Robert Wagner's calibration (only ions)
            roberts_corr = 0.713*dp_ion**0.120
            neg_ions = neg_ions / roberts_corr
            pos_ions = pos_ions / roberts_corr

            # If all data is NaNs then skip
            if (np.all(np.isnan(neg_ions)) | np.all(np.isnan(pos_ions))):
                db.update({'ion_error': 'All ion data are NaNs'}, check.timestamp==x['timestamp'])
                continue

            # Integrate total ion number concentrations
            total_neg_ions = np.nansum(neg_ions*dlogdp_ion,axis=1)[np.newaxis].T
            total_pos_ions = np.nansum(pos_ions*dlogdp_ion,axis=1)[np.newaxis].T

            # Get the utc offset
            if ions.index[0].utcoffset()==None:
                utc_offset_ions = 0
            else:
                utc_offset_ions = ions.index[0].utcoffset().total_seconds()/(60.0*60.0)

            time_ions = np.array([datetime2datenum(x) for x in ions.index])[np.newaxis].T

            # Construct the header
            ion_header = np.insert(dp_ion*1e-9,0,(utc_offset_ions,0))[np.newaxis]

            # Construct the sum-files
            negions = np.concatenate((ion_header,np.concatenate((time_ions,total_neg_ions,neg_ions),axis=1)))
            posions = np.concatenate((ion_header,np.concatenate((time_ions,total_pos_ions,pos_ions),axis=1)))

            # Save the sum matrices using the standard names
            np.savetxt(save_path+model+'n'+x['timestamp']+'nds.sum',negions)
            np.savetxt(save_path+model+'p'+x['timestamp']+'nds.sum',posions)
            
            # Update the database
            db.update({'processed_neg_ion_file': save_path+model+'n'+x['timestamp']+'nds.sum',
                       'processed_pos_ion_file': save_path+model+'p'+x['timestamp']+'nds.sum'}, check.timestamp==x['timestamp'])
            db.update({'ion_error':''}, check.timestamp==x['timestamp']) 
           
        except Exception as error_msg:
            db.update({'ion_error': str(error_msg)}, check.timestamp==x['timestamp'])
            continue
 
    # Iterate through the unprocessed particle data excluding today
    for x in iter(db.search( ((check.timestamp==last_particle_day) & 
                             (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str)) |
                             (check.diagnostics.exists() & check.particles.exists() &
                             ~check.processed_neg_particle_file.exists() &
                             ~check.processed_pos_particle_file.exists() &
                             (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str))                             
                           )):

        try:

            # find delimiter
            delimiter=find_delimiter(x['particles'])

            print('processing %s' % os.path.split(x['particles'])[1])
 
            # Read the particle data
            particles = pd.read_table(x['particles'],
                                      header=None,
                                      sep=delimiter,
                                      comment='#',
                                      engine='python',
                                      error_bad_lines=False)

            # Read the diagnostic data
            records = pd.read_table(x['diagnostics'],
                                    header=None,
                                    sep=delimiter,
                                    comment='#',
                                    engine='python',
                                    error_bad_lines=False)

            # Remove rows with too few fields 
            particles = particles[particles.count(1)>len(particles.columns)/3]
            records = records[records.count(1)>len(records.columns)/3]

            # Remove duplicate rows based on first row
            particles = particles.drop_duplicates(subset=0)
            records = records.drop_duplicates(subset=0)

            # Set the first row as the header and remove it from the actual data
            particles.columns = particles.iloc[0,:]
            particles = particles.drop(0)
            records.columns = records.iloc[0,:]
            records = records.drop(0)

            particles = particles.set_index(particles.columns[0])
            particles.index = [parse(y) for y in particles.index]

            # Get the list of diagnostic parameters
            particle_columns = particles.columns
            inverter_reso = int((len(particle_columns)-2)/4)
     
            dp_par_inv = 2.0*np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",t)[0]) for t in particle_columns[2:2+inverter_reso]])

            neg_particles = particles.iloc[:,2:2+inverter_reso].astype(float).interpolate().values
            pos_particles = particles.iloc[:,2+2*inverter_reso:2+3*inverter_reso].astype(float).interpolate().values

            neg_particles = average_dp(neg_particles,dp_par_inv)
            pos_particles = average_dp(pos_particles,dp_par_inv)

            # Index records by the operation mode
            records = records.set_index('opmode')

            # And extract the particle records
            particle_records = records.loc['particles'].set_index(records.columns[0])
            particle_records.index = [parse(y) for y in particle_records.index]

            # Match the times in records to those in spectra
            particle_records = particle_records.reindex(index=particles.index, method='nearest')

            # List all variable names in the diagnostic file
            diag_params = list(records)

            temperature_name=''
            for temp_name in possible_temperature_names:
                if temp_name in diag_params:
                    temperature_name = temp_name
                    break

            # No temperature data, use standard conditions
            if temperature_name=='':
                temp_particles = temp_ref*np.ones(neg_particles.shape[0])
            # Read the temperature data
            else:
                temp_particles = 273.15 + particle_records[temperature_name].astype(float).interpolate().values.flatten()


            pressure_name = ''
            for pres_name in possible_pressure_names:
                if pres_name in diag_params:
                    pressure_name = pres_name
                    break           

            # Extract pressure data from the diagnostic files
            if pressure_name=='':
                pres_particles = pres_ref*np.ones(neg_particles.shape[0])
            else:
                pres_particles = 100.0 * particle_records[pressure_name].astype(float).interpolate().values.flatten()


            # Correct the number concentrations to standard conditions (optional)
            if sealevel_correction:
                stp_corr_particles = (pres_ref*temp_particles)/(temp_ref*pres_particles)
                neg_particles = (stp_corr_particles*neg_particles.T).T
                pos_particles = (stp_corr_particles*pos_particles.T).T

            sampleflow_name = []
            for flow_name in possible_sampleflow_names:
                if flow_name in diag_params:
                    sampleflow_name.append(flow_name)

            # Extract sample flow rate data from the diagnostic files
            try:
                flow_particles = particle_records[sampleflow_name].astype(float).sum(axis=1).interpolate().values.flatten()
            except ValueError:
                flow_particles = particle_records[sampleflow_name].astype(float).interpolate().values.flatten()

            # Test if the sampleflow is in cm3/s (old models) or l/min and possibly convert to l/min
            if flow_particles[0]>100:
                flow_particles = (flow_particles/1000.0) * 60.0
            else:
                pass

            throughput_particles = tubeloss(dp_par*1e-9,flow_particles*1.667e-5,pipelength,temp_particles,pres_particles)

            neg_particles = neg_particles / throughput_particles
            pos_particles = pos_particles / throughput_particles

            # If all data is NaNs then skip
            if (np.all(np.isnan(neg_particles)) | np.all(np.isnan(pos_particles))):
                db.update({'particle_error': 'All particle data are NaNs'}, check.timestamp==x['timestamp'])
                continue

            # Integrate total number concentrations
            total_neg_particles = np.nansum(neg_particles*dlogdp_par,axis=1)[np.newaxis].T
            total_pos_particles = np.nansum(pos_particles*dlogdp_par,axis=1)[np.newaxis].T

            # Get utc offset in hours
            if particles.index[0].utcoffset()==None:
                utc_offset_particles = 0
            else:
                utc_offset_particles = particles.index[0].utcoffset().total_seconds()/(60.0*60.0)

            # Define time vectors
            time_particles = np.array([datetime2datenum(x) for x in particles.index])[np.newaxis].T
            
            # Construct the header
            particle_header = np.insert(dp_par*1e-9,0,(utc_offset_particles,0))[np.newaxis]

            # Construct the sum-files
            negparticles = np.concatenate((particle_header,np.concatenate\
                                                ((time_particles,total_neg_particles,neg_particles),axis=1)))
            posparticles = np.concatenate((particle_header,np.concatenate\
                                                ((time_particles,total_pos_particles,pos_particles),axis=1)))

            # Save the sum matrices using the standard names
            np.savetxt(save_path+model+'n'+x['timestamp']+'np.sum',negparticles)
            np.savetxt(save_path+model+'p'+x['timestamp']+'np.sum',posparticles)

            # Update the database
            db.update({'processed_neg_particle_file': save_path+model+'n'+x['timestamp']+'np.sum',
                       'processed_pos_particle_file': save_path+model+'p'+x['timestamp']+'np.sum'}, check.timestamp==x['timestamp'])
            db.update({'particle_error':''},check.timestamp==x['timestamp'])

        except Exception as error_msg:
            db.update({'particle_error': str(error_msg)},check.timestamp==x['timestamp'])
            continue


def nais_plotter(config_file):
    """ Function to plot the processed NAIS data """

    warnings.filterwarnings("ignore")

    today_dt = datetime.today()
    today = today_dt.strftime('%Y%m%d')
    # Check that config file exists
    if os.path.isfile(config_file)==False:
        raise Exception('"%s" does not exist' % config_file)
    else:
      # Read from the configuration file
      with open(config_file,'r') as stream:
          try:
              config = yaml.load(stream)
              location = config['location']
              fig_path = config['figure_folder']
              database = config['database_file']
              start_date = config['start_date']
              if 'end_date' in config:
                    end_date = config['end_date']
              else:
                    end_date = today
          except Exception as error_msg:
              raise Exception(error_msg)

    try:
       start_dt = pd.to_datetime(start_date)
       end_dt = pd.to_datetime(end_date)
    except:
       raise Exception('bad start_date or end_date')

    # Test if figure path exists
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    # Test if the database file does not exist
    if os.path.isfile(database)==False:
        raise Exception('"%s" does not exist' % database)
    # Otherwise initialize the database
    else:
        db = TinyDB(database)
        check = Query()

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    fig_path = os.path.abspath(fig_path) + '/'

    plt.style.use('dark_background')

    fontsize = 14
    plt.rcParams.update({'font.size': fontsize,
                         'axes.titlesize': fontsize,
                         'axes.labelsize': fontsize,
                         'xtick.labelsize': fontsize,
                         'ytick.labelsize': fontsize,
                         'figure.titlesize': fontsize,
                         'legend.fontsize': fontsize})

    # Find last day with figures
    processed_ion_days = db.search( check.ion_figure.exists() )
    if len(processed_ion_days)!=0:
      last_ion_day=np.max([datetime.strptime(x['timestamp'],'%Y%m%d') for x in processed_ion_days]).strftime('%Y%m%d')
    else:
      last_ion_day='empty'

    processed_particle_days=db.search( check.particle_figure.exists() )
    if len(processed_particle_days)!=0:
      last_particle_day=np.max([datetime.strptime(x['timestamp'],'%Y%m%d') for x in processed_particle_days]).strftime('%Y%m%d')
    else:
      last_particle_day='empty'

    # Iterate through the unprocessed ion data excluding today
    for x in iter(db.search(    (~check.ion_figure.exists() &
                                 check.processed_neg_ion_file.exists() &
                                 check.processed_pos_ion_file.exists() & 
                                 (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str)) |
                                 ((check.timestamp==last_ion_day) & 
                                 (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str) ))):
    
           
        try:

            print('Plotting %s and %s' % (os.path.split(x['processed_neg_ion_file'])[1],
                                          os.path.split(x['processed_pos_ion_file'])[1]))

            pos_ion_data = np.loadtxt(x['processed_pos_ion_file'])
            neg_ion_data = np.loadtxt(x['processed_neg_ion_file'])

            model = re.findall('(.*)p[0-9]',os.path.split(x['processed_pos_ion_file'])[1])[0]

            fig,ax = plt.subplots(2,1,figsize=(7,7.5),dpi=100)
            ax = ax.flatten()

            plot_sumfile(ax[0],pos_ion_data,clim=(10,10000))
            plot_sumfile(ax[1],neg_ion_data,clim=(10,10000))
                       
            ax[0].set_title('Positive ions')
            ax[1].set_title('Negative ions')

            plt.tight_layout(rect=[0, 0.0, 1, 0.96])

            fig.suptitle(model + ' ' + x['timestamp'] + '\n' + location, y=1.0)

            plt.savefig(fig_path+model+'_ions_'+ x['timestamp'] +'.png',dpi=100)

            db.update({'ion_figure': fig_path+model+'_ions_'+ x['timestamp'] +'.png'}, check.timestamp==x['timestamp'])

            plt.close()

            db.update({'ion_fig_error':''},check.timestamp==x['timestamp'])

        except Exception as error_msg:
            db.update({'ion_fig_error': str(error_msg)}, check.timestamp==x['timestamp'])
            continue


    for x in iter(db.search(    (~check.particle_figure.exists() &
                                 check.processed_neg_particle_file.exists() &
                                 check.processed_pos_particle_file.exists() & 
                                 (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str)) |
                                 ((check.timestamp==last_particle_day) & 
                                 (check.timestamp>=start_date_str) & (check.timestamp<=end_date_str)) )):
 
        try:

            print('Plotting %s and %s' % (os.path.split(x['processed_neg_particle_file'])[1],
                                          os.path.split(x['processed_pos_particle_file'])[1]))

            pos_particle_data = np.loadtxt(x['processed_pos_particle_file'])
            neg_particle_data = np.loadtxt(x['processed_neg_particle_file'])

            model = re.findall('(.*)p[0-9]',os.path.split(x['processed_pos_particle_file'])[1])[0]

            fig,ax = plt.subplots(2,1,figsize=(7,7.5),dpi=100)
            ax = ax.flatten()

            plot_sumfile(ax[0],pos_particle_data,clim=(10,100000))
            plot_sumfile(ax[1],neg_particle_data,clim=(10,100000))
                       
            ax[0].set_title('Particles (positive polarity)')
            ax[1].set_title('Particles (negative polarity)')

            plt.tight_layout(rect=[0, 0.0, 1, 0.96])

            fig.suptitle(model + ' ' + x['timestamp'] + '\n' + location, y=1.0)

            plt.savefig(fig_path+model+'_particles_'+ x['timestamp'] +'.png',dpi=100)

            db.update({'particle_figure': fig_path+model+'_particles_'+ x['timestamp'] +'.png'}, check.timestamp==x['timestamp'])

            plt.close()

            db.update({'particle_fig_error':''},check.timestamp==x['timestamp'])

        except Exception as error_msg:
            db.update({'particle_fig_error': str(error_msg)}, check.timestamp==x['timestamp'])
            continue



