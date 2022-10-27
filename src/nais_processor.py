import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from datetime import date, datetime, timedelta
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
import time
import json
import aerosol_functions as af

__pdoc__ = {
    'tubeloss': False,
    'average_mob': False,
    'average_dp': False,
    'find_diagnostic_names': False,
    'process_data': False,
    'correct_data': False,
    'data_cleaner': False,
    'corona_limit': False,
    'corona_ion_cleaner': False
}

# The final geometric mean diameters of diameter and mobility bins
dp_ion = np.array([7.86360416e-10, 9.08232168e-10, 1.04902018e-09, 1.21167006e-09,
1.39958930e-09, 1.61672083e-09, 1.86762862e-09, 2.15759741e-09,
2.49274932e-09, 2.88018000e-09, 3.32811839e-09, 3.84611427e-09,
4.44525917e-09, 5.13844742e-09, 5.94068566e-09, 6.86946146e-09,
7.94518431e-09, 9.19171623e-09, 1.06370142e-08, 1.23139134e-08,
1.42610904e-08, 1.65242568e-08, 1.91576555e-08, 2.22259544e-08,
2.58066722e-08, 2.99933244e-08, 3.48995548e-08, 4.06646353e-08])*1e9

dp_par = np.array([7.498942093324539870e-01,8.659643233600640144e-01,
9.999999999999980016e-01,1.154781984689456031e+00,1.333521432163321974e+00,
1.539926526059490097e+00,1.778279410038920094e+00,2.053525026457140079e+00,
2.371373705661659947e+00,2.738419634264360081e+00,3.162277660168379967e+00,
3.651741272548380213e+00,4.216965034285819591e+00,4.869675251658620141e+00,
5.623413251903479626e+00,6.493816315762099833e+00,7.498942093324560076e+00,
8.659643233600640144e+00,1.000000000000000000e+01,1.154781984689457985e+01,
1.333521432163323972e+01,1.539926526059490008e+01,1.778279410038922137e+01,
2.053525026457139901e+01,2.371373705661660125e+01,2.738419634264360170e+01,
3.162277660168379967e+01,3.651741272548380124e+01,4.216965034285819769e+01])

mob_ion = np.array([3.162277660168379937e-04,2.371373705661659990e-04,
1.778279410038920258e-04,1.333521432163320159e-04,1.000000000000000048e-04,
7.498942093324559917e-05,5.623413251903490022e-05,4.216965034285820205e-05,
3.162277660168380208e-05,2.371373705661660125e-05,1.778279410038919852e-05,
1.333521432163319990e-05,1.000000000000000082e-05,7.498942093324561442e-06,
5.623413251903490361e-06,4.216965034285830030e-06,3.162277660168380038e-06,
2.371373705661659871e-06,1.778279410038920148e-06,1.333521432163330027e-06,
1.000000000000000167e-06,7.498942093324570124e-07,5.623413251903499890e-07,
4.216965034285829924e-07,3.162277660168379721e-07,2.371373705661660136e-07,
1.778279410038920042e-07,1.333521432163329868e-07])*1e4

mob_ion_geomeans=np.array([2.73841963e-04, 2.05352503e-04, 1.53992653e-04, 1.15478198e-04,
8.65964323e-05, 6.49381632e-05, 4.86967525e-05, 3.65174127e-05,
2.73841963e-05, 2.05352503e-05, 1.53992653e-05, 1.15478198e-05,
8.65964323e-06, 6.49381632e-06, 4.86967525e-06, 3.65174127e-06,
2.73841963e-06, 2.05352503e-06, 1.53992653e-06, 1.15478198e-06,
8.65964323e-07, 6.49381632e-07, 4.86967525e-07, 3.65174127e-07,
2.73841963e-07, 2.05352503e-07, 1.53992653e-07])*1e4

dp_par_geomeans=np.array([0.80584219,  0.93057204,  1.07460783,  1.24093776,  1.43301257,
1.6548171 ,  1.91095297,  2.20673407,  2.54829675,  2.94272718,
3.39820833,  3.92418976,  4.53158364,  5.23299115,  6.0429639 ,
6.97830585,  8.05842188,  9.30572041, 10.74607828, 12.40937761,
14.3301257 , 16.548171  , 19.10952975, 22.06734069, 25.48296748,
29.42727176, 33.98208329, 39.24189758])

dlogmob_ion=np.array([0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
0.125])

dlogdp_ion = np.array([0.06257524, 0.0625811 , 0.06259375, 0.06260838, 0.06262533,
0.06264495, 0.06266769, 0.06269404, 0.06272461, 0.06276008,
0.06280128, 0.06284916, 0.06290487, 0.06296974, 0.06304539,
0.0631337 , 0.06323696, 0.06335788, 0.06349974, 0.0636665 ,
0.06386292, 0.06409481, 0.06436924, 0.06469482, 0.06508209,
0.06554394, 0.06609614, 0.06639699])

dlogdp_par=np.array([0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
0.0625, 0.0625, 0.0625, 0.0625, 0.0625])

filename_formats = [
["%Y-%m-%d.ions.nds","%Y-%m-%d.particles.nds","%Y-%m-%d.log"],
["%Y%m%d-block-ions.spectra","%Y%m%d-block-particles.spectra","%Y%m%d-block.records"],
["%Y%m%d-block-ions.spectra","%Y%m%d-block-particles.spectra","%Y%m%d-block.diagnostics"]]

possible_sampleflow_names1 = [
"sampleflow",
"Flowaer"]

possible_sampleflow_names2 = [
"pos_sampleflow.mean",
"neg_sampleflow.mean",
"pos_sampleflow",
"neg_sampleflow"
]

possible_temperature_names = [
"temperature.mean",
"temperature",
"temp"]

possible_pressure_names = [
"baro.mean",
"baro"]

# Define standard conditions
temp_ref = 273.15 # K, 0C
pres_ref = 101325.0 # Pa, 1atm

def make_config():
    """
    Make a configuration file for processing NAIS data

    Running `make_config()` asks information about the
    measurement and the data, then writes the configuration
    file to the specified location.

    """

    # Collect data from the user
    print()
    print("Enter absolute path to configuration file.")
    print("For example: /home/user/campaign.yml")
    while True:
        config_file = input("> ")
        if len(config_file)>0:
            if not config_file.lower().endswith(".yml"):
                config_file = config_file + ".yml"
            break
        else:
            continue

    print()
    print("Enter absolute path(s) to raw data folder(s). Separate multiple paths with comma.")
    print("For example: /home/user/data/2021,/home/user/data/2022")
    while True:
        user_input = input("> ")
        if len(user_input)>0:
            load_path=user_input.split(",")
            break
        else:
            continue

    print()
    print("Enter absolute path to processed data folder.")
    print("For example: /home/user/campaign")
    while True:
        save_path = input("> ")
        if len(save_path)>0:
            break
        else:
            continue

    print()
    print("Enter start of measurement (YYYY-MM-DD)")
    print("Leave empty if you want to start from the earliest date")
    while True:
        start_date = input("> ")
        if len(start_date)==0:
            break
        try:
            start_dt = pd.to_datetime(start_date)
            break
        except:
            continue

    print()
    print("Enter end of measurement (YYYY-MM-DD)")
    print("If empty processor assumes current day")
    while True:
        end_date = input("> ")
        if len(end_date)==0:
            break
        try:
            end_dt = pd.to_datetime(end_date)
            break
        except:
            continue

    print()
    print("Enter absolute path to database file")
    print("For example: /home/user/campaign.json")
    while True:
        database_file = input("> ")
        if len(database_file)>0:
            if not database_file.lower().endswith(".json"):
                database_file = database_file + ".json"
            break
        else:
            continue

    print()
    print("Reprocess on re-run (True/False)")
    while True:
        allow_reprocessing = input("> ")
        if ((allow_reprocessing=='True') or (allow_reprocessing=='False')):
            if (allow_reprocessing=='True'):
                allow_reprocessing=True
            else:
                allow_reprocessing=False
            break
        else:
            continue

    print()
    print("Enter measurement location")
    print("For example: Helsinki, Kumpula")
    location = input("> ")

    ############## CLEANINIG
    print()
    print("Apply data cleaning procedures (True/False)")
    print("Attempt to remove corona ions and electrometer noise from data")
    while True:
        apply_cleaning = input("> ")
        if ((apply_cleaning=='True') or (apply_cleaning=='False')):
            if (apply_cleaning=='True'):
                apply_cleaning=True
            else:
                apply_cleaning=False
                remove_corona_ions=False
                remove_noisy_electrometers=False
                reclean=False
            break
        else:
            continue

    if apply_cleaning:
        print()
        print("Remove corona charger ions from particle data (True/False)")
        while True:
            remove_corona_ions = input("> ")
            if ((remove_corona_ions=='True') or (remove_corona_ions=='False')):
                if remove_corona_ions=='True':
                    remove_corona_ions=True
                else:
                    remove_corona_ions=False
                break
            else:
                continue

        print()
        print("Remove noisy electrometer data (True/False)")
        while True:
            remove_noisy_electrometers = input("> ")
            if ((remove_noisy_electrometers=='True') or (remove_noisy_electrometers=='False')):
                if remove_noisy_electrometers=='True':
                    remove_noisy_electrometers=True
                else:
                    remove_noisy_electrometers=False
                break
            else:
                continue

    ################### CORRECTIONS
    print()
    print("Apply corrections (True/False)")
    print("Requires a NAIS with temperature and pressure sensors.")
    while True:
        apply_corrections = input("> ")
        if ((apply_corrections=='True') or (apply_corrections=='False')):
            if apply_corrections=='True':
                apply_corrections=True
            else:
                apply_corrections=False
                sealevel_correction=False
                recorrect=False
                inlet_length = 0.0
            break
        else:
            continue

    if apply_corrections:
        print()
        print("Length of the inlet in meters")
        while True:
            inlet_length = input("> ")
            try:
                inlet_length = float(inlet_length)
                break
            except:
                continue

        print()
        print("Correct concentrations to sealevel conditions (True/False)")
        while True:
            sealevel_correction = input("> ")
            if sealevel_correction=='True':
                sealevel_correction=True
                break
            elif sealevel_correction=='False':
                sealevel_correction=False
                break
            else:
                continue

    print()
    print("Configuration saved to: %s"%config_file)
    print()

    # Make a dictionary out of the user input
    config_info = {
        "data_folder": load_path,
        "processed_folder": save_path,
        "start_date": start_date,
        "end_date": end_date,
        "database_file": database_file,
        "location": location,
        "inlet_length": inlet_length,
        "apply_corrections":apply_corrections,
        "sealevel_correction": sealevel_correction,
        "remove_corona_ions": remove_corona_ions,
        "remove_noisy_electrometers": remove_noisy_electrometers,
        "allow_reprocess": allow_reprocessing,
        "apply_cleaning": apply_cleaning
    }

    # Save the config file
    with open(config_file,"w") as cf:
        yaml.dump(config_info,cf)


# Inlet losses
################################################################################
def tubeloss(dpp,pflow,plength,temp,press):
    """ Laminar diffusion losses in circular straight conduit """
    DPP,TEMP = np.meshgrid(dpp,temp)
    DPP,PRESS = np.meshgrid(dpp,press)
    DPP,PFLOW = np.meshgrid(dpp,pflow)
    rmuu = np.pi*af.particle_diffusivity(DPP,TEMP,PRESS)*plength/PFLOW
    pene = np.nan*np.ones(rmuu.shape)
    cond1=rmuu<0.02
    cond2=rmuu>=0.02
    pene[cond1] = 1. - 2.56*rmuu[cond1]**(2./3.) + 1.2*rmuu[cond1]+0.177*rmuu[cond1]**(4./3.)
    pene[cond2] = 0.819*np.exp(-3.657*rmuu[cond2]) + 0.097*np.exp(-22.3*rmuu[cond2]) + 0.032*np.exp(-57.0*rmuu[cond2])
    return pene

# Read raw data file into a dataframe
################################################################################
def read_file(fn):
    """
    Read NAIS raw data file into a pandas DataFrame

    Parameters
    ----------

    fn : str
        Raw data filename

    Returns
    -------

    pandas.DataFrame
        Contents of the file in dataframe

    """

    with open(fn,'r') as f:

        header_found = False
        data_matrix=[]
        lines = f.read().splitlines()
        
        for line in lines:

             # Skip empty and comments
             if (len(line)==0):
                 continue

             if (line[0]=='#'):
                 continue

             # Test if it is a header
             elif (header_found==False):
                 if "opmode" in line:
                     delimiter = re.search('(.)opmode',line).group(1)
                     header = line.split(delimiter)
                     number_of_columns = len(header)
                     header_found = True
                     continue
                 else:
                     continue 
             else:
                 data_line = line.split(delimiter)
                 
                 if ((len(data_line)==number_of_columns) & ("opmode" not in data_line)):
                     data_matrix.append(data_line)
                 continue

    if len(data_matrix)==0:
        return None

    else:
        # Convert anything that can be converted to float and the rest is coerced to NaNs
        df = pd.DataFrame(columns = header, data = data_matrix)
        df.iloc[:,3:] = df.iloc[:,3:].apply(pd.to_numeric, errors='coerce').astype(float)

        # Establish begin_time (first column) as index
        df = df.set_index(df.columns[0])
        df.index = pd.to_datetime(df.index)
        
        return df

# Average data into the standard size bins
################################################################################
def average_mob(y,h):
    data = pd.DataFrame([])
    
    for i in range(0,len(mob_ion_geomeans)):
        if i==0:
            y_block = y.iloc[:,h>mob_ion_geomeans[i]]
        else:
            y_block = y.iloc[:,((h>mob_ion_geomeans[i]) & (h<=mob_ion_geomeans[i-1]))]

        data[i] = y_block.median(axis=1)

    y_block = y.iloc[:,h<=mob_ion_geomeans[i]]
    data[i+1] = y_block.mean(axis=1)

    return data

def average_dp(y,h):
    data = pd.DataFrame([])

    for i in range(0,len(dp_par_geomeans)):
        if i==0:
            y_block = y.iloc[:,h<dp_par_geomeans[i]]
        else:
            y_block = y.iloc[:,((h<dp_par_geomeans[i]) & (h>=dp_par_geomeans[i-1]))]

        data[i] = y_block.median(axis=1)

    y_block = y.iloc[:,h>=dp_par_geomeans[i]]
    data[i+1] = y_block.mean(axis=1)

    return data

# Find diagnostic names
################################################################################
def find_diagnostic_names(diag_params):

    sampleflow_name=None
    temperature_name=None
    pressure_name=None

    for temp_name in possible_temperature_names:
         if temp_name in diag_params:
             temperature_name = temp_name
             break

    for pres_name in possible_pressure_names:
        if pres_name in diag_params:
            pressure_name = pres_name
            break

    # try single flow sensor
    for flow_name in possible_sampleflow_names1:
        if flow_name in diag_params:
            sampleflow_name = flow_name
            break

    if sampleflow_name is None:
        # try two flow sensors
        sf_name = []
        for flow_name in possible_sampleflow_names2:
            if flow_name in diag_params:
                sf_name.append(flow_name)
        if len(sf_name)==2:
            sampleflow_name=sf_name

    return temperature_name, pressure_name, sampleflow_name

# Process the data (convert to dndlogdp & corrections)
################################################################################
def process_data(df,mode):

    if (df is None):
        return None, None

    elif not df.index.to_series().is_monotonic_increasing:
        return None, None

    else:
        df_columns = df.columns
        df_inverter_reso = int((len(df_columns)-2)/4)

        neg_df = df.iloc[:,2:2+df_inverter_reso]
        pos_df = df.iloc[:,2+2*df_inverter_reso:2+3*df_inverter_reso]

        if mode=="ions":
            mob_ion_inv = np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0])
                                    for y in df_columns[2:2+df_inverter_reso]])

            neg_df = average_mob(neg_df,mob_ion_inv)
            pos_df = average_mob(pos_df,mob_ion_inv)

            # Convert to number size distributions
            neg_df = neg_df * dlogmob_ion / dlogdp_ion
            pos_df = pos_df * dlogmob_ion / dlogdp_ion

        if mode=="particles":
            dp_par_inv = 2.0*np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0])
                                       for y in df_columns[2:2+df_inverter_reso]])
        
            neg_df = average_dp(neg_df,dp_par_inv)
            pos_df = average_dp(pos_df,dp_par_inv)

        # Construct the headers
        if mode=="ions":
            df_header = dp_ion*1e-9
        if mode=="particles":
            df_header = dp_par*1e-9

        negdf = pd.DataFrame(columns=df_header, index=df.index, data=neg_df.values)
        posdf = pd.DataFrame(columns=df_header, index=df.index, data=pos_df.values)

        negdf.index.name = "Time"
        posdf.index.name= "Time"

        if negdf.isna().all().all():
            negdf = None
        if posdf.isna().all().all():
            posdf = None

        return negdf, posdf

def correct_data(
    df,
    rec,
    mode,
    do_sealevel_corr,
    pipe_length):

    if ((rec is None) or (df is None)):
        return None

    else:        
        # Extract the records that match the mode
        if mode=="ions":
            df_rec = rec[rec.opmode=='ions']
        if mode=="particles":
            df_rec = rec[rec.opmode=='particles']

        if not df_rec.index.to_series().is_monotonic_increasing:
            return None
        
        df_rec = df_rec.reindex(df.index,method="nearest")

        # Check that the relevant diagnostic data is found
        t_name,p_name,sf_name = find_diagnostic_names(list(df_rec))
        if ((t_name is not None) & 
            (p_name is not None) &
            (sf_name is not None)):
            pass
        else:
            return None
    
        # Temperature
        t_df = 273.15 + pd.DataFrame(df_rec[t_name].astype(float))

        # Pressure
        p_df = 100.0 * pd.DataFrame(df_rec[p_name].astype(float))
    
        # Sampleflow
        if len(sf_name)==2:
            flow_df = pd.DataFrame(df_rec[sf_name].sum(axis=1,min_count=2).astype(float))
        else:
            flow_df = pd.DataFrame(df_rec[sf_name].astype(float))
    
        # Test if the sampleflow is in cm3/s (old models) or 
        # l/min and possibly convert to l/min
        if (np.nanmedian(flow_df)>300):
            flow_df = (flow_df/1000.0) * 60.0
        else:
            pass
    
        # If all parameters are NaN
        # e.g. sensor is broken
        if (flow_df.isna().all().all() |
            p_df.isna().all().all() |
            t_df.isna().all().all()):
            return None
    
        # Sanity check the values
        t_df = t_df.where(((t_df>=223.)|(t_df<=353.)),np.nan)
        p_df = p_df.where(((p_df>=37000.)|(p_df<=121000.)),np.nan)
        flow_df = flow_df.where(((flow_df>=48.)|(flow_df<=60.)),np.nan)
    
        # Correct the number concentrations to standard conditions (optional)
        if (do_sealevel_corr):
            stp_corr_df = (pres_ref*t_df.values)/(temp_ref*p_df.values)
            df = stp_corr_df * df
       
        # Diffusion loss correction
        if mode=="ions":
            throughput = tubeloss(dp_ion*1e-9,flow_df.values*1.667e-5,pipe_length,t_df.values,p_df.values)
        if mode=="particles":
            throughput = tubeloss(dp_par*1e-9,flow_df.values*1.667e-5,pipe_length,t_df.values,p_df.values)
        
        df = df / throughput
    
        # Robert Wagner's calibration (only ions)
        if mode=="ions":
            roberts_corr = 0.713*dp_ion**0.120
            df = df / roberts_corr
    
        return df

# Data clean-up
################################################################################
def data_cleaner(df):
    """ Returns a cleaned data array and portion of data removed """

    if df is None:
        return None

    # Rolling time window
    reso_in_seconds = (df.index[1]-df.index[0]).seconds
    small_window = int((10.*60.)/(reso_in_seconds))
    large_window = int((24.*60.*60.)/(reso_in_seconds))

    # Calculate standard deviation in 10 min segments
    df2=df.rolling(small_window, min_periods=int((small_window+1.)/2.), center=True).std()

    # In a bigger window (24 hours) calculate the 75th quantile of the standard deviations
    # (semi)continous noise causes higher values compared to normal and rare sudden changes in conc
    df3=df2.rolling(large_window, min_periods=int((large_window+1.)/2.), center=True).quantile(0.75)

    # a Good threshold for significant noise seems to be 5 times the median background
    threshold = 5*np.nanmedian(df3)

    df4 = df.where(df3 < threshold, np.nan)

    return df4

def corona_limit(df):
    """ Find corona ion upper limit using maximum concentration difference """

    # Only consider likely limit range
    lower = 1.5e-9
    upper = 5.0e-9
    c = (lower <= df.columns.values) & (upper >= df.columns.values)
    df2 = df.loc[:, c]

    # Find maximum difference between size bin medians
    return df2.columns.values[df2.median().diff().abs().argmax()]

def corona_ion_cleaner(df):
    """ Return a cleaned data array and ratio of data removed """

    if df is None:
        return None

    corona_lim = corona_limit(df)
    df2 = df.copy()

    # Set values below corona ion limit to NaNs
    df2.iloc[:,df2.columns.values<=corona_lim]=np.nan

    return df2

def check_config(f):
    """
    Check that config file is ok

    Parameters
    ----------

    f : `str`
        full path to the configuration file

    Returns
    -------

    boolean
        `True` if file is OK
        `False` if file is not OK

    """

    if not os.path.isfile(f):
        print("Config not found")
        return False

    print(f)

    with open(f,'r') as stream:
        try:
            config = yaml.safe_load(stream)
            load_path = config['data_folder']
            save_path = config['processed_folder']
            start_date = config['start_date']
            database = config['database_file']
            location = config['location']
            end_date = config['end_date']
            allow_reprocess = config["allow_reprocess"]
            pipelength = config['inlet_length']
            sealevel_correction = config['sealevel_correction']
            apply_corrections = config['apply_corrections']
            apply_cleaning=config["apply_cleaning"]
            remove_noisy_electrometers = config["remove_noisy_electrometers"]
            remove_corona_ions = config["remove_corona_ions"]
        except:
            print("Config badly formatted")
            return False

    try:
      db = TinyDB(database)
      check = Query()
    except:
        print("Could not init DB")
        return False

    if start_date=='':
        pass
    elif isinstance(start_date,date):
        pass
    else:
        print("Bad start_date")
        return False

    if end_date=='':
        pass
    elif isinstance(end_date,date):
        pass
    else:
        print("Bad end_date")
        return False

    if os.path.exists(save_path)==False:
        print("Bad save path")
        return False

    for x in load_path:
        if os.path.exists(x)==False:
            print("Bad load path")
            return False

    if (allow_reprocess==True) or (allow_reprocess==False):
        pass
    else:
        print("Bad allow_reprocess")
        return False

    if (remove_corona_ions==True) or (remove_corona_ions==False):
        pass
    else:
        print("Bad remove_corona_ions")
        return False

    if (remove_noisy_electrometers==True) or (remove_noisy_electrometers==False):
        pass
    else:
        print("Bad remove_noisy_electrometers")
        return False

    if (sealevel_correction==True) or (sealevel_correction==False):
        pass
    else:
        print("Bad sealevel_correction")
        return False

    if (apply_cleaning==True) or (apply_cleaning==False):
        pass
    else:
        print("Bad apply_cleaning")
        return False

    if (apply_corrections==True) or (apply_corrections==False):
        pass
    else:
        print("Bad apply_corrections")
        return False

    try:
        float(pipelength)
    except:
        print("Bad inlet_length")
        return False

    return True

def nais_processor(config_file):
    """ Function that is called to processes data from the NAIS

    Parameters
    ----------

    config_file : str
        full path to the configuration file

    """

    ################# READING CONFIG
    if not check_config(config_file):
        return 

    # Today
    today_dt = date.today()

    with open(config_file,'r') as stream:
        config = yaml.safe_load(stream)
        load_path = config['data_folder']
        save_path = config['processed_folder']
        start_date = config['start_date']
        database = config['database_file']
        location = config['location']
        end_date = config['end_date']
        allow_reprocess = config["allow_reprocess"]
        pipelength = config['inlet_length']
        sealevel_correction = config['sealevel_correction']
        apply_corrections = config['apply_corrections']
        apply_cleaning=config["apply_cleaning"]
        remove_noisy_electrometers = config["remove_noisy_electrometers"]
        remove_corona_ions = config["remove_corona_ions"]


    ##################### UPDATING DATABASE
    if start_date=='':
        start_date = date(2000,1,1)
    if end_date=='':
        end_date = today_dt

    db = TinyDB(database)
    check = Query()

    start_dt=pd.to_datetime(start_date)
    end_dt=pd.to_datetime(end_date)

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    # list existing dates
    list_of_existing_dates = [x["timestamp"] for x in db.search(check.diagnostics.exists())]

    if len(list_of_existing_dates)==0:
        print("building database...")
        list_of_datetimes = pd.date_range(start=start_date_str, end=end_date_str)
    else:
        last_existing_date = sorted(list_of_existing_dates)[-1]
        list_of_datetimes = pd.date_range(start=last_existing_date, end=end_date_str)
    
    # Add unprocessed datafiles to the database
    for x in list_of_datetimes:
        if (x.strftime("%Y%m%d") in list_of_existing_dates):
            continue
        else:
            files_found=False
            for z in load_path:
                for y in filename_formats:

                    ion_fn = os.path.join(z,x.strftime(y[0]))
                    particle_fn = os.path.join(z,x.strftime(y[1]))
                    diagnostic_fn = os.path.join(z,x.strftime(y[2]))

                    if ( (os.path.exists(ion_fn) | # ions
                         os.path.exists(particle_fn)) & # particles
                         os.path.exists(diagnostic_fn) # diagnostics
                       ):

                        dtstr = x.strftime("%Y%m%d")

                        db.insert(
                            {"timestamp":dtstr,
                            "diagnostics":diagnostic_fn}
                            )

                        if os.path.exists(ion_fn):
                            db.update(
                                {"ions":ion_fn},
                                check.timestamp==dtstr)

                        if os.path.exists(particle_fn):
                            db.update(
                                {"particles":particle_fn},
                                check.timestamp==dtstr)

                        files_found=True
                        break

                if files_found:
                    break

    # From the database find the last day with processed data
    processed_days = db.search(
        check.processed_neg_ion_file.exists() |
        check.processed_pos_ion_file.exists() |
        check.processed_neg_particle_file.exists() |
        check.processed_pos_particle_file.exists())

    if len(processed_days)!=0:
        last_day=np.max([datetime.strptime(x["timestamp"],"%Y%m%d") for x in processed_days]).strftime("%Y%m%d")
    else:
        last_day=None

    ############## PROCESSING
    # reprocess data in db
    if allow_reprocess:
        iterator1 = iter(db.search(
         (check.diagnostics.exists() &
          (check.ions.exists() |
          check.particles.exists()) &
          (check.timestamp>=start_date_str) &
          (check.timestamp<=end_date_str))))
    else:
        iterator1 = iter(db.search(
            ((check.timestamp==last_day) &
             (check.timestamp>=start_date_str) &
             (check.timestamp<=end_date_str)) |
            (check.diagnostics.exists() &
             (check.ions.exists() |
             check.particles.exists()) &
             ~check.processed_neg_ion_file.exists() &
             ~check.processed_pos_ion_file.exists() &
             ~check.processed_neg_particle_file.exists() &
             ~check.processed_pos_particle_file.exists() &
             (check.timestamp>=start_date_str) &
             (check.timestamp<=end_date_str))))

    for x in iterator1:

        print("processing %s (%s)" % (x["timestamp"],location))

        ions_exist=bool(db.search(
            check.ions.exists() &
            (check.timestamp==x["timestamp"])))

        particles_exist=bool(db.search(
            check.particles.exists() &
            (check.timestamp==x["timestamp"])))

        records = read_file(x["diagnostics"])

        # ions
        if ions_exist:

            ions = read_file(x["ions"])

            negion_datamatrix,posion_datamatrix = process_data(ions,"ions")

            if apply_corrections:
                negion_datamatrix = correct_data(
                       negion_datamatrix,
                       records,
                       "ions",
                       sealevel_correction,
                       pipelength)
    
                posion_datamatrix = correct_data(
                       posion_datamatrix,
                       records,
                       "ions",
                       sealevel_correction,
                       pipelength)

            if apply_cleaning:
                if remove_noisy_electrometers:
                    negion_datamatrix = data_cleaner(negion_datamatrix)
                    posion_datamatrix = data_cleaner(posion_datamatrix)

            if (negion_datamatrix is not None):
                my_save_path_neg=os.path.join(save_path,"NAISn"+x["timestamp"]+"nds.sum")
                negion_datamatrix.to_csv(my_save_path_neg)
                db.update({"processed_neg_ion_file": my_save_path_neg},
                    check.timestamp==x["timestamp"])

            if (posion_datamatrix is not None):
                my_save_path_pos=os.path.join(save_path,"NAISp"+x["timestamp"]+"nds.sum")
                posion_datamatrix.to_csv(my_save_path_pos)
                db.update({"processed_pos_ion_file": my_save_path_pos},
                    check.timestamp==x["timestamp"])

        # particles
        if particles_exist:

            # Process particles
            particles = read_file(x["particles"])

            negpar_datamatrix,pospar_datamatrix = process_data(particles,"particles")

            if apply_corrections:
                negpar_datamatrix = correct_data(
                       negpar_datamatrix,
                       records,
                       "particles",
                       sealevel_correction,
                       pipelength)
    
                pospar_datamatrix = correct_data(
                       pospar_datamatrix,
                       records,
                       "particles",
                       sealevel_correction,
                       pipelength)

            if apply_cleaning:
                if remove_corona_ions:
                    negpar_datamatrix = corona_ion_cleaner(negpar_datamatrix)
                    pospar_datamatrix = corona_ion_cleaner(pospar_datamatrix)
    
                if remove_noisy_electrometers:
                    negpar_datamatrix = data_cleaner(negpar_datamatrix)
                    negpar_datamatrix = data_cleaner(pospar_datamatrix)
 
            if (pospar_datamatrix is not None):
                my_save_path_neg=os.path.join(save_path,"NAISn"+x["timestamp"]+"np.sum")
                negpar_datamatrix.to_csv(my_save_path_neg)
                db.update({"processed_neg_particle_file": my_save_path_neg},
                    check.timestamp==x["timestamp"])

            if (negpar_datamatrix is not None):
                my_save_path_pos=os.path.join(save_path,"NAISp"+x["timestamp"]+"np.sum")
                pospar_datamatrix.to_csv(my_save_path_pos)
                db.update({"processed_pos_particle_file": my_save_path_pos},
                    check.timestamp==x["timestamp"])

    print("Done!")


def combine_databases(database_list, combined_database):
    """Combine JSON databases

    If the measurement setup changes one may have to use multiple configuration files
    which results in multiple databases. With this function you can combine the databases
    into a single database after processing.

    Parameters
    ----------

    database_list : `str`
        List of full paths to databases that should be combined

        First database should have the earliest data, second database
        the second earliest and so on

    combined_database : `str`
        full path to combined database
    
    """

    DB = {}
    i = 0

    for database in database_list:

        fid=open(database)

        database_json=json.load(fid)

        for key in database_json["_default"]:
            DB[i] = database_json["_default"][key]
            i=i+1

    with open(combined_database, "w") as f:
        json.dump({"_default":DB},f)

def combine_spectra(
    database_file,
    begin_time,
    end_time,
    spectrum_type="negion",
    reso=60):
    """
    Combine (cleaned) processed particle or ion data from some time range

    Parameters
    ----------

    database_file : `str`
        full path to database_file

    begin_time : `str`
        time zone aware iso formatted time string

        For example `"2013-01-02 15:00:00+02:00"`

    end_time : `str`
        time zone aware iso formatted time string

        For example `"2013-01-03 17:00:00+02:00"`

    spectrum_type : `str`
        negative ions `negion` (default)

        positive ions `posion`

        negative particles `negpar`

        positive particles `pospar`

    reso : `int`
        desired resolution given in minutes

    Returns
    -------

    pandas.DataFrame
        Combined aerosol number size distribution in the given 
        time interval

    """

    db = TinyDB(database_file)
    check = Query()

    begin_dt=pd.to_datetime(begin_time)
    end_dt=pd.to_datetime(end_time)

    begin_date=begin_dt.strftime("%Y%m%d")
    end_date=end_dt.strftime("%Y%m%d")

    if spectrum_type=="negpar":
        iterator = iter(db.search(
            (check.processed_neg_particle_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date)))
        db_entry = "processed_neg_particle_file"
    elif spectrum_type=="pospar":
        iterator = iter(db.search(
            (check.processed_pos_particle_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date)))
        db_entry = "processed_pos_particle_file"
    elif spectrum_type=="negion":
        iterator = iter(db.search(
            (check.processed_neg_ion_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date)))
        db_entry = "processed_neg_ion_file"
    elif spectrum_type=="posion":
        iterator = iter(db.search(
            (check.processed_pos_ion_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date)))
        db_entry = "processed_pos_ion_file"
    else:
        print("ERROR: %s is not valid 'spectrum_type'" % spectrum_type)
        return

    filenames = [x[db_entry] for x in iterator]

    df = af.stack_data(filenames,begin_time,end_time,reso)

    return df
