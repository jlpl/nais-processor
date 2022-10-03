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
import time
import json

__pdoc__ = {
    'visc': False,
    'rlambda': False,
    'cunn': False,
    'diffuus': False,
    'tubeloss': False,
    'datetime2datenum': False,
    'datenum2datetime': False,
    'log_ticks': False,
    'average_mob': False,
    'average_dp': False,
    'find_diagnostic_names': False,
    'process_data': False,
    'average_dp': False,
    'average_dp': False,
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
['%Y-%m-%d.ions.nds','%Y-%m-%d.particles.nds','%Y-%m-%d.log'],
['%Y%m%d-block-ions.spectra','%Y%m%d-block-particles.spectra','%Y%m%d-block.records'],
['%Y%m%d-block-ions.spectra','%Y%m%d-block-particles.spectra','%Y%m%d-block.diagnostics']]

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
'temp']

possible_pressure_names = [
'baro.mean',
'baro']

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
    print("Enter name of configuration file.")
    print("E.g. ./configs/campaign.yml")
    while True:
        config_file = input("> ")
        if len(config_file)>0:
            if not config_file.lower().endswith(".yml"):
                config_file = config_file + ".yml"
            break
        else:
            continue

    print()
    print("Give path(s) to raw data. If multiple paths give them as comma separated list.")
    print("E.g. /data/nais/2021,/data/nais/2022")
    while True:
        user_input = input("> ")
        if len(user_input)>0:
            load_path=user_input.split(",")
            break
        else:
            continue

    print()
    print("Path to where processed data is saved.")
    print("E.g. ./data/campaign/processed")
    while True:
        save_path = input("> ")
        if len(save_path)>0:
            break
        else:
            continue

    print()
    print("Path to where figures are saved. Leave empty if no figures.")
    print("E.g. ./data/campaign/figures")
    fig_path = input("> ")

    print()
    print("Start of measurement (YYYY-MM-DD)")
    while True:        
        start_date = input("> ")
        try:
            start_dt = pd.to_datetime(start_date)
            break
        except:
            continue

    print()
    print("End of measurement (YYYY-MM-DD)")
    print("If empty processor assumes current day, use for continuous measurements.")
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
    print("Enter name of database file")
    print("E.g. ./logs/campaign.json")
    while True:
        database_file = input("> ")
        if len(database_file)>0:
            if not database_file.lower().endswith(".json"):
                database_file = database_file + ".json"
            break
        else:
            continue

    print()
    print("Allow reprocessing (True/False)")
    print("Overwrites already processed datafiles in the database when running the processor again.")
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
    print("Measurement location")
    print("E.g. Helsinki, Kumpula")
    location = input("> ")

    print()
    print("Apply corrections to data? (True/False)")
    print("Requires a NAIS with temperature and pressure sensors.")
    while True:
        apply_corrections = input("> ")
        if ((apply_corrections=='True') or (apply_corrections=='False')):
            if apply_corrections=='True':
                apply_corrections=True
            else:
                apply_corrections=False
                sealevel_correction=False
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
        print("Correct concentrations to sealevel conditions? (True/False)")
        while True:
            sealevel_correction = input("> ")
            if ((sealevel_correction=='True') or (sealevel_correction=='False')):
                if sealevel_correction=='True':
                    sealevel_correction=True
                else:
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
        "figure_folder": fig_path,
        "start_date": start_date,
        "end_date": end_date,
        "database_file": database_file,
        "location": location,
        "inlet_length": inlet_length,
        "apply_corrections":apply_corrections,
        "sealevel_correction": sealevel_correction,
        "allow_reprocess": allow_reprocessing 
    }

    # Save the config file
    with open(config_file,"w") as cf:
        yaml.dump(config_info,cf)


# Inlet losses
################################################################################
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
    pene[cond1] = 1. - 2.56*rmuu[cond1]**(2./3.) + 1.2*rmuu[cond1]+0.177*rmuu[cond1]**(4./3.)
    pene[cond2] = 0.819*np.exp(-3.657*rmuu[cond2]) + 0.097*np.exp(-22.3*rmuu[cond2]) + 0.032*np.exp(-57.0*rmuu[cond2])
    return pene


# Time format conversion
################################################################################
def datetime2datenum(dt):
    mdn = dt + timedelta(days = 366)
    frac = (dt-datetime(dt.year,dt.month,dt.day,0,0,0,tzinfo=dt.tzinfo)).seconds \
           / (24.0 * 60.0 * 60.0)
    return mdn.toordinal() + frac

def datenum2datetime(datenum):
    return datetime.fromordinal(int(datenum)) + timedelta(days=datenum%1) - timedelta(days = 366)

# Plotting
################################################################################
def log_ticks():
    x=np.arange(1,10)
    y=np.arange(-10,-4).astype(float)
    tt=[]
    tts=[]
    for j in y:
        for i in x:
            tt.append(np.log10(np.round(i*10**j,int(np.abs(j)))))
            if i==1:
                tts.append("10$^{%d}$"%j)
            else:
                tts.append('')
    
    tt=np.array(tt)
    return tt,tts
 
def plot_sumfile(handle,v,clims=(10,100000),hour_step=2,date_formatter="%Y-%m-%d %H:%M"):
    """ 
    Plot sum-formatted aerosol number size distribution 

    Parameters
    ----------

    handle : matplotlib.Axes.axes
        Axes object on which to plot

    v : 2-d np.array
        Aerosol size distribution matrix

    clims : list/tuple, size 2
        Data limits for the color

    hour_step : int
        Resolution on the time axis given in hours

    date_formatter : str
        Date format string

    """

    time = v[1:,0]
    dp = v[0,2:]
    data = v[1:,2:]
 
    tt,tts = log_ticks()

    norm = colors.LogNorm(vmin=clims[0],vmax=clims[1])
    color_ticks = LogLocator(subs=range(10))
        
    handle.set_yticks(tt)
    handle.set_yticklabels(tts)
    
    t1=dts.date2num(datenum2datetime(time.min()))
    t2=dts.date2num(datenum2datetime(time.max()))
    
    img = handle.imshow(
        np.flipud(data.T),
        origin="upper",
        aspect="auto",
        interpolation="hanning",
        cmap="turbo",
        norm=norm,
        extent=(t1,t2,np.log10(dp.min()),np.log10(dp.max()))
    )
    
    handle.xaxis.set_major_locator(dts.HourLocator(interval=hour_step))
    handle.xaxis.set_major_formatter(dts.DateFormatter(date_formatter))
    plt.setp(handle.get_xticklabels(), rotation=80)
    
    box = handle.get_position()
    c_handle = plt.axes([box.x0*1.025 + box.width * 1.025, box.y0, 0.01, box.height])
    cbar = plt.colorbar(img,cax=c_handle,ticks=color_ticks)

    handle.set_ylabel('Dp, [m]')
    handle.set_xlabel('UTC'+'%+d'%v[0,0]+', [h]')
    cbar.set_label('dN/dlogDp, [cm-3]')

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
        
    pandas DataFrame
             
    """
    
    with open(fn) as f:

        header_found = False
        data_matrix=[]
        line = f.readline()
        while line:

             line = line.rstrip("\n")

             # Skip completely empty line
             if len(line)==0:
                 line = f.readline()
                 continue

             # Skip commented line
             if line[0]=='#':
                 line = f.readline()
                 continue

             # Extract header line
             if header_found==False:
                 delimiter = re.search('(.)opmode',line).group(1)
                 header = line.split(delimiter)
                 number_of_columns = len(header)
                 header_found = True
                 line = f.readline()
                 continue

             data_line = line.split(delimiter)

             # Add data line to data matrix only if correct number of columns
             # and not a header line if header is inserted mid file 
             # when the NAIS is restarted
             if ((len(data_line)==number_of_columns) & ("opmode" not in data_line)):
                 data_matrix.append(data_line)

             line = f.readline()
                 
    # Construct data frame
    # In the data convert anything that can be converted to float and the rest is coerced to NaNs
    df = pd.DataFrame(columns = header, data = data_matrix)
    df.iloc[:,3:] = df.iloc[:,3:].apply(pd.to_numeric, errors='coerce').astype(float)

    # Convert infinities to nan also in case any were present for some reason
    df.iloc[:,3:] = df.iloc[:,3:].replace([np.inf,-np.inf],np.nan)
              
    return df

# Average data into the standard size bins
################################################################################
def average_mob(y,h):
    data = np.nan*np.ones((y.shape[0],len(mob_ion)))

    for i in range(0,len(mob_ion_geomeans)):
        if i==0:
            y_block = y[:,h>mob_ion_geomeans[i]]
        else:
            y_block = y[:,((h>mob_ion_geomeans[i]) & (h<=mob_ion_geomeans[i-1]))]

        data[:,i] = np.nanmean(y_block,axis=1)
        
    y_block = y[:,h<=mob_ion_geomeans[i]]
    data[:,i+1] = np.nanmean(y_block,axis=1)

    return data

def average_dp(y,h):
    data = np.nan*np.ones((y.shape[0],len(dp_par)))

    for i in range(0,len(dp_par_geomeans)):
        if i==0:
            y_block = y[:,h<dp_par_geomeans[i]]
        else:
            y_block = y[:,((h<dp_par_geomeans[i]) & (h>=dp_par_geomeans[i-1]))]

        data[:,i] = np.nanmean(y_block,axis=1)
 
    y_block = y[:,h>=dp_par_geomeans[i]]
    data[:,i+1] = np.nanmean(y_block,axis=1)

    return data

def find_diagnostic_names(diag_params):
    
    sampleflow_name=None
    sampleflow_names=None
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

    sf_name = []
    for flow_name in possible_sampleflow_names:
        if flow_name in diag_params:
            sf_name.append(flow_name)

    if len(sf_name)==2:
        sampleflow_names = sf_name
    if len(sf_name)==1:
        sampleflow_name = sf_name

    return temperature_name, pressure_name, sampleflow_names, sampleflow_name

# Process the data (convert to dndlogdp & corrections)
################################################################################
def process_data(
    df,
    rec,
    mode,
    apply_corr,
    sealevel_corr,
    pipel):

    try:
        df = df.set_index(df.columns[0])
        df.index = [parse(y) for y in df.index]
        df_columns = df.columns
        df_inverter_reso = int((len(df_columns)-2)/4)
    
        # get the number densities and interpolate out the nans
        neg_df = df.iloc[:,2:2+df_inverter_reso+1].astype(float)
        neg_df = neg_df.interpolate().values
        pos_df = df.iloc[:,2+2*df_inverter_reso:2+3*df_inverter_reso+1].astype(float)
        pos_df = pos_df.interpolate().values
    
        if mode=="ions":
            mob_ion_inv = np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0]) 
                                    for y in df_columns[2:2+df_inverter_reso+1]])
            neg_df = average_mob(neg_df,mob_ion_inv)
            pos_df = average_mob(pos_df,mob_ion_inv)
    
            # get the number size distributions
            neg_df = neg_df * dlogmob_ion / dlogdp_ion
            pos_df = pos_df * dlogmob_ion / dlogdp_ion
         
        if mode=="particles":
            dp_par_inv = 2.0*np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0]) 
                                       for y in df_columns[2:2+df_inverter_reso+1]])
            neg_df = average_dp(neg_df,dp_par_inv)
            pos_df = average_dp(pos_df,dp_par_inv)
    
        # If all data is NaNs then skip
        if (np.all(np.isnan(neg_df)) | 
            np.all(np.isnan(pos_df)) ):
            return None
    
        if apply_corr:
            rec = rec.set_index('opmode')
    
            # If relevant diagnostic data to make corrections does not exist: 
            # return nothing
            t_name,p_name,sf_names,sf_name = find_diagnostic_names(list(rec))
            if ((t_name is not None) & (p_name is not None) & 
                ((sf_names is not None) | (sf_name is not None))):
                pass
            else:
                return None
    
            # Then extract the records that match the polarity
            if mode=="ions":
                df_rec = rec.loc['ions'].set_index(rec.columns[0])
            if mode=="particles":
                df_rec = rec.loc['particles'].set_index(rec.columns[0])
    
            df_rec.index = [parse(y) for y in df_rec.index]            
    
            # Match records to data according to time
            df_rec = df_rec.reindex(index=df.index,method='nearest')
    
            # Temperature
            t_df = 273.15 + df_rec[t_name].astype(float)
            t_df = t_df.interpolate().values.flatten()
    
            # Pressure
            p_df = 100.0 * df_rec[p_name].astype(float)
            p_df = p_df.interpolate().values.flatten()
           
            # Sampleflow
            if sf_names is not None:
                flow_df = df_rec[sf_names].astype(float).sum(axis=1)         
                flow_df = flow_df.interpolate().values.flatten()
            if sf_name is not None:
                flow_df = df_rec[sf_name].astype(float)
                flow_df = flow_df.interpolate().values.flatten()            
    
            # If any of the relevant data are all NaNs return nothing
            if (np.all(np.isnan(flow_df)) |
                np.all(np.isnan(p_df)) |
                np.all(np.isnan(t_df))):
                return None
    
            # Test if the sampleflow is in cm3/s (old models) or l/min and possibly convert to l/min
            if flow_df[0]>100:
                flow_df = (flow_df/1000.0) * 60.0
            else:
                pass
        
            # Correct the number concentrations to standard conditions (optional)
            if (sealevel_corr):
                stp_corr_df = (pres_ref*t_df)/(temp_ref*p_df)
                neg_df = (stp_corr_df*neg_df.T).T
                pos_df = (stp_corr_df*pos_df.T).T
        
            # Diffusion loss correction
            if mode=="ions":
                throughput_df = tubeloss(dp_ion*1e-9,flow_df*1.667e-5,pipel,t_df,p_df)
            if mode=="particles":
                throughput_df = tubeloss(dp_par*1e-9,flow_df*1.667e-5,pipel,t_df,p_df)
            neg_df = neg_df / throughput_df
            pos_df = pos_df / throughput_df
        
            # Robert Wagner's calibration (only ions)
            if mode=="ions":
                roberts_corr = 0.713*dp_ion**0.120
                neg_df = neg_df / roberts_corr
                pos_df = pos_df / roberts_corr
            
        # CREATE FINAL DATA MATRICES
    
        # Integrate total number concentrations
        if mode=="ions":
            total_neg_df = np.nansum(neg_df*dlogdp_ion,axis=1)[np.newaxis].T
            total_pos_df = np.nansum(pos_df*dlogdp_ion,axis=1)[np.newaxis].T
        if mode=="particles":
            total_neg_df = np.nansum(neg_df*dlogdp_par,axis=1)[np.newaxis].T
            total_pos_df = np.nansum(pos_df*dlogdp_par,axis=1)[np.newaxis].T      
    
        # Get the utc offset
        if df.index[0].utcoffset()==None:
            utc_offset = 0
        else:
            utc_offset = df.index[0].utcoffset().total_seconds()/(60.0*60.0)
    
        time_df = np.array([datetime2datenum(x) for x in df.index])[np.newaxis].T
     
        # Construct the headers
        if mode=="ions":
            df_header = np.insert(dp_ion*1e-9,0,(utc_offset,0))[np.newaxis]
        if mode=="particles": 
            df_header = np.insert(dp_par*1e-9,0,(utc_offset,0))[np.newaxis]
        
        # Construct the sum-files
        negdf = np.concatenate((df_header,np.concatenate((time_df,total_neg_df,neg_df),axis=1)))
        posdf = np.concatenate((df_header,np.concatenate((time_df,total_pos_df,pos_df),axis=1)))
    
        return [negdf,posdf]
    
    except:
        return None

def nais_processor(config_file):
    """ Function that is called to processes data from the NAIS

    Parameters
    ----------

    config_file : str
        full path to the configuration file    

    """

    # Find out today
    today_dt = datetime.today()
    today = today_dt.strftime('%Y%m%d')

    # Check that the config file exists
    if os.path.isfile(config_file)==False:
        print('"%s" does not exist' % config_file)
        return
    else:
        # Try to parse the config file
        with open(config_file,'r') as stream:
            try:
                config = yaml.safe_load(stream)

                load_path = config['data_folder']        
                save_path = config['processed_folder']
                start_date = config['start_date']
                database = config['database_file']
                location = config['location']
                end_date = config['end_date']
                ignore_db = config["allow_reprocess"]
                if len(end_date)==0:
                    end_date = today
                pipelength = config['inlet_length']
                sealevel_correction = config['sealevel_correction']
                apply_corrections = config['apply_corrections']
            except:
                print("Something went wrong with parsing %s",config_file)
                return

    # Then check if you can initialize the database
    try:
      db = TinyDB(database)
      check = Query()
    except:
        print("Could not initialize database")
        return
   
    # Test if the configuration information is valid
    try:
        float(pipelength)
    except:
        print('"%s" must be a number' % pipelength)
        return

    # Test if start and end dates are valid
    try:
       start_dt = pd.to_datetime(start_date)
       end_dt = pd.to_datetime(end_date)
    except:
       print('bad start_date or end_date')
       return

    # Check if given data folders exist
    for x in load_path:
        if os.path.exists(x):
            continue
        else:
            print("At least one data folder does not exist")
            return

    # Check for save path and create folders if they do not exist.
    if os.path.exists(save_path):
        pass
    else:
        print("Save path does not exist")
        return

    print("Configuration file: %s" % config_file)

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    # Convert load and save paths to absolute paths
    load_path = [os.path.abspath(x) + '/' for x in load_path]
    save_path = os.path.abspath(save_path) + '/'

    # Make a list of datetimes that the config file is interested in
    list_of_datetimes = pd.date_range(start=start_date_str, end=end_date_str)

    # list existing dates
    list_of_existing_dates = [x['timestamp'] for x in db.search(check.diagnostics.exists())]
    earliest_existing_date = np.min(pd.to_datetime(list_of_existing_dates))

    # Add unprocessed datafiles to the database
    for x in list_of_datetimes:
        if ((x.strftime('%Y%m%d') in list_of_existing_dates) | 
            (x < earliest_existing_date)):
            continue
        else:
            files_found=False
            for z in load_path:
                for y in filename_formats:

                    if ( (os.path.exists(z+x.strftime(y[0])) | # ions
                         os.path.exists(z+x.strftime(y[1]))) & # particles
                         os.path.exists(z+x.strftime(y[2])) # diagnostics
                       ):

                        db.insert(
                            {'timestamp':x.strftime('%Y%m%d'),
                            'diagnostics':z+x.strftime(y[2])}
                            )

                        if os.path.exists(z+x.strftime(y[0])):
                            db.update(
                                {'ions':z+x.strftime(y[0])},                               
                                check.timestamp==x.strftime('%Y%m%d'))

                        if os.path.exists(z+x.strftime(y[1])):
                            db.update(
                                {'particles':z+x.strftime(y[1])},                               
                                check.timestamp==x.strftime('%Y%m%d'))

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
        last_day=np.max([datetime.strptime(x['timestamp'],'%Y%m%d') for x in processed_days]).strftime('%Y%m%d')
    else:
        last_day=None
    
    # decide your data iterator

    # reprocess data in db
    if ignore_db:
        iterator =  iter(db.search( 
        ((check.timestamp==last_day) & 
         (check.timestamp>=start_date_str) &
         (check.timestamp<=end_date_str)) |
         (check.diagnostics.exists() &
          (check.ions.exists() |
          check.particles.exists()) &
          (check.timestamp>=start_date_str) &
          (check.timestamp<=end_date_str))))

    # do not reprocess data in db
    else:
        iterator =  iter(db.search( 
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
    
    for x in iterator: 
 
        print('processing %s' % x['timestamp'])

        # Read the diagnostics
        records = read_file(x['diagnostics'])

        ions_exist=db.search(
            check.ions.exists() & 
            (check.timestamp==x['timestamp']))
        particles_exist=db.search(
            check.particles.exists() & 
            (check.timestamp==x['timestamp']))


        # ions
        if ions_exist:

            ions = read_file(x['ions'])
            ion_datamatrices = process_data(
                 ions,
                 records,
                 "ions",
                 apply_corrections,
                 sealevel_correction,
                 pipelength)

            if ion_datamatrices is not None:
                
                # Save the sum matrices using the standard names
                np.savetxt(save_path+'NAISn'+x['timestamp']+'nds.sum',ion_datamatrices[0])
                np.savetxt(save_path+'NAISp'+x['timestamp']+'nds.sum',ion_datamatrices[1])
            
                # Update the database
                db.update(
                    {'processed_neg_ion_file': save_path+'NAISn'+x['timestamp']+'nds.sum',
                    'processed_pos_ion_file': save_path+'NAISp'+x['timestamp']+'nds.sum'},
                    check.timestamp==x['timestamp'])

        # particles
        if particles_exist:

            # Process particles
            particles = read_file(x['particles'])
            particle_datamatrices = process_data(
                particles,
                records,
                "particles",
                apply_corrections,
                sealevel_correction,
                pipelength)

            if particle_datamatrices is not None:

                # Save the sum matrices using the standard names
                np.savetxt(save_path+'NAISn'+x['timestamp']+'np.sum',particle_datamatrices[0])
                np.savetxt(save_path+'NAISp'+x['timestamp']+'np.sum',particle_datamatrices[1])
            
                # Update the database
                db.update(
                    {'processed_neg_particle_file': save_path+'NAISn'+x['timestamp']+'np.sum',
                    'processed_pos_particle_file': save_path+'NAISp'+x['timestamp']+'np.sum'},
                    check.timestamp==x['timestamp'])

    print("Done!")



def do_daily_figs(config_file):
    """ Make daily plots of NAIS data into the fig folder

    Parameters
    ----------

    config_file : str
        full path to the configuration file

    """

    # Find out today
    today_dt = datetime.today()
    today = today_dt.strftime('%Y%m%d')

    # Check that the config file exists
    if os.path.isfile(config_file) == False:
        print('"%s" does not exist' % config_file)
        return
    else:
        # Try to parse the config file
        with open(config_file,'r') as stream:
            try:
                config = yaml.safe_load(stream)
                save_path = config['processed_folder']
                database = config['database_file']
                location = config['location']
                fig_path = config['figure_folder']
                ignore_db = config['allow_reprocess']
                if len(fig_path)==0:
                    fig_path = None
            except:
                print("Something went wrong with parsing %s",config_file)
                return

    # Check if you can initialize the database
    try:
      db = TinyDB(database)
      check = Query()
    except:
        print("Could not initialize database")
        return

    # Check if processed data path exists
    if not os.path.exists(save_path):
        print("Path to processed data is invalid")
        return
        
    # Check the fig path exists
    if fig_path is not None:
        if not os.path.exists(fig_path):
            print('Figure path does not exist')
            return
    else:
        print("figure path not given")
        return 

    fig_path = os.path.abspath(fig_path) + '/'

    # Define some plotting styles
    plt.style.use('dark_background')

    fontsize = 14
    plt.rcParams.update({'font.size': fontsize,
                         'axes.titlesize': fontsize,
                         'axes.labelsize': fontsize,
                         'xtick.labelsize': fontsize,
                         'ytick.labelsize': fontsize,
                         'figure.titlesize': fontsize,
                         'legend.fontsize': fontsize})
 
    # From the database find the last day with processed data
    processed_days = db.search( 
        check.processed_neg_ion_file.exists() |
        check.processed_pos_ion_file.exists() |
        check.processed_neg_particle_file.exists() |
        check.processed_pos_particle_file.exists())

    if len(processed_days)!=0:
        last_day=np.max([datetime.strptime(x['timestamp'],'%Y%m%d') for x in processed_days]).strftime('%Y%m%d')
    else:
        last_day=None

    if ignore_db:
        iterator = iter(db.search((     
          (check.processed_neg_ion_file.exists() &
          check.processed_pos_ion_file.exists() &
          (check.timestamp==last_day)) |
          (check.processed_neg_particle_file.exists() &
          check.processed_pos_particle_file.exists() &
          (check.timestamp==last_day)) |
          (check.processed_neg_ion_file.exists() &
          check.processed_pos_ion_file.exists()) |
          (check.processed_neg_particle_file.exists() &
          check.processed_pos_particle_file.exists()))))

    else:
        iterator = iter(db.search((     
          (check.processed_neg_ion_file.exists() &
          check.processed_pos_ion_file.exists() &
          (check.timestamp==last_day)) |
          (check.processed_neg_particle_file.exists() &
          check.processed_pos_particle_file.exists() &
          (check.timestamp==last_day)) |
          (check.processed_neg_ion_file.exists() &
          check.processed_pos_ion_file.exists() &
          ~check.ion_figure.exists()) |
          (check.processed_neg_particle_file.exists() &
          check.processed_pos_particle_file.exists() &
          ~check.particle_figure.exists()))))
 

    for x in iterator: 

        print('plotting %s' % x['timestamp'])

        ions_exist=db.search(
            check.processed_neg_ion_file.exists() &
            check.processed_pos_ion_file.exists() &
            (check.timestamp==x['timestamp']))
        particles_exist=db.search(
            check.processed_neg_particle_file.exists() &
            check.processed_pos_particle_file.exists() &
            (check.timestamp==x['timestamp']))

        if ions_exist:
            negion=np.loadtxt(x["processed_neg_ion_file"])
            posion=np.loadtxt(x["processed_pos_ion_file"])
            fig,ax = plt.subplots(2,1,figsize=(8,7))
            ax = ax.flatten()
            plot_sumfile(ax[0],posion,clims=(10,10000),hour_step=1,date_formatter="%H:%M")
            plot_sumfile(ax[1],negion,clims=(10,10000),hour_step=1,date_formatter="%H:%M")
            ax[0].set_xticklabels([])
            ax[0].set_xlabel('')
            ax[0].set_title('Negative ions',loc="left")
            ax[1].set_title('Positive ions',loc="left")
            fig.suptitle(x['timestamp'] + ' ' + location)
            plt.savefig(fig_path+'NAIS_ions_'+ x['timestamp'] +'.png',dpi=100,bbox_inches='tight')
            db.update({'ion_figure': fig_path+'NAIS_ions_'+ x['timestamp'] +'.png'}, check.timestamp==x['timestamp'])
            plt.close()
 
        if particles_exist:
            negpar=np.loadtxt(x["processed_neg_particle_file"])
            pospar=np.loadtxt(x["processed_pos_particle_file"])
            fig,ax = plt.subplots(2,1,figsize=(8,7))
            ax = ax.flatten()
            plot_sumfile(ax[0],pospar,clims=(10,100000),hour_step=1,date_formatter="%H:%M")
            plot_sumfile(ax[1],negpar,clims=(10,100000),hour_step=1,date_formatter="%H:%M")
            ax[0].set_xticklabels([])
            ax[0].set_xlabel('')
            ax[0].set_title('Particles (positive polarity)',loc="left")
            ax[1].set_title('Particles (negative polarity)',loc="left")
            fig.suptitle(x['timestamp'] + ' ' + location)
            plt.savefig(fig_path+'NAIS_particles_'+ x['timestamp'] +'.png',dpi=100,bbox_inches='tight')
            db.update({'particle_figure': fig_path+'NAIS_particles_'+x['timestamp'] +'.png'}, check.timestamp==x['timestamp'])
            plt.close()

    print("Done!")

def combine_spectra(config_file,begin_time,end_time,spectrum_type="negion"):
    """
    Combine processed particle or ion data from some time range 
    
    Parameters
    ----------

    config_file : str
        full path to configuration file

    begin_time : str
        date string for begin time

    end_time : str
        date string for end time

    spectrum_type : str
        negative ions `negion` (default)

        positive ions `posion`

        negative particles `negpar`

        positive particles `pospar`


    Returns
    -------

    2-d array
        combined sum file between begin_time
        and end_time for the specified spectrum_type

    """

    if os.path.isfile(config_file) == False:
        print('"%s" does not exist' % config_file)
        return
    else:
        with open(config_file,'r') as stream:
            try:
                config = yaml.safe_load(stream)
                save_path = config['processed_folder']
                database = config['database_file']
            except:
                print("Something went wrong with parsing %s",config_file)
                return

    try:
      db = TinyDB(database)
      check = Query()
    except:
        print("Could not initialize database")
        return

    begin_dt=pd.to_datetime(begin_time)
    end_dt=pd.to_datetime(end_time)

    begin_dnum=datetime2datenum(begin_dt)
    end_dnum=datetime2datenum(end_dt)

    begin_date=begin_dt.strftime("%Y%m%d")
    end_date=end_dt.strftime("%Y%m%d")

    if spectrum_type=="negpar":
        iterator = iter(db.search(
            (check.processed_neg_particle_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date)))
    elif spectrum_type=="pospar":
        iterator = iter(db.search(
            (check.processed_pos_particle_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date)))
    elif spectrum_type=="negion":
        iterator = iter(db.search(
            (check.processed_neg_ion_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date)))
    elif spectrum_type=="posion":
        iterator = iter(db.search(
            (check.processed_pos_ion_file.exists()) &
            (check.timestamp>=begin_date) &
            (check.timestamp<=end_date))) 
    else:
        print("ERROR: %s is not valid 'spectrum_type'" % spectrum_type)
        return

    iter_num=1
    for x in iterator:
        if spectrum_type=="negpar":
            spectrum = np.loadtxt(x["processed_neg_particle_file"])
        if spectrum_type=="pospar":
            spectrum = np.loadtxt(x["processed_pos_particle_file"])
        if spectrum_type=="negion":
            spectrum = np.loadtxt(x["processed_neg_ion_file"])
        if spectrum_type=="posion":
            spectrum = np.loadtxt(x["processed_pos_ion_file"])

        data = spectrum[1:,:]

        if (iter_num==1):
            header = np.expand_dims(spectrum[0,:],axis=0)
            combined_spectrum = data
            iter_num = iter_num+1
        else:
            combined_spectrum = np.vstack((combined_spectrum,data))
            
    
    if iter_num==1:
        print("No data found")
        return
    else:
        findex = np.argwhere(
            (combined_spectrum[:,0]>=begin_dnum) & 
            (combined_spectrum[:,0]<=end_dnum)).flatten()

        return np.vstack((header,combined_spectrum[findex,:]))








