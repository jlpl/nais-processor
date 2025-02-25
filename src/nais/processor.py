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
import aerosol.functions as af
from scipy.interpolate import interp1d
import xarray as xr

# CONSTANTS
LEN_DP = 55
DP_STANDARD = np.array([
8.02879995e-10, 8.62828998e-10, 9.27254239e-10, 9.96153364e-10,
1.07017200e-09, 1.15119523e-09, 1.23835276e-09, 1.33004880e-09,
1.42853463e-09, 1.53530865e-09, 1.65006335e-09, 1.77383533e-09,
1.90689151e-09, 2.04892193e-09, 2.20153115e-09, 2.36722904e-09,
2.54539816e-09, 2.73599893e-09, 2.94087199e-09, 3.16011296e-09,
3.39569827e-09, 3.65377117e-09, 3.93145759e-09, 4.22385724e-09,
4.53800393e-09, 4.87892846e-09, 5.24546548e-09, 5.64126564e-09,
6.06693116e-09, 6.52192045e-09, 7.01103164e-09, 7.54292682e-09,
8.11517448e-09, 8.72847179e-09, 9.38811851e-09, 1.00955239e-08,
1.08562331e-08, 1.16916292e-08, 1.25913098e-08, 1.35414036e-08,
1.45631880e-08, 1.56757975e-08, 1.68734090e-08, 1.81717471e-08,
1.95699870e-08, 2.10714554e-08, 2.26881210e-08, 2.44557438e-08,
2.63610813e-08, 2.84161425e-08, 3.06314124e-08, 3.30248335e-08,
3.56052673e-08, 3.84642036e-08, 4.15526991e-08])

DP_STANDARD_NM = DP_STANDARD*1e9

MOB_STANDARD = np.array([
3.16000000e+00, 2.73662239e+00, 2.37000000e+00, 2.05390772e+00,
1.78000000e+00, 1.53862197e+00, 1.33000000e+00, 1.15324482e+00,
1.00000000e+00, 8.66015299e-01, 7.50000000e-01, 6.49221498e-01,
5.62000000e-01, 4.86987328e-01, 4.22000000e-01, 3.65167264e-01,
3.16000000e-01, 2.73658306e-01, 2.37000000e-01, 2.05387387e-01,
1.78000000e-01, 1.53859150e-01, 1.33000000e-01, 1.15321945e-01,
1.00000000e-01, 8.65992819e-02, 7.50000000e-02, 6.49201803e-02,
5.62000000e-02, 4.86970409e-02, 4.22000000e-02, 3.65152221e-02,
3.16000000e-02, 2.73645321e-02, 2.37000000e-02, 2.05376145e-02,
1.78000000e-02, 1.53848960e-02, 1.33000000e-02, 1.15313400e-02,
1.00000000e-02, 8.65916493e-03, 7.50000000e-03, 6.49134354e-03,
5.62000000e-03, 4.86911920e-03, 4.22000000e-03, 3.65099699e-03,
3.16000000e-03, 2.73599507e-03, 2.37000000e-03, 2.05336060e-03,
1.78000000e-03, 1.53812246e-03, 1.33000000e-03])

DLOGDP_STANDARD = np.array([
0.0312741 , 0.0312741 , 0.03120074, 0.03112738, 0.03141139,
0.0316954 , 0.0313593 , 0.0310232 , 0.03116406, 0.03130492,
0.0313588 , 0.03141268, 0.03130605, 0.03119943, 0.03135743,
0.03151544, 0.03143782, 0.0313602 , 0.03129334, 0.03122649,
0.03151937, 0.03181225, 0.03148394, 0.03115564, 0.03130761,
0.03145958, 0.03152605, 0.03159252, 0.03149947, 0.03140643,
0.0315822 , 0.03175797, 0.03169916, 0.03164035, 0.03159532,
0.0315503 , 0.03187309, 0.03219588, 0.03189432, 0.03159277,
0.03178298, 0.03197319, 0.03208352, 0.03219385, 0.03214892,
0.032104  , 0.03234317, 0.03258235, 0.03259212, 0.0326019 ,
0.03263776, 0.03267361, 0.03310805, 0.03354249, 0.03354249])

DLOGMOB_STANDARD = np.array([
0.0624749 , 0.06246942, 0.06231702, 0.06216423, 0.0627246 ,
0.06328425, 0.06260524, 0.0619259 , 0.06219808, 0.06246946,
0.06256644, 0.06266258, 0.06243771, 0.06221205, 0.06251308,
0.06281283, 0.06264173, 0.06246953, 0.06231759, 0.06216436,
0.06272553, 0.06328441, 0.06260579, 0.06192607, 0.06219918,
0.06246967, 0.06256765, 0.06266282, 0.06243888, 0.06221233,
0.06251488, 0.06281316, 0.06264342, 0.06246991, 0.06231955,
0.0621648 , 0.06272879, 0.06328495, 0.06260779, 0.06192666,
0.0622031 , 0.06247038, 0.06257201, 0.06266366, 0.06244327,
0.06221329, 0.06252155, 0.06281432, 0.06264986, 0.06247124,
0.06232717, 0.06216632, 0.06274126, 0.06328681, 0.06298631])

TOTAL_SAMPLEFLOW_NAMES = [
"sampleflow",
"Flowaer"]

POS_SAMPLEFLOW_NAMES = [
"pos_sampleflow.mean",
"pos_sampleflow"]

NEG_SAMPLEFLOW_NAMES = [
"neg_sampleflow.mean",
"neg_sampleflow"]

TEMPERATURE_NAMES = [
"temperature.mean",
"temperature",
"temp"]

PRESSURE_NAMES = [
"baro.mean",
"baro"]

RELHUM_NAMES = [
"relhum.mean",
"relhum"
]

DILUTION_FLOW_NAMES = [
"diluter_sample_flow_rate.mean",
"diluter_flow.mean"]

# Standard conditions
TEMP_REF = 293.15
PRES_REF = 101325.0

from nais import __version__
version=__version__

def make_config_template(file_name):
    """
    Make a configuration file template

    Parameters
    ----------

    file_name : str
        full path to configuration file. For example `/home/user/config.yml`

    """
    
    with open(file_name,"w") as f:
        f.write("measurement_location: # Name of the measurement site\n")
        f.write("id: # short id that identifies the measurement\n")
        f.write("description: # Additional description of the measurement\n")
        f.write("instrument_model: # e.g. NAIS-5-33\n")
        f.write("longitude: # decimal degrees west/east = -/+ (float) or null\n")
        f.write("latitude: # decimal degrees south/north = -/+ (float) or null\n")
        f.write("data_folder: # Full paths to raw data folders\n")
        f.write("- # Data folder 1\n")
        f.write("- # Data folder 2, and so on...\n")
        f.write("processed_folder: # Full path to folder where procesed data is saved\n")
        f.write("database_file: # Full path to database file (will be created on first run)\n")
        f.write("start_date: # Format: yyyy-mm-dd\n")
        f.write("end_date: # Format: yyyy-mm-dd or null for current day\n")
        f.write("inlet_length: # length of inlet in meters (float)\n")
        f.write("do_inlet_loss_correction: # true or false\n")
        f.write("convert_to_standard_conditions: # true or false\n")
        f.write("do_wagner_ion_mode_correction: # true or false\n")
        f.write("remove_corona_ions: # true or false\n")
        f.write("allow_reprocess: # true or false\n")
        f.write("redo_database: # true or false\n")
        f.write("fill_temperature: # null or temperature in K\n")
        f.write("fill_pressure: # null or pressure in Pa\n")
        f.write("fill_flowrate: # null or flow rate in lpm\n")
        f.write("dilution_on: # true or false (is the integrated dilution system used)\n")
        f.write('file_format: # 1s, 10s or block\n')
        f.write('resolution: # processed data time resolution (pandas time offset string), e.g. 5min')


def check_config_file(config_file):
    """ Check goodness of configuration file
    
    Parameters
    ----------
    
    config_file : str
        full path to configuration file

    """

    # Load configuration
    with open(config_file,'r') as stream:
        config = yaml.safe_load(stream)
        
        # Read in the configuration file
        load_path = config['data_folder']
        save_path = config['processed_folder']
        start_date = config['start_date']
        database = config['database_file']
        location = config['measurement_location']
        description = config['description']
        ide = config["id"]
        instrument_model = config['instrument_model']
        longitude = config["longitude"]
        latitude = config["latitude"]
        end_date = config['end_date']
        allow_reprocess = config["allow_reprocess"]
        redo_database = config["redo_database"]
        pipelength = config['inlet_length']
        do_inlet_loss_correction = config['do_inlet_loss_correction']
        convert_to_standard_conditions = config['convert_to_standard_conditions']
        do_wagner_ion_mode_correction = config["do_wagner_ion_mode_correction"]
        remove_charger_ions = config["remove_corona_ions"]
        file_format = config["file_format"]
        resolution = config["resolution"]
        fill_temperature = config["fill_temperature"]
        fill_pressure = config["fill_pressure"]
        fill_flowrate = config["fill_flowrate"]
        dilution_on = config["dilution_on"]
        
    # Check the config file
    assert isinstance(start_date,date)
    assert ((end_date is None) or isinstance(end_date,date))
    assert os.path.exists(save_path)
    assert all([os.path.exists(x) for x in load_path])
    assert isinstance(allow_reprocess,bool)
    assert isinstance(remove_charger_ions,bool)
    assert isinstance(convert_to_standard_conditions,bool)
    assert isinstance(do_wagner_ion_mode_correction,bool)
    assert isinstance(do_inlet_loss_correction,bool)
    assert isinstance(pipelength,float)
    assert (isinstance(longitude,float) or (longitude is None))
    assert (isinstance(latitude,float) or (latitude is None))
    assert (isinstance(fill_temperature,float) or (fill_temperature is None))
    assert (isinstance(fill_pressure,float) or (fill_pressure is None))
    assert (isinstance(fill_flowrate,float) or (fill_flowrate is None))
    assert isinstance(dilution_on,bool)
    assert isinstance(redo_database,bool)
    pd.tseries.frequencies.to_offset(resolution)


def read_raw(file_name,file_type,timestamp,resolution_str):
    with open(file_name,'r') as f:
        header_found = False
        data_matrix = []
        flag_explanations = []
        lines = f.read().splitlines()
        
        for line in lines:
            # Skip empty and comments
            if (len(line)==0):
                continue
            # Collect a list of flags and skip comments
            if line[:6]=="# flag":
                # The flags should have yaml format
                try:
                    diagnostic_comment_yaml = yaml.safe_load(line[7:].rstrip('\r\n'))
                    flag_name = list(diagnostic_comment_yaml.keys())[0]
                    flag_message = diagnostic_comment_yaml[flag_name]["message"]
                    flag_explanations.append([flag_name,flag_message])
                except:
                    continue
            elif (line[0]=='#'):
                continue
            else:
                pass
            # Test if it is a header
            if (header_found==False):
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

    # If no data was read we return nothing
    if len(data_matrix)==0:
        if file_type=="records":
            return None,None,None,None,None,None
        if file_type=="spectra":
            return None

    else:
        # Construct dataframes
        df = pd.DataFrame(columns = header, data = data_matrix, dtype = str)
        flag_explanations = pd.DataFrame(columns=["flag","message"], data=flag_explanations, dtype=str)
        
        # Remove duplicate flags (may happen due to restarts)        
        flag_explanations = flag_explanations[~flag_explanations["message"].duplicated()]
       
        data_tz = pd.to_datetime(df[df.columns[0]].loc[0]).tz

        # Transform time strings to utc aware datetime objects
        begin_time = pd.to_datetime(df[df.columns[0]].values, format="mixed", utc=True).tz_convert('UTC')
        end_time = pd.to_datetime(df[df.columns[1]].values, format="mixed", utc=True).tz_convert('UTC')

        center_time = begin_time + (end_time - begin_time)/2.
        df.index = center_time

        df = df.sort_index() # sort just in case
        df = df[~df.index.duplicated(keep='first')] # remove duplicates just in case
        
        # Still if something is wrong, just bail out 
        if not df.index.is_monotonic_increasing:
            if file_type=="records":
                return None,None,None,None,None,None
            if file_type=="spectra":
                return None
        
        # Define a standard time index for all data of the day
        if data_tz is None:
            standard_start_time = pd.to_datetime(timestamp).tz_localize('UTC')
        else:
            standard_start_time = pd.to_datetime(timestamp).tz_localize(data_tz).tz_convert("UTC")
        
        standard_end_time = standard_start_time + pd.Timedelta(days=1)
        standard_time = pd.date_range(start=standard_start_time, end=standard_end_time, freq=resolution_str, inclusive="left")
        
        if file_type=="records":
            # Extract records for ions, particles and offset
            ion_records_and_flags = df[df.opmode=='ions']
            particle_records_and_flags = df[df.opmode=='particles']
            offset_records_and_flags = df[df.opmode=='offset']
            
            if ion_records_and_flags.empty==False:
                ion_records = ion_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
                ion_flags = ion_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')
                
                ion_records = ion_records.resample(resolution_str).median()
                ion_flags = ion_flags.resample(resolution_str).apply(lambda col: ''.join(col))

                ion_records = ion_records.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))
                ion_flags = ion_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))

            else:
                ion_records = None
                ion_flags = None
                
            if particle_records_and_flags.empty == False:
                particle_records = particle_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
                particle_flags = particle_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')

                particle_records = particle_records.resample(resolution_str).median()
                particle_flags = particle_flags.resample(resolution_str).apply(lambda col: ''.join(col))

                particle_records = particle_records.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))
                particle_flags = particle_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))

            else:
                particle_records = None
                particle_flags = None
                
            if offset_records_and_flags.empty==False:
                offset_records = offset_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
                offset_flags = offset_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')

                offset_records = offset_records.resample(resolution_str).median()
                offset_flags = offset_flags.resample(resolution_str).apply(lambda col: ''.join(col))

                offset_records = offset_records.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))
                offset_flags = offset_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))

            else:
                offset_records = None
                offset_flags = None
                
            if flag_explanations.empty==True:
                flag_explanations = None

            return ion_records, particle_records, ion_flags, particle_flags, offset_flags, flag_explanations
        
        if file_type=="spectra":
            spectra = df.iloc[:,3:].apply(pd.to_numeric, errors='coerce').astype(float)

            spectra = spectra.resample(resolution_str).median()

            spectra = spectra.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))
             
            return spectra

def regrid_columns(data, old_columns, new_columns):
    data_out = interp1d(old_columns,data,bounds_error=False)(new_columns)
    return data_out

def raw2sum(spectra, mode):
    if (spectra is None):
        return None, None

    else:
        spectra_columns = spectra.columns
        spectra_inverter_reso = int(len(spectra_columns)/4)
        
        neg_arr = spectra.iloc[:,:spectra_inverter_reso].values
        pos_arr = spectra.iloc[:,2*spectra_inverter_reso:3*spectra_inverter_reso].values

        if mode == "ions":
            mob_ion_column = np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0])
                for y in spectra_columns[:spectra_inverter_reso]])

            # mobilities need to be flipped so they are monotonically increasing
            neg_arr = regrid_columns(np.fliplr(neg_arr),np.flip(mob_ion_column),np.flip(MOB_STANDARD))
            pos_arr = regrid_columns(np.fliplr(pos_arr),np.flip(mob_ion_column),np.flip(MOB_STANDARD))
            
            # Convert to number size distributions
            neg_arr = np.fliplr(neg_arr) * DLOGMOB_STANDARD / DLOGDP_STANDARD
            pos_arr = np.fliplr(pos_arr) * DLOGMOB_STANDARD / DLOGDP_STANDARD
            
        if mode == "particles":
            dp_par_column = 2.0*np.array([float(re.findall(r"[-+]?\d*\.\d+|\d+",y)[0])
                for y in spectra_columns[:spectra_inverter_reso]])*1e-9
            
            neg_arr = regrid_columns(neg_arr,dp_par_column,DP_STANDARD)
            pos_arr = regrid_columns(pos_arr,dp_par_column,DP_STANDARD)
        
        negdf = pd.DataFrame(columns=DP_STANDARD, index=spectra.index, data=neg_arr)
        posdf = pd.DataFrame(columns=DP_STANDARD, index=spectra.index, data=pos_arr)
        
        return negdf, posdf

def find_diagnostic_names(diag_params):
    sampleflow_name=None
    neg_sampleflow_name=None
    pos_sampleflow_name=None
    temperature_name=None
    pressure_name=None
    relhum_name=None
    dilution_flow_name=None

    for temp_name in TEMPERATURE_NAMES:
        if temp_name in diag_params:
            temperature_name = temp_name
            break

    for pres_name in PRESSURE_NAMES:
        if pres_name in diag_params:
            pressure_name = pres_name
            break

    for rh_name in RELHUM_NAMES:
        if rh_name in diag_params:
            relhum_name = rh_name
            break

    for flow_name in TOTAL_SAMPLEFLOW_NAMES:
        if flow_name in diag_params:
            sampleflow_name = flow_name
            break

    for pos_flow_name in POS_SAMPLEFLOW_NAMES:
        if pos_flow_name in diag_params:
            pos_sampleflow_name = pos_flow_name
            break

    for neg_flow_name in NEG_SAMPLEFLOW_NAMES:
        if neg_flow_name in diag_params:
            neg_sampleflow_name = neg_flow_name
            break

    for dilflow_name in DILUTION_FLOW_NAMES:
        if dilflow_name in diag_params:
            dilution_flow_name = dilflow_name 
            break

    return temperature_name, pressure_name, sampleflow_name, pos_sampleflow_name, neg_sampleflow_name, dilution_flow_name, relhum_name

def get_diagnostic_data(
    records,
    fill_pressure,
    fill_temperature,
    fill_flowrate):

    if records is None:
        return None, None, None, None, None, False, False, False
    else:
        (temperature_name,
            pressure_name,
            sampleflow_name,
            pos_sampleflow_name,
            neg_sampleflow_name,
            dilution_flow_name,
            relhum_name) = find_diagnostic_names(list(records))

        if temperature_name is not None:
            temperature = 273.15 + records[temperature_name].astype(float)
            temperature_filled = False
            # Values may be missing: e.g. sensor is broken
            if (temperature.isna().all() and (fill_temperature is not None)):
                temperature = pd.Series(index = records.index, dtype=float)
                temperature[:] = fill_temperature
                temperature_filled = True
        elif fill_temperature is not None:
            temperature = pd.Series(index = records.index, dtype=float)
            temperature[:] = fill_temperature
            temperature_filled = True
        else:
            temperature = None
            temperature_filled=False
    
        if pressure_name is not None:
            pressure = 100.0 * records[pressure_name].astype(float)
            pressure_filled = False
            if (pressure.isna().all() and (fill_pressure is not None)):
                pressure = pd.Series(index = pressure.index, dtype=float)
                pressure[:] = fill_pressure
                pressure_filled = True
        elif fill_pressure is not None:
            pressure = pd.Series(index = records.index, dtype=float)
            pressure[:] = fill_pressure
            pressure_filled=True
        else:
            pressure = None
            pressure_filled=False

        if relhum_name is not None:
            relhum = records[relhum_name].astype(float)
            if (relhum.isna().all()):
                relhum = None
        else:
            relhum = None

        if sampleflow_name is not None:
            sampleflow = records[sampleflow_name].astype(float)
            sampleflow_filled=False
            if (sampleflow.isna().all() and (fill_flowrate is not None)):
                sampleflow = pd.Series(index = records.index, dtype=float)
                sampleflow[:] = fill_flowrate
                sampleflow_filled=True
        elif ((neg_sampleflow_name is not None) and (pos_sampleflow_name is not None)):
            neg_sampleflow = records[neg_sampleflow_name].astype(float)
            pos_sampleflow = records[pos_sampleflow_name].astype(float)
            sampleflow = neg_sampleflow + pos_sampleflow
            sampleflow_filled=False
            if (sampleflow.isna().all() and (fill_flowrate is not None)):
                sampleflow = pd.Series(index = records.index, dtype=float)
                sampleflow[:] = fill_flowrate
                sampleflow_filled=True
        elif fill_flowrate is not None:
            sampleflow = pd.Series(index = records.index, dtype=float)
            sampleflow[:] = fill_flowrate
            sampleflow_filled=True
        else:
            sampleflow = None
            sampleflow_filled=False
    
        # Convert from cm3/s to lpm
        if (np.nanmedian(sampleflow)>300):
            sampleflow = (sampleflow/1000.0) * 60.0
        else:
            pass
    
        if dilution_flow_name is not None:
            dilution_flow = records[dilution_flow_name].astype(float)
            if dilution_flow.isna().all():
                dilution_flow = None
        else:
            dilution_flow = None
    
        # Sanity check the values
        if temperature is not None:
            temperature = temperature.where(((temperature>=223.)&(temperature<=353.)),np.nan)
        
        if pressure is not None:
            pressure = pressure.where(((pressure>=37000.)&(pressure<=121000.)),np.nan)
        
        if sampleflow is not None:
            sampleflow = sampleflow.where(((sampleflow>=48.)&(sampleflow<=65.)),np.nan)
        
        return temperature, pressure, sampleflow, dilution_flow, relhum, temperature_filled, pressure_filled, sampleflow_filled

def bring_to_sealevel(
    spectra,
    temperature,
    pressure):

    if ((spectra is None) or (temperature is None) or (pressure is None)):
        return None
    else:
        stp_corr_factor = (PRES_REF*temperature.values)/(TEMP_REF*pressure.values)
        spectra = stp_corr_factor.reshape(-1,1) * spectra
        
        return spectra

def correct_inlet_losses(
    spectra,
    pipe_length,
    temperature,
    pressure,
    sampleflow):

    if ((spectra is None) or 
        (temperature is None) or 
        (pressure is None) or 
        (sampleflow is None)):
        return None
    else:
        throughput = af.tubeloss(
            pd.Series(DP_STANDARD),
            sampleflow,
            pipe_length,
            temperature,
            pressure)
        spectra = spectra / throughput.values
        return spectra

def wagner_ion_mode_correction(spectra):
    if spectra is None:
        return None
    else:
        roberts_corr = 0.713*DP_STANDARD_NM**0.120
        spectra = spectra / roberts_corr
        
        return spectra

def dilution_correction(
    spectra,
    dilution_flow,
    sampleflow):

    if ((spectra is None) or
        (dilution_flow is None) or
        (sampleflow is None)):
        return None
    else:
        dilution_factor = sampleflow.values/dilution_flow.values
        spectra  = spectra * dilution_factor.reshape(-1,1)
        return spectra
    
def flags2polarity(
    flags_spectra,
    flags_offset,
    flag_explanations):

    if ((flags_spectra is None) or 
        (flags_offset is None) or
        (flag_explanations is None)):
        
        return None, None
    
    else:
        
        LEN_FLAG = len(flag_explanations["flag"].values)
        LEN_TIME = len(flags_spectra.index)
        
        flags_pos_spectra = pd.DataFrame(index = flags_spectra.index, 
                                         columns = flag_explanations["message"].values,
                                         data = np.zeros((LEN_TIME,LEN_FLAG)),dtype=int)
        flags_neg_spectra = pd.DataFrame(index = flags_spectra.index, 
                                         columns = flag_explanations["message"].values,
                                         data = np.zeros((LEN_TIME,LEN_FLAG)),dtype=int)
        
        for flag, message in zip(flag_explanations["flag"],flag_explanations["message"]):
                        
            spectra_idx = flags_spectra[(flags_spectra==flag).sum(axis=1)>0].index
            offset_idx = flags_offset[(flags_offset==flag).sum(axis=1)>0].index
            combined_idx = spectra_idx.join(offset_idx,how="outer")
                        
            # Message/flag only concerns positive polarity of the given spectrum
            if ("+" in message):                
                flags_pos_spectra.loc[combined_idx,message] = 1 
                                    
            # M essage/flag only concerns negative polarity
            elif ("−" in message):   
                flags_neg_spectra.loc[combined_idx,message] = 1
            
            # Message/flag concerns both polarities
            else:
                flags_pos_spectra.loc[combined_idx,message] = 1
                flags_neg_spectra.loc[combined_idx,message] = 1
            
        return flags_neg_spectra, flags_pos_spectra

def remove_corona_ions(spectra):

    if (spectra is None):
        return None

    # Only consider likely limit range
    lower = 1.5e-9
    upper = 5.0e-9
    c = (lower <= spectra.columns.values) & (upper >= spectra.columns.values)
    spectra2 = spectra.loc[:, c]
    
    # If all the values are NA return spectra
    # else ifnd the cutoff and set data below it to NA

    if spectra2.isna().all().all():
        return spectra
    else:
        corona_lim = spectra2.columns.values[spectra2.median().diff().abs().argmax()]
    
        # Set values below corona ion limit to NaNs
        spectra.iloc[:,spectra.columns.values<=corona_lim]=np.nan

        return spectra

def save_as_netcdf(
    netcdf_save_path,
    negion_data,
    posion_data,
    negpar_data,
    pospar_data,
    negion_flags,
    posion_flags,
    negpar_flags,
    pospar_flags ,
    flag_explanations,
    temp_data,
    pres_data,
    rh_data,
    sampleflow_data,
    measurement_info):
    
    # If eveything is None, then save nothing
    if ((negion_data is None) &
        (posion_data is None) &
        (negpar_data is None) &
        (pospar_data is None)):
        return False
    else:    
        # Something is not None so let's extract the time and
        # generate the coords
        for spectra in [negion_data,posion_data,negpar_data,pospar_data]:
            if spectra is not None:
                time = spectra.index.values
                LEN_TIME = len(time)
                nan_data = np.nan*np.ones((LEN_TIME,LEN_DP))
                nan_env_data = np.nan*np.ones(LEN_TIME)
                break

        ds = xr.Dataset()
        ds = ds.assign_coords(
            coords = {
                "time": time,
                "diameter": DP_STANDARD
            }
        )
        
        ds.diameter.attrs["units"] = "m"
        ds.diameter.attrs["description"] = "Bin geometric mean diameters"
        ds.time.attrs["timezone"] = "utc"
                
        if negion_data is not None:
            ds = ds.assign(neg_ions=(("time","diameter"),negion_data.values))
        else:
            ds = ds.assign(neg_ions=(("time","diameter"),nan_data))
            
        if posion_data is not None:
            ds = ds.assign(pos_ions=(("time","diameter"),posion_data.values))
        else:
            ds = ds.assign(pos_ions=(("time","diameter"),nan_data))
            
        if negpar_data is not None:
            ds = ds.assign(neg_particles=(("time","diameter"),negpar_data.values))
        else:
            ds = ds.assign(neg_particles=(("time","diameter"),nan_data))
            
        if pospar_data is not None:
            ds = ds.assign(pos_particles=(("time","diameter"),pospar_data.values))
        else:
            ds = ds.assign(pos_particles=(("time","diameter"),nan_data))
            
        ds.neg_ions.attrs["units"] = "cm-3"
        ds.neg_ions.attrs["description"] = "Negative ion number-size distribution (dN/dlogDp)"
        
        ds.pos_ions.attrs["units"] = "cm-3"
        ds.pos_ions.attrs["description"] = "Positive ion number-size distribution (dN/dlogDp)"
        
        ds.neg_particles.attrs["units"] = "cm-3"
        ds.neg_particles.attrs["description"] = "Negative particle number-size distribution (dN/dlogDp)"
        
        ds.pos_particles.attrs["units"] = "cm-3"
        ds.pos_particles.attrs["description"] = "Positive particle number-size distribution (dN/dlogDp)"

        # Add environmental sensor data
        if temp_data is not None:
            ds = ds.assign(temperature=(("time",), temp_data.values))   
        else:
            ds = ds.assign(temperature=(("time",), nan_env_data))

        ds.temperature.attrs["units"] = "K"
        ds.temperature.attrs["description"] = "Sample air temperature at the inlet"

        if pres_data is not None:
            ds = ds.assign(pressure=(("time",), pres_data.values))   
        else:
            ds = ds.assign(pressure=(("time",), nan_env_data))

        ds.pressure.attrs["units"] = "Pa"
        ds.pressure.attrs["description"] = "Sample air pressure at the inlet"

        if rh_data is not None:
            ds = ds.assign(relhum=(("time",), rh_data.values))   
        else:
            ds = ds.assign(relhum=(("time",), nan_env_data))

        ds.relhum.attrs["units"] = "%"
        ds.relhum.attrs["description"] = "Sample air relative humidity at the inlet"

        if sampleflow_data is not None:
            ds = ds.assign(sample_flow=(("time",), sampleflow_data.values))   
        else:
            ds = ds.assign(sample_flow=(("time",), nan_env_data))

        ds.sample_flow.attrs["units"] = "L/min"
        ds.sample_flow.attrs["description"] = "Flow rate at the outlet"

        # Add flag explanations if they exist
        if flag_explanations is not None:
            ds = ds.assign_coords(flag = flag_explanations["message"].values)
            nan_flags = np.zeros((LEN_TIME,ds.flag.size)).astype(int)
        else:
            ds = ds.assign_coords(flag = [""])
            nan_flags = np.zeros((LEN_TIME,1)).astype(int)

        # If flags exist add them     
        if ((negion_flags is not None) and (flag_explanations is not None)):
            ds = ds.assign(neg_ion_flags=(("time","flag"),negion_flags.values))
        else:
            ds = ds.assign(neg_ion_flags=(("time","flag"),nan_flags))
            
        if ((posion_flags is not None) and (flag_explanations is not None)):
            ds = ds.assign(pos_ion_flags=(("time","flag"),posion_flags.values))
        else:
            ds = ds.assign(pos_ion_flags=(("time","flag"),nan_flags))
            
        if ((negpar_flags is not None) and (flag_explanations is not None)):
            ds = ds.assign(neg_particle_flags=(("time","flag"),negpar_flags.values))
        else:
            ds = ds.assign(neg_particle_flags=(("time","flag"),nan_flags))
            
        if ((pospar_flags is not None) and (flag_explanations is not None)):
            ds = ds.assign(pos_particle_flags=(("time","flag"),pospar_flags.values))
        else:
            ds = ds.assign(pos_particle_flags=(("time","flag"),nan_flags))
            
        # Add measurement info
        ds = ds.assign_attrs(measurement_info)
        
        ds.to_netcdf(netcdf_save_path)
        
        return True
       
def nais_processor(config_file):
    """ Processes NAIS data
    
    Parameters
    ----------
    
    config_file : str
        full path to configuration file

    """

    # Load configuration
    with open(config_file,'r') as stream:
        config = yaml.safe_load(stream)
        
        # Read in the configuration file
        load_path = config['data_folder']
        ide = config["id"]
        save_path = config['processed_folder']
        start_date = config['start_date']
        database = config['database_file']
        location = config['measurement_location']
        description = config['description']
        instrument_model = config['instrument_model']
        longitude = config["longitude"]
        latitude = config["latitude"]
        end_date = config['end_date']
        allow_reprocess = config["allow_reprocess"]
        redo_database = config["redo_database"]
        pipelength = config['inlet_length']
        do_inlet_loss_correction = config['do_inlet_loss_correction']
        convert_to_standard_conditions = config['convert_to_standard_conditions']
        do_wagner_ion_mode_correction = config["do_wagner_ion_mode_correction"]
        remove_charger_ions = config["remove_corona_ions"]
        file_format = config["file_format"]
        resolution = config["resolution"]
        fill_temperature = config["fill_temperature"]
        fill_pressure = config["fill_pressure"]
        fill_flowrate = config["fill_flowrate"]
        dilution_on = config["dilution_on"]
        
    # Check the config file
    assert isinstance(start_date,date)
    assert ((end_date is None) or isinstance(end_date,date))
    assert os.path.exists(save_path)
    assert all([os.path.exists(x) for x in load_path])
    assert isinstance(allow_reprocess,bool)
    assert isinstance(remove_charger_ions,bool)
    assert isinstance(convert_to_standard_conditions,bool)
    assert isinstance(do_wagner_ion_mode_correction,bool)
    assert isinstance(do_inlet_loss_correction,bool)
    assert isinstance(pipelength,float)
    assert (isinstance(longitude,float) or (longitude is None))
    assert (isinstance(latitude,float) or (latitude is None))
    assert (isinstance(fill_temperature,float) or (fill_temperature is None))
    assert (isinstance(fill_pressure,float) or (fill_pressure is None))
    assert (isinstance(fill_flowrate,float) or (fill_flowrate is None))
    assert isinstance(dilution_on,bool)
    assert isinstance(redo_database,bool)
    pd.tseries.frequencies.to_offset(resolution)
    
    # Extract relevant info for metadata from the config
    measurement_info = {
        'measurement_location':location,
        'id':ide,
        'description':description,
        'instrument_model':instrument_model,
        'start_date': start_date.strftime("%Y-%m-%d"),
        'end_date': str(end_date) if end_date is None else end_date.strftime("%Y-%m-%d"),
        'longitude':str(longitude),
        'latitude':str(latitude),
        'inlet_length':pipelength,
        'do_inlet_loss_correction':str(do_inlet_loss_correction),
        'convert_to_standard_conditions':str(convert_to_standard_conditions),
        "do_wagner_ion_mode_correction":str(do_wagner_ion_mode_correction),
        "remove_corona_ions":str(remove_charger_ions),
        "fill_temperature":str(fill_temperature),
        "fill_pressure":str(fill_pressure),
        "fill_flowrate":str(fill_flowrate),
        "dilution_on":str(dilution_on),
        "resolution":resolution,
        "nais_processor_version":version,
        "date_processed":date.today().strftime("%Y-%m-%d")
    }    

    end_date = date.today() if end_date is None else end_date

    if redo_database:
        try:
            os.remove(database)
        except:
            pass

    db = TinyDB(database)
    check = Query()

    filename_formats = [
    ["%Y-%m-%d.ions.nds","%Y-%m-%d.particles.nds","%Y-%m-%d.log"],
    [f"%Y%m%d-{file_format}-ions.spectra",f"%Y%m%d-{file_format}-particles.spectra",f"%Y%m%d-{file_format}.records"],
    [f"%Y%m%d-{file_format}-ions.spectra",f"%Y%m%d-{file_format}-particles.spectra",f"%Y%m%d-{file_format}.diagnostics"]]

    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    # List existing dates based on if diagnostic file was found
    list_of_existing_dates = [x["timestamp"] for x in db.search(check.diagnostics.exists())]

    # Make a list of all dates in the range
    list_of_datetimes = pd.date_range(start=start_date_str, end=end_date_str)

    preferred_path = load_path[0]
    preferred_format = filename_formats[0]

    # Add unprocessed datafiles to the database
    for x in list_of_datetimes:
        if (x.strftime("%Y%m%d") in list_of_existing_dates):
            continue
        else:
            files_found=False
            for i in range(len(load_path)):
                for j in range(len(filename_formats)):

                    ion_filename = os.path.join(load_path[i],x.strftime(filename_formats[j][0]))
                    particle_filename = os.path.join(load_path[i],x.strftime(filename_formats[j][1]))
                    diagnostic_filename = os.path.join(load_path[i],x.strftime(filename_formats[j][2]))

                    if ( (os.path.exists(ion_filename) | # ions
                         os.path.exists(particle_filename)) & # particles
                         os.path.exists(diagnostic_filename) # diagnostics
                       ):

                        date_str = x.strftime("%Y%m%d")

                        db.insert(
                            {"timestamp":date_str,
                            "diagnostics":diagnostic_filename}
                            )

                        if os.path.exists(ion_filename):
                            db.update(
                                {"ions":ion_filename},
                                check.timestamp==date_str)

                        if os.path.exists(particle_filename):
                            db.update(
                                {"particles":particle_filename},
                                check.timestamp==date_str)

                        files_found=True

                        if (load_path[i]!=preferred_path):
                            preferred_path = load_path.pop(i)
                            load_path.insert(0,preferred_path)

                        if (filename_formats[j]!=preferred_format):
                            preferred_format = filename_formats.pop(j)
                            filename_formats.insert(0,preferred_format)

                        print(f"Added {x.strftime('%Y%m%d')} to database ({location}) ...")
                        
                        break

                if files_found:
                    break

    # From the database find the last day with processed data
    processed_days = db.search(check.processed_file.exists())

    if len(processed_days)!=0:
        last_day=np.max([datetime.strptime(x["timestamp"],"%Y%m%d") 
            for x in processed_days]).strftime("%Y%m%d")
    else:
        allow_reprocess=True

    if allow_reprocess:
        database_iterator = iter(db.search(
            ((check.diagnostics.exists()) &
            (check.ions.exists() | check.particles.exists()) &
            (~check.processed_file.exists()) &
            (check.timestamp>=start_date_str) &
            (check.timestamp<=end_date_str))))
    else:
        database_iterator = iter(db.search(
            (check.diagnostics.exists() &
            (check.ions.exists() | check.particles.exists()) &
            (check.timestamp>=last_day) &
            (check.timestamp>=start_date_str) &
            (check.timestamp<=end_date_str))))

    for x in database_iterator:

        print("Processing %s (%s) ..." % (x["timestamp"], location))

        ions_exist = bool(db.search(
            check.ions.exists() & (check.timestamp==x["timestamp"])))

        particles_exist = bool(db.search(
            check.particles.exists() & (check.timestamp==x["timestamp"])))
       
        (ion_records, 
            particle_records, 
            ion_flags, 
            particle_flags, 
            offset_flags, 
            flag_explanations) = read_raw(x["diagnostics"],"records",x["timestamp"],resolution)

        # ions
        if ions_exist:

            ions = read_raw(x["ions"],"spectra",x["timestamp"],resolution)

            negion_datamatrix, posion_datamatrix = raw2sum(ions,"ions")
            negion_flags, posion_flags = flags2polarity(ion_flags, offset_flags, flag_explanations)
           
            # Get diagnostic data for corrections and conversions
            #if (convert_to_standard_conditions or do_inlet_loss_correction or dilution_on):
            (temperature_ion,
                pressure_ion,
                sampleflow_ion,
                dilution_flow_ion,
                relhum_ion,
                temperature_ion_filled,
                pressure_ion_filled,
                sampleflow_ion_filled) = get_diagnostic_data(
                    ion_records,
                    fill_pressure,
                    fill_temperature,
                    fill_flowrate)

            if convert_to_standard_conditions:
                negion_datamatrix = bring_to_sealevel(
                    negion_datamatrix,
                    temperature_ion,
                    pressure_ion
                    )
                posion_datamatrix = bring_to_sealevel(
                    posion_datamatrix,
                    temperature_ion,
                    pressure_ion
                    )
             
            if do_inlet_loss_correction:
                negion_datamatrix = correct_inlet_losses(
                    negion_datamatrix,
                    pipelength,
                    temperature_ion,
                    pressure_ion,
                    sampleflow_ion
                    )
                posion_datamatrix = correct_inlet_losses(
                    posion_datamatrix,
                    pipelength,
                    temperature_ion,
                    pressure_ion,
                    sampleflow_ion
                    )
            
            if do_wagner_ion_mode_correction:
                negion_datamatrix = wagner_ion_mode_correction(negion_datamatrix)
                posion_datamatrix = wagner_ion_mode_correction(posion_datamatrix)
                
            if dilution_on:
                negion_datamatrix = dilution_correction(negion_datamatrix,dilution_flow_ion,sampleflow_ion)
                posion_datamatrix = dilution_correction(posion_datamatrix,dilution_flow_ion,sampleflow_ion)
        
        else:
            temperature_ion = None
            pressure_ion = None
            sampleflow_ion = None
            relhum_ion = None
            negion_datamatrix, posion_datamatrix = None, None
            negion_flags, posion_flags = None, None
                    
        # Particles
        if particles_exist:

            particles = read_raw(x["particles"],"spectra",x["timestamp"],resolution)            
            negpar_datamatrix, pospar_datamatrix = raw2sum(particles,"particles")
            negpar_flags, pospar_flags = flags2polarity(particle_flags, offset_flags, flag_explanations)

            ## Get diagnostic data for corrections and conversions
            #if (convert_to_standard_conditions or do_inlet_loss_correction or dilution_on):
            (temperature_particle,
                pressure_particle,
                sampleflow_particle,
                dilution_flow_particle,
                relhum_particle,
                temperature_particle_filled,
                pressure_particle_filled,
                sampleflow_particle_filled) = get_diagnostic_data(
                    particle_records,
                    fill_pressure,
                    fill_temperature,
                    fill_flowrate)

            if convert_to_standard_conditions:
                negpar_datamatrix = bring_to_sealevel(
                    negpar_datamatrix,
                    temperature_particle,
                    pressure_particle)
                pospar_datamatrix = bring_to_sealevel(
                    pospar_datamatrix,
                    temperature_particle,
                    pressure_particle)
                
            if do_inlet_loss_correction:
                negpar_datamatrix = correct_inlet_losses(
                    negpar_datamatrix,
                    pipelength,
                    temperature_particle,
                    pressure_particle,
                    sampleflow_particle)
                pospar_datamatrix = correct_inlet_losses(
                    pospar_datamatrix,
                    pipelength,
                    temperature_particle,
                    pressure_particle,
                    sampleflow_particle)    
                
            if remove_charger_ions:  
                negpar_datamatrix = remove_corona_ions(negpar_datamatrix)
                pospar_datamatrix = remove_corona_ions(pospar_datamatrix)

            if dilution_on:
                negpar_datamatrix = dilution_correction(negpar_datamatrix,dilution_flow_particle,sampleflow_particle)
                pospar_datamatrix = dilution_correction(pospar_datamatrix,dilution_flow_particle,sampleflow_particle)

        else:
            temperature_particle = None
            pressure_particle = None
            sampleflow_particle = None
            relhum_particle = None
            negpar_datamatrix, pospar_datamatrix = None, None
            negpar_flags, pospar_flags = None, None


        # Both ion and particle data exists
        if ((temperature_ion is not None) & (temperature_particle is not None)):
            if ((not temperature_ion_filled) & (not temperature_particle_filled)): 
                temperature_data = (temperature_ion + temperature_particle)/2.
            elif (temperature_ion_filled & (not temperature_particle_filled)):
                temperature_data = temperature_particle
            elif ((not temperature_ion_filled) & temperature_particle_filled):
                temperature_data = temperature_ion
            else:
                temperature_data = temperature_ion
        # Only ion data exists
        elif (temperature_ion is not None):
            temperature_data = temperature_ion
        # Only particle data exists
        elif (temperature_particle is not None):
            temperature_data = temperature_particle
        else:
            temperature_data = None


        # Both ion and particle data exists
        if ((pressure_ion is not None) & (pressure_particle is not None)):
            if ((not pressure_ion_filled) & (not pressure_particle_filled)): 
                pressure_data = (pressure_ion + pressure_particle)/2.
            elif (pressure_ion_filled & (not pressure_particle_filled)):
                pressure_data = pressure_particle
            elif ((not pressure_ion_filled) & pressure_particle_filled):
                pressure_data = pressure_ion
            else:
                pressure_data = pressure_ion
        # Only ion data exists
        elif (pressure_ion is not None):
            pressure_data = pressure_ion
        # Only particle data exists
        elif (pressure_particle is not None):
            pressure_data = pressure_particle
        else:
            pressure_data = None

        # Both ion and particle data exists
        if ((relhum_ion is not None) & (relhum_particle is not None)):
            relhum_data = (relhum_ion + relhum_particle)/2.
        # Only ion data exists
        elif (relhum_ion is not None):
            relhum_data = relhum_ion
        # Only particle data exists
        elif (relhum_particle is not None):
            relhum_data = relhum_particle
        else:
            relhum_data = None


        # Both ion and particle data exists
        if ((sampleflow_ion is not None) & (sampleflow_particle is not None)):
            if ((not sampleflow_ion_filled) & (not sampleflow_particle_filled)): 
                sampleflow_data = (sampleflow_ion + sampleflow_particle)/2.
            elif (sampleflow_ion_filled & (not sampleflow_particle_filled)):
                sampleflow_data = sampleflow_particle
            elif ((not sampleflow_ion_filled) & sampleflow_particle_filled):
                sampleflow_data = sampleflow_ion
            else:
                sampleflow_data = sampleflow_ion
        # Only ion data exists
        elif (sampleflow_ion is not None):
            sampleflow_data = sampleflow_ion
        # Only particle data exists
        elif (sampleflow_particle is not None):
            sampleflow_data = sampleflow_particle
        else:
            sampleflow_data = None

        my_save_path = os.path.join(save_path,"NAIS_"+x["timestamp"]+".nc")
        
        saved = save_as_netcdf(
            my_save_path,
            negion_datamatrix,
            posion_datamatrix,
            negpar_datamatrix,
            pospar_datamatrix,
            negion_flags,
            posion_flags,
            negpar_flags,
            pospar_flags,
            flag_explanations,
            temperature_data,
            pressure_data,
            relhum_data,
            sampleflow_data,
            measurement_info
        )
        
        if saved:
            db.update({"processed_file": my_save_path},check.timestamp==x["timestamp"])
            
    print("Done!")
