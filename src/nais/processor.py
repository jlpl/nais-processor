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
LEN_DP = 40
LEN_TIME = 289

DP_STANDARD = np.array([
5.00000000e-10, 5.61009227e-10, 6.29462706e-10, 7.06268772e-10,
7.92446596e-10, 8.89139705e-10, 9.97631157e-10, 1.11936057e-09,
1.25594322e-09, 1.40919147e-09, 1.58113883e-09, 1.77406695e-09,
1.99053585e-09, 2.23341796e-09, 2.50593617e-09, 2.81170663e-09,
3.15478672e-09, 3.53972892e-09, 3.97164117e-09, 4.45625469e-09,
5.00000000e-09, 5.61009227e-09, 6.29462706e-09, 7.06268772e-09,
7.92446596e-09, 8.89139705e-09, 9.97631157e-09, 1.11936057e-08,
1.25594322e-08, 1.40919147e-08, 1.58113883e-08, 1.77406695e-08,
1.99053585e-08, 2.23341796e-08, 2.50593617e-08, 2.81170663e-08,
3.15478672e-08, 3.53972892e-08, 3.97164117e-08, 4.45625469e-08])

DP_STANDARD_NM = 1e9*np.array([
5.00000000e-10, 5.61009227e-10, 6.29462706e-10, 7.06268772e-10,
7.92446596e-10, 8.89139705e-10, 9.97631157e-10, 1.11936057e-09,
1.25594322e-09, 1.40919147e-09, 1.58113883e-09, 1.77406695e-09,
1.99053585e-09, 2.23341796e-09, 2.50593617e-09, 2.81170663e-09,
3.15478672e-09, 3.53972892e-09, 3.97164117e-09, 4.45625469e-09,
5.00000000e-09, 5.61009227e-09, 6.29462706e-09, 7.06268772e-09,
7.92446596e-09, 8.89139705e-09, 9.97631157e-09, 1.11936057e-08,
1.25594322e-08, 1.40919147e-08, 1.58113883e-08, 1.77406695e-08,
1.99053585e-08, 2.23341796e-08, 2.50593617e-08, 2.81170663e-08,
3.15478672e-08, 3.53972892e-08, 3.97164117e-08, 4.45625469e-08])

MOB_STANDARD = np.array([
7.81537869e+00, 6.20905133e+00, 4.93298365e+00, 3.91926324e+00,
3.11394395e+00, 2.47417363e+00, 1.96591248e+00, 1.56212080e+00,
1.24131927e+00, 9.86445511e-01, 7.83945709e-01, 6.23053061e-01,
4.95214643e-01, 3.93636140e-01, 3.12920199e-01, 2.48779161e-01,
1.97806874e-01, 1.57297423e-01, 1.25101151e-01, 9.95102784e-02,
7.91680500e-02, 6.29965615e-02, 5.01394299e-02, 3.99162509e-02,
3.17864225e-02, 2.53204062e-02, 2.01768979e-02, 1.60846922e-02,
1.28282759e-02, 1.02363833e-02, 8.17290349e-03, 6.52965824e-03,
5.22066324e-03, 4.17756954e-03, 3.34604174e-03, 2.68288051e-03,
2.15373646e-03, 1.73129367e-03, 1.39382651e-03, 1.12405285e-03])

DLOGDP_STANDARD = np.array([
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

DLOGMOB_STANDARD = np.array([
0.09992478, 0.09992478, 0.09991557, 0.09990524, 0.09989364,
0.09988062, 0.09986599, 0.09984956, 0.09983112, 0.0998104 ,
0.09978711, 0.09976095, 0.09973156, 0.09969851, 0.09966137,
0.0996196 , 0.09957262, 0.09951977, 0.0994603 , 0.09939336,
0.09931799, 0.0992331 , 0.09913745, 0.09902965, 0.09890809,
0.09877098, 0.09861627, 0.09844164, 0.09824446, 0.09802175,
0.09777017, 0.09748592, 0.09716477, 0.09680199, 0.09639233,
0.09593002, 0.09540877, 0.09482182, 0.09416202, 0.09342199])

FILENAME_FORMATS = [
["%Y-%m-%d.ions.nds","%Y-%m-%d.particles.nds","%Y-%m-%d.log"],
["%Y%m%d-block-ions.spectra","%Y%m%d-block-particles.spectra","%Y%m%d-block.records"],
["%Y%m%d-block-ions.spectra","%Y%m%d-block-particles.spectra","%Y%m%d-block.diagnostics"]]

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

DILUTION_FLOW_NAMES = [
"diluter_sample_flow_rate.mean"]

# Standard conditions
TEMP_REF = 273.15
PRES_REF = 101325.0

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
        f.write("longitude: # decimal degrees west/east = -/+ (float)\n")
        f.write("latitude: # decimal degrees south/north = -/+ (float)\n")
        f.write("data_folder: # Full paths to raw data folders\n")
        f.write("- # Data folder 1\n")
        f.write("- # Data folder 2, and so on...\n")
        f.write("processed_folder: # Full path to folder where procesed data is saved\n")
        f.write("database_file: # Full path to database file (will be created on first run) \n")
        f.write("start_date: # Format: yyyy-mm-dd\n")
        f.write("end_date: # Format: yyyy-mm-dd or '' for current day\n")
        f.write("inlet_length: # length of inlet in meters (float)\n")
        f.write("do_inlet_loss_correction: # true or false\n")
        f.write("convert_to_standard_conditions: # true or false\n")
        f.write("do_wagner_ion_mode_correction: # true or false\n")
        f.write("remove_corona_ions: # true or false\n")
        f.write("allow_reprocess: # true or false\n")
        f.write("use_fill_values: # true or false\n")
        f.write("fill_temperature: # temperature in K (float)\n")
        f.write("fill_pressure: # pressure in Pa (float)\n")
        f.write("fill_flowrate: # flow rate in lpm (float)\n")
        f.write("dilution_on: # true or false (is the integrated dilution system used)\n")

def read_raw(file_name,file_type,timestamp):
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
        
        # Construct time index, this will fail for example in the rare case if time zone was changed
        try:
            begin_time = pd.DatetimeIndex(df[df.columns[0]])
            end_time = pd.DatetimeIndex(df[df.columns[1]])
            center_time = begin_time + (end_time - begin_time)/2.
            df.index = center_time
            df = df.sort_index() # sort just in case
            df = df[~df.index.duplicated(keep='first')] # remove duplicates just in case
            
            if df.index[0].tz is None:
                df.index = df.index.tz_localize("UTC")
            else:
                pass
        
        except:
            if file_type=="records":
                return None,None,None,None,None,None
            if file_type=="spectra":
                return None
        
        # Still if something is wrong, just bail out 
        if not df.index.is_monotonic_increasing:
            if file_type=="records":
                return None,None,None,None,None,None
            if file_type=="spectra":
                return None
        
        # Define a standard time index for all data of the day
        standard_start_time = pd.to_datetime(timestamp)
        standard_end_time = standard_start_time + pd.Timedelta(days=1)
        standard_time = pd.date_range(start=standard_start_time, end=standard_end_time, freq="5min").tz_localize(df.index[0].tz)
        
        assert (len(standard_time)==LEN_TIME)
        
        if file_type=="records":
            # Extract records for ions, particles and offset
            ion_records_and_flags = df[df.opmode=='ions']
            particle_records_and_flags = df[df.opmode=='particles']
            offset_records_and_flags = df[df.opmode=='offset']
            
            if ion_records_and_flags.empty==False:
                ion_records_and_flags = ion_records_and_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta("5min"))
                ion_records_and_flags.index = ion_records_and_flags.index.tz_convert('UTC')
                ion_records = ion_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
                ion_flags = ion_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')
            else:
                ion_records = None
                ion_flags = None
                
            if particle_records_and_flags.empty == False:
                particle_records_and_flags = particle_records_and_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta("5min"))
                particle_records_and_flags.index = particle_records_and_flags.index.tz_convert('UTC')
                particle_records = particle_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
                particle_flags = particle_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')
            else:
                particle_records = None
                particle_flags = None
                
            if offset_records_and_flags.empty==False:
                offset_records_and_flags = offset_records_and_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta("5min"))
                offset_records_and_flags.index = offset_records_and_flags.index.tz_convert('UTC')
                offset_records = offset_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
                offset_flags = offset_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')
            else:
                offset_records = None
                offset_flags = None
                
            if flag_explanations.empty==True:
                flag_explanations = None

            return ion_records, particle_records, ion_flags, particle_flags, offset_flags, flag_explanations
        
        if file_type=="spectra":
            spectra = df.iloc[:,3:].apply(pd.to_numeric, errors='coerce').astype(float)
            spectra = spectra.reindex(standard_time, method="nearest", tolerance=pd.Timedelta("5min"))
            spectra.index = spectra.index.tz_convert('UTC')
            
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
    dilution_flow_name=None

    for temp_name in TEMPERATURE_NAMES:
        if temp_name in diag_params:
            temperature_name = temp_name
            break

    for pres_name in PRESSURE_NAMES:
        if pres_name in diag_params:
            pressure_name = pres_name
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

    return temperature_name, pressure_name, sampleflow_name, pos_sampleflow_name, neg_sampleflow_name, dilution_flow_name

def get_diagnostic_data(
    records,
    use_fill_values,
    fill_pressure,
    fill_temperature,
    fill_flowrate):

    if records is None:
        return None, None, None, None

    else:        
        # Check that the relevant diagnostic data is found
        (temperature_name,
        pressure_name,
        sampleflow_name,
        pos_sampleflow_name,
        neg_sampleflow_name,
        dilution_flow_name) = find_diagnostic_names(list(records))

        if temperature_name is not None:
            temperature = 273.15 + records[temperature_name].astype(float)
            # Values may be missing: e.g. sensor is broken
            if (temperature.isna().all() and use_fill_values):
                temperature = pd.Series(index = records.index, dtype=float)
                temperature[:] = fill_temperature
        elif use_fill_values:
            temperature = pd.Series(index = records.index, dtype=float)
            temperature[:] = fill_temperature
        else:
            temperature = None

        if pressure_name is not None:
            pressure = 100.0 * records[pressure_name].astype(float)
            if (pressure.isna().all() and use_fill_values):
                pressure = pd.Series(index = pressure.index, dtype=float)
                pressure[:] = fill_pressure
        elif use_fill_values:
            pressure = pd.Series(index = records.index, dtype=float)
            pressure[:] = fill_pressure
        else:
            pressure = None

        if sampleflow_name is not None:
            sampleflow = records[sampleflow_name].astype(float)
            if (sampleflow.isna().all() and use_fill_values):
                sampleflow = pd.Series(index = records.index, dtype=float)
                sampleflow[:] = fill_flowrate
        elif ((neg_sampleflow_name is not None) and (pos_sampleflow_name is not None)):
            neg_sampleflow = records[neg_sampleflow_name].astype(float)
            pos_sampleflow = records[pos_sampleflow_name].astype(float)
            sampleflow = neg_sampleflow + pos_sampleflow
            if (sampleflow.isna().all() and use_fill_values):
                sampleflow = pd.Series(index = records.index, dtype=float)
                sampleflow[:] = fill_flowrate
        elif use_fill_values:
            sampleflow = pd.Series(index = records.index, dtype=float)
            sampleflow[:] = fill_flowrate
        else:
            sampleflow = None

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
        
        return temperature, pressure, sampleflow, dilution_flow

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
            DP_STANDARD,
            sampleflow.values*1.667e-5,
            pipe_length,
            temperature.values,
            pressure.values)
        spectra = spectra / throughput
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
    
    # Find maximum difference between size bin medians
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
                nan_data = np.nan*np.ones((LEN_TIME,LEN_DP))
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
        save_path = config['processed_folder']
        start_date = config['start_date']
        database = config['database_file']
        location = config['measurement_location']
        longitude = config["longitude"]
        latitude = config["latitude"]
        end_date = config['end_date']
        allow_reprocess = config["allow_reprocess"]
        pipelength = config['inlet_length']
        do_inlet_loss_correction = config['do_inlet_loss_correction']
        convert_to_standard_conditions = config['convert_to_standard_conditions']
        do_wagner_ion_mode_correction = config["do_wagner_ion_mode_correction"]
        remove_charger_ions = config["remove_corona_ions"]
        use_fill_values = config["use_fill_values"]
        if use_fill_values:
            fill_temperature = config["fill_temperature"]
            fill_pressure = config["fill_pressure"]
            fill_flowrate = config["fill_flowrate"]
        else:
            fill_temperature = 1.0
            fill_pressure = 1.0
            fill_flowrate = 1.0
        dilution_on = config["dilution_on"]
        
    # Check the config file
    assert isinstance(start_date,date)
    assert (end_date=='' or isinstance(end_date,date))
    assert os.path.exists(save_path)
    assert all([os.path.exists(x) for x in load_path])
    assert isinstance(allow_reprocess,bool)
    assert isinstance(remove_charger_ions,bool)
    assert isinstance(convert_to_standard_conditions,bool)
    assert isinstance(do_wagner_ion_mode_correction,bool)
    assert isinstance(do_inlet_loss_correction,bool)
    assert isinstance(pipelength,float)
    assert isinstance(use_fill_values,bool)
    assert isinstance(longitude,float)
    assert isinstance(latitude,float)
    assert isinstance(fill_temperature,float)
    assert isinstance(fill_pressure,float)
    assert isinstance(fill_flowrate,float)
    assert isinstance(dilution_on,bool)
    
    # Extract relevant info for metadata from the config
    measurement_info = {
        'measurement_location':location,
        'longitude':longitude,
        'latitude':latitude,
        'inlet_length':pipelength,
        'do_inlet_loss_correction':str(do_inlet_loss_correction),
        'convert_to_standard_conditions':str(convert_to_standard_conditions),
        "do_wagner_ion_mode_correction":str(do_wagner_ion_mode_correction),
        "remove_corona_ions":str(remove_charger_ions),
        "use_fill_values":str(use_fill_values),
        "fill_temperature":fill_temperature,
        "fill_pressure":fill_pressure,
        "fill_flowrate":fill_flowrate,
        "dilution_on":str(dilution_on),
    }    

    end_date = date.today() if end_date=='' else end_date

    db = TinyDB(database)
    check = Query()

    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    # List existing dates based on if diagnostic file was found
    list_of_existing_dates = [x["timestamp"] for x in db.search(check.diagnostics.exists())]

    if len(list_of_existing_dates)==0:
        print("Building database...")
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
                for y in FILENAME_FORMATS:

                    ion_filename = os.path.join(z,x.strftime(y[0]))
                    particle_filename = os.path.join(z,x.strftime(y[1]))
                    diagnostic_filename = os.path.join(z,x.strftime(y[2]))

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
                        break

                if files_found:
                    break

    # From the database find the last day with processed data
    processed_days = db.search(check.processed_file.exists())

    if len(processed_days)!=0:
        last_day=np.max([datetime.strptime(x["timestamp"],"%Y%m%d") 
            for x in processed_days]).strftime("%Y%m%d")
    else:
        last_day=end_date_str

    if allow_reprocess:
        database_iterator = iter(db.search(
            (check.diagnostics.exists() &
            (check.ions.exists() | check.particles.exists()) &
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

        print("Processing %s (%s)" % (x["timestamp"], location))

        ions_exist = bool(db.search(
            check.ions.exists() & (check.timestamp==x["timestamp"])))

        particles_exist = bool(db.search(
            check.particles.exists() & (check.timestamp==x["timestamp"])))
       
        (ion_records, 
            particle_records, 
            ion_flags, 
            particle_flags, 
            offset_flags, 
            flag_explanations) = read_raw(x["diagnostics"],"records",x["timestamp"])

        # ions
        if ions_exist:

            ions = read_raw(x["ions"],"spectra",x["timestamp"])                        
            negion_datamatrix, posion_datamatrix = raw2sum(ions,"ions")
            negion_flags, posion_flags = flags2polarity(ion_flags, offset_flags, flag_explanations)
                        
            # Get diagnostic data for corrections and conversions
            if (convert_to_standard_conditions or do_inlet_loss_correction or dilution_on):
                temperature_ion,pressure_ion,sampleflow_ion,dilution_flow_ion = get_diagnostic_data(
                    ion_records,
                    use_fill_values,
                    fill_pressure,
                    fill_temperature,
                    fill_flowrate)

            if convert_to_standard_conditions:
                negion_datamatrix = bring_to_sealevel(
                    negion_datamatrix,
                    temperature_ion,
                    pressure_ion)
                posion_datamatrix = bring_to_sealevel(
                    posion_datamatrix,
                    temperature_ion,
                    pressure_ion)
                
            if do_inlet_loss_correction:
                negion_datamatrix = correct_inlet_losses(
                    negion_datamatrix,
                    pipelength,
                    temperature_ion,
                    pressure_ion,
                    sampleflow_ion)
                posion_datamatrix = correct_inlet_losses(
                    posion_datamatrix,
                    pipelength,
                    temperature_ion,
                    pressure_ion,
                    sampleflow_ion)
                
            if do_wagner_ion_mode_correction:
                negion_datamatrix = wagner_ion_mode_correction(negion_datamatrix)
                posion_datamatrix = wagner_ion_mode_correction(posion_datamatrix)
                
            if dilution_on:
                negion_datamatrix = dilution_correction(negion_datamatrix,dilution_flow_ion,sampleflow_ion)
                posion_datamatrix = dilution_correction(posion_datamatrix,dilution_flow_ion,sampleflow_ion)

        else:
            negion_datamatrix, posion_datamatrix = None, None
            negion_flags, posion_flags = None, None
                    
        # Particles
        if particles_exist:

            particles = read_raw(x["particles"],"spectra",x["timestamp"])            
            negpar_datamatrix, pospar_datamatrix = raw2sum(particles,"particles")
            negpar_flags, pospar_flags = flags2polarity(particle_flags, offset_flags, flag_explanations)

            # Get diagnostic data for corrections and conversions
            if (convert_to_standard_conditions or do_inlet_loss_correction or dilution_on):
                temperature_particle,pressure_particle,sampleflow_particle,dilution_flow_particle = get_diagnostic_data(
                    particle_records,
                    use_fill_values,
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
            negpar_datamatrix, pospar_datamatrix = None, None
            negpar_flags, pospar_flags = None, None
                
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
            measurement_info
        )
        
        if saved:
            db.update({"processed_file": my_save_path},check.timestamp==x["timestamp"])
            
    print("Done!")
