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

# Names of relevant quantities found in the CIC files
POS_SAMPLEFLOW_NAMES = ["red_flow_rate.mean"]
NEG_SAMPLEFLOW_NAMES = ["blue_flow_rate.mean"]

POS_TEMPERATURE_NAMES = ["red_temperature.mean"]           
NEG_TEMPERATURE_NAMES = ["blue_temperature.mean"]

POS_PRESSURE_NAMES = ["red_pressure.mean"]     
NEG_PRESSURE_NAMES = ["blue_pressure.mean"]

NEG_CLUSTER_ION_MOB_NAMES = ["neg_avg_cluster_mob"]

NEG_MOB_1_NAMES = ["neg_mob_1"]
NEG_MOB_2_NAMES = ["neg_mob_2"]
NEG_MOB_3_NAMES = ["neg_mob_3"]

NEG_CONC_1_NAMES = ["neg_conc_1"]
NEG_CONC_2_NAMES = ["neg_conc_2"]
NEG_CONC_3_NAMES = ["neg_conc_3"]

POS_CLUSTER_ION_MOB_NAMES = ["pos_avg_cluster_mob"]

POS_MOB_1_NAMES = ["pos_mob_1"]
POS_MOB_2_NAMES = ["pos_mob_2"]
POS_MOB_3_NAMES = ["pos_mob_3"]

POS_CONC_1_NAMES = ["pos_conc_1"]
POS_CONC_2_NAMES = ["pos_conc_2"]
POS_CONC_3_NAMES = ["pos_conc_3"]

NAME_LIST = [
NEG_SAMPLEFLOW_NAMES,
POS_SAMPLEFLOW_NAMES,
NEG_TEMPERATURE_NAMES,
POS_TEMPERATURE_NAMES,
NEG_PRESSURE_NAMES,
POS_PRESSURE_NAMES,
NEG_CLUSTER_ION_MOB_NAMES,
POS_CLUSTER_ION_MOB_NAMES,
NEG_CONC_1_NAMES,
POS_CONC_1_NAMES,
NEG_CONC_2_NAMES,
POS_CONC_2_NAMES,
NEG_CONC_3_NAMES,
POS_MOB_3_NAMES,
NEG_MOB_1_NAMES,
POS_MOB_1_NAMES,
NEG_MOB_2_NAMES,
POS_MOB_2_NAMES,
NEG_MOB_3_NAMES,
POS_MOB_3_NAMES
]

# Standard conditions
TEMP_REF = 293.15
PRES_REF = 101325.0

def make_cic_config_template(file_name):
    """  
    Make a configuration file template

    Parameters
    ----------

    file_name : str
        full path to configuration file. For example `/home/user/cic_config.yml`

    """
    
    with open(file_name,"w") as f:
        f.write("measurement_location: # Name of the measurement site\n")
        f.write("description: # Additional description of the measurement\n")
        f.write("instrument_model: # e.g. CIC-2-1")
        f.write("longitude: # decimal degrees west/east = -/+ (float)\n")
        f.write("latitude: # decimal degrees south/north = -/+ (float)\n")
        f.write("data_folder: # Full paths to raw data folders\n")
        f.write("- # Data folder 1\n")
        f.write("- # Data folder 2, and so on...\n")
        f.write("processed_folder: # Full path to folder where procesed data is saved\n")
        f.write("database_file: # Full path to database file (will be created on first run)\n")
        f.write("start_date: # Format: yyyy-mm-dd\n")
        f.write("end_date: # Format: yyyy-mm-dd or '' for current day\n")
        f.write("inlet_length: # length of inlet in meters (float)\n")
        f.write("do_inlet_loss_correction: # true or false\n")
        f.write("convert_to_standard_conditions: # true or false\n")
        f.write("allow_reprocess: # true or false\n")
        f.write("redo_database: # true or false\n")
        f.write('file_format: # 1s, 10s or block\n')
        f.write('resolution: # processed data time resolution (pandas time offset string), e.g. 5min')

def read_raw(file_name,timestamp,resolution_str):
    with open(file_name,'r') as f:
        header_found = False
        data_matrix = []
        flag_explanations = []
        lines = f.read().splitlines()
        
        for line in lines:
            # Skip empty
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
        return None,None,None,None

    else:
        # Construct dataframes
        df = pd.DataFrame(columns = header, data = data_matrix, dtype = str)
        flag_explanations = pd.DataFrame(columns=["flag","message"], data=flag_explanations, dtype=str)
        
        # Remove duplicate flags (may happen due to restarts)        
        flag_explanations = flag_explanations[~flag_explanations["message"].duplicated()]
       
        # Contruct the datetime index
        # there can be timezone change

        begin_time = []
        end_time = []

        for bt, et in zip(df[df.columns[0]], df[df.columns[1]]):
            # Computer in utc time
            data_tz = pd.to_datetime(bt).tz
            if data_tz is None:
                begin_time.append(pd.to_datetime(bt).tz_localize("UTC"))
                end_time.append(pd.to_datetime(et).tz_localize("UTC"))
            # Computer not in utc time
            else:
                begin_time.append(pd.to_datetime(bt).tz_convert("UTC"))
                end_time.append(pd.to_datetime(et).tz_convert("UTC"))

        begin_time = pd.DatetimeIndex(begin_time)
        end_time = pd.DatetimeIndex(end_time)
        center_time = begin_time + (end_time - begin_time)/2.
        df.index = center_time

        df = df.sort_index() # sort just in case
        df = df[~df.index.duplicated(keep='first')] # remove duplicates just in case
        
        # Still if something is wrong, just bail out 
        if not df.index.is_monotonic_increasing:
            return None,None,None,None
        
        # Define a standard time index for all data of the day
        if data_tz is None:
            standard_start_time = pd.to_datetime(timestamp).tz_localize('UTC')
        else:
            standard_start_time = pd.to_datetime(timestamp).tz_localize(data_tz).tz_convert("UTC")
        standard_end_time = standard_start_time + pd.Timedelta(days=1)
        standard_time = pd.date_range(start=standard_start_time, end=standard_end_time, freq=resolution_str)
        
        # Extract records for ions and offset
        ion_records_and_flags = df[df.opmode=='ions']
        offset_records_and_flags = df[df.opmode=='offset']
        
        if ion_records_and_flags.empty==False:
            ion_records_and_flags = ion_records_and_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))
            ion_records = ion_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
            ion_flags = ion_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')
        else:
            ion_records = None
            ion_flags = None

        if offset_records_and_flags.empty==False:
            offset_records_and_flags = offset_records_and_flags.reindex(standard_time, method="nearest", tolerance=pd.Timedelta(resolution_str))
            offset_records = offset_records_and_flags.iloc[:,3:-1].apply(pd.to_numeric, errors='coerce').astype(float)
            offset_flags = offset_records_and_flags.iloc[:,-1].fillna('').str.split("!",expand=True).fillna('')
        else:
            offset_records = None
            offset_flags = None

        if flag_explanations.empty==True:
            flag_explanations = None

        return ion_records, ion_flags, offset_flags, flag_explanations

def find_names(params):
    neg_sampleflow_name=None
    pos_sampleflow_name=None
    neg_temperature_name=None
    pos_temperature_name=None
    neg_pressure_name=None
    pos_pressure_name=None
    neg_cluster_mob_name=None
    pos_cluster_mob_name=None
    neg_conc_1_name=None
    neg_conc_2_name=None
    neg_conc_3_name=None
    pos_conc_1_name=None
    pos_conc_2_name=None
    pos_conc_3_name=None
    neg_mob_1_name=None
    neg_mob_2_name=None
    neg_mob_3_name=None
    pos_mob_1_name=None
    pos_mob_2_name=None
    pos_mob_3_name=None

    generic_names = [
    "neg_sampleflow_name",
    "pos_sampleflow_name",
    "neg_temperature_name",
    "pos_temperature_name",
    "neg_pressure_name",
    "pos_pressure_name",
    "neg_cluster_mob_name",
    "pos_cluster_mob_name",
    "neg_conc_1_name",
    "pos_conc_1_name",
    "neg_conc_2_name",
    "pos_conc_2_name",
    "neg_conc_3_name",
    "pos_conc_3_name",
    "neg_mob_1_name",
    "pos_mob_1_name",
    "neg_mob_2_name",
    "pos_mob_2_name",
    "neg_mob_3_name",
    "pos_mob_3_name"
    ]

    name_values = [
    neg_sampleflow_name,
    pos_sampleflow_name,
    neg_temperature_name,
    pos_temperature_name,
    neg_pressure_name,
    pos_pressure_name,
    neg_cluster_mob_name,
    pos_cluster_mob_name,
    neg_conc_1_name,
    pos_conc_1_name,
    neg_conc_2_name,
    pos_conc_2_name,
    neg_conc_3_name,
    pos_conc_3_name,
    neg_mob_1_name,
    pos_mob_1_name,
    neg_mob_2_name,
    pos_mob_2_name,
    neg_mob_3_name,
    pos_mob_3_name
    ]

    name_dict = {}

    for i in range(len(NAME_LIST)):
        for test_name in NAME_LIST[i]:
            if test_name in params:
                name_values[i] = test_name
                break

    for name,value in zip(generic_names,name_values):
        name_dict[name] = value

    return name_dict

def get_data(records):

    if records is None:
        return None

    else:        
        # Check that the relevant diagnostic data is found
        data_names = find_names(list(records))

        # TEMPERATURE (K)
        if (data_names["neg_temperature_name"] is not None):
             neg_temperature = 273.15 + records[data_names["neg_temperature_name"]].astype(float)    
        else:
             neg_temperature = None

        # TEMPERATURE (K)
        if (data_names["pos_temperature_name"] is not None):
             pos_temperature = 273.15 + records[data_names["pos_temperature_name"]].astype(float)
        else:
             pos_temperature = None

        # PRESSURE (Pa)
        if (data_names["neg_pressure_name"] is not None):
             neg_pressure = 100.0 * records[data_names["neg_pressure_name"]].astype(float)
        else:
             neg_pressure = None

        if (data_names["pos_pressure_name"] is not None):
             pos_pressure = 100.0 * records[data_names["pos_pressure_name"]].astype(float)
        else:
             pos_pressure = None

        # SAMPLEFLOW (lpm)
        if (data_names["neg_sampleflow_name"] is not None):
             neg_sampleflow = records[data_names["neg_sampleflow_name"]].astype(float)
        else:
             neg_sampleflow = None
        
        if (data_names["pos_sampleflow_name"] is not None):
             pos_sampleflow = records[data_names["pos_sampleflow_name"]].astype(float)
        else:
             pos_sampleflow = None

        # MOBILITIES
        if ((data_names["neg_cluster_mob_name"] is not None) and
            (data_names["neg_mob_1_name"] is not None) and
            (data_names["neg_mob_2_name"] is not None) and
            (data_names["neg_mob_3_name"] is not None)
            ):
            neg_cluster_mob = records[data_names["neg_cluster_mob_name"]].astype(float) # cm2/sV
            neg_mob_1 = records[data_names["neg_mob_1_name"]].astype(float) # cm2/sV
            neg_mob_2 = records[data_names["neg_mob_2_name"]].astype(float) # cm2/sV
            neg_mob_3 = records[data_names["neg_mob_3_name"]].astype(float) # cm2/sV
        else:
            neg_cluster_mob = None
            neg_mob_1 = None
            neg_mob_2 = None
            neg_mob_3 = None

        if ((data_names["pos_cluster_mob_name"] is not None) and
            (data_names["pos_mob_1_name"] is not None) and
            (data_names["pos_mob_2_name"] is not None) and
            (data_names["pos_mob_3_name"] is not None)
            ):
            pos_cluster_mob = records[data_names["pos_cluster_mob_name"]].astype(float) # cm2/sV
            pos_mob_1 = records[data_names["pos_mob_1_name"]].astype(float) # cm2/sV
            pos_mob_2 = records[data_names["pos_mob_2_name"]].astype(float) # cm2/sV
            pos_mob_3 = records[data_names["pos_mob_3_name"]].astype(float) # cm2/sV
        else:
            pos_cluster_mob = None
            pos_mob_1 = None
            pos_mob_2 = None
            pos_mob_3 = None

        # CONCENTRATIONS
        if ((data_names["neg_conc_1_name"] is not None) and 
            (data_names["neg_conc_2_name"] is not None) and
            (data_names["neg_conc_3_name"] is not None)):
            neg_conc_1 = records[data_names["neg_conc_1_name"]].astype(float)
            neg_conc_2 = records[data_names["neg_conc_2_name"]].astype(float)
            neg_conc_3 = records[data_names["neg_conc_3_name"]].astype(float)
            neg_cluster_conc = (neg_conc_1 + neg_conc_2) - neg_conc_3
        else:
            neg_conc_1 = None
            neg_conc_2 = None
            neg_conc_3 = None
            neg_cluster_conc = None

        if ((data_names["pos_conc_1_name"] is not None) and 
            (data_names["pos_conc_2_name"] is not None) and
            (data_names["pos_conc_3_name"] is not None)):
            pos_conc_1 = records[data_names["pos_conc_1_name"]].astype(float)
            pos_conc_2 = records[data_names["pos_conc_2_name"]].astype(float)
            pos_conc_3 = records[data_names["pos_conc_3_name"]].astype(float)
            pos_cluster_conc = (pos_conc_1 + pos_conc_2) - pos_conc_3
        else:
            pos_conc_1 = None 
            pos_conc_2 = None 
            pos_conc_3 = None 
            pos_cluster_conc = None

        # Sanity check the values
        if neg_temperature is not None:
            neg_temperature = neg_temperature.where(((neg_temperature>=223.)&(neg_temperature<=353.)),np.nan)

        if pos_temperature is not None:
            pos_temperature = pos_temperature.where(((pos_temperature>=223.)&(pos_temperature<=353.)),np.nan)

        if neg_pressure is not None:
            neg_pressure = neg_pressure.where(((neg_pressure>=37000.)&(neg_pressure<=121000.)),np.nan) 

        if pos_pressure is not None:
            pos_pressure = pos_pressure.where(((pos_pressure>=37000.)&(pos_pressure<=121000.)),np.nan) 
        
        if neg_sampleflow is not None:
            neg_sampleflow = neg_sampleflow.where(((neg_sampleflow>=5.)&(neg_sampleflow<=65.)),np.nan)

        if pos_sampleflow is not None:
            pos_sampleflow = pos_sampleflow.where(((pos_sampleflow>=5.)&(pos_sampleflow<=65.)),np.nan)


        # Make a dictionary out of the values
        data = {
            "neg_temperature": neg_temperature,
            "neg_pressure": neg_pressure,
            "pos_temperature": pos_temperature,
            "pos_pressure": pos_pressure,
            "neg_sampleflow": neg_sampleflow,
            "pos_sampleflow": neg_sampleflow,
            "neg_cluster_mob": neg_cluster_mob,
            "pos_cluster_mob": pos_cluster_mob,
            "neg_cluster_conc": neg_cluster_conc,
            "pos_cluster_conc": pos_cluster_conc,
            "neg_conc_1": neg_conc_1,
            "neg_conc_2": neg_conc_2,
            "neg_conc_3": neg_conc_3,
            "pos_conc_1": pos_conc_1,
            "pos_conc_2": pos_conc_2,
            "pos_conc_3": pos_conc_3,
            "neg_mob_1": neg_mob_1,
            "neg_mob_2": neg_mob_2,
            "neg_mob_3": neg_mob_3,
            "pos_mob_1": pos_mob_1,
            "pos_mob_2": pos_mob_2,
            "pos_mob_3": pos_mob_3
            }

        return data


def bring_to_sealevel(data):

    if (data is None):
        return None

    # Negative polarity
    if ((data["neg_cluster_conc"] is None) or 
        (data["neg_conc_1"] is None) or
        (data["neg_conc_2"] is None) or
        (data["neg_conc_3"] is None) or 
        (data["neg_temperature"] is None) or 
        (data["neg_pressure"] is None)):
        data["neg_cluster_conc"] = None
        data["neg_conc_1"]=None 
        data["neg_conc_2"]=None
        data["neg_conc_3"]=None
    else:
        stp_corr_factor = (PRES_REF*data["neg_temperature"].values)/(TEMP_REF*data["neg_pressure"].values)
        data["neg_cluster_conc"] = stp_corr_factor * data["neg_cluster_conc"]
        data["neg_conc_1"] = stp_corr_factor * data["neg_conc_1"]
        data["neg_conc_2"] = stp_corr_factor * data["neg_conc_2"]
        data["neg_conc_3"] = stp_corr_factor * data["neg_conc_3"]

    # Positive polarity
    if ((data["pos_cluster_conc"] is None) or 
        (data["pos_conc_1"] is None) or
        (data["pos_conc_2"] is None) or
        (data["pos_conc_3"] is None) or 
        (data["pos_temperature"] is None) or 
        (data["pos_pressure"] is None)):
        data["pos_cluster_conc"] = None
        data["pos_conc_1"]=None 
        data["pos_conc_2"]=None
        data["pos_conc_3"]=None
    else:
        stp_corr_factor = (PRES_REF*data["pos_temperature"].values)/(TEMP_REF*data["pos_pressure"].values)
        data["pos_cluster_conc"] = stp_corr_factor * data["pos_cluster_conc"]
        data["pos_conc_1"] = stp_corr_factor * data["pos_conc_1"]
        data["pos_conc_2"] = stp_corr_factor * data["pos_conc_2"]
        data["pos_conc_3"] = stp_corr_factor * data["pos_conc_3"]

    return data


def correct_inlet_losses(data,pipe_length):

    if (data is None):
        return None

    elif ((data["neg_cluster_conc"] is None) or
        (data["neg_conc_1"] is None) or
        (data["neg_conc_2"] is None) or
        (data["neg_conc_3"] is None) or
        (data["neg_temperature"] is None) or 
        (data["neg_pressure"] is None) or 
        (data["neg_sampleflow"] is None) or
        (data["neg_cluster_mob"] is None) or
        (data["neg_mob_1"] is None) or
        (data["neg_mob_2"] is None) or
        (data["neg_mob_3"] is None)
        ):
        data["neg_cluster_conc"] = None
        data["neg_conc_1"] = None
        data["neg_conc_2"] = None
        data["neg_conc_3"] = None
    else:
        avg_cluster_dp = af.mob2diam(data["neg_cluster_mob"].mean(),data["neg_temperature"].mean(),data["neg_pressure"].mean())
        avg_dp_1 = af.mob2diam(data["neg_mob_1"].mean(),data["neg_temperature"].mean(),data["neg_pressure"].mean())
        avg_dp_2 = af.mob2diam(data["neg_mob_2"].mean(),data["neg_temperature"].mean(),data["neg_pressure"].mean())
        avg_dp_3 = af.mob2diam(data["neg_mob_3"].mean(),data["neg_temperature"].mean(),data["neg_pressure"].mean())

        dp = pd.Series([avg_cluster_dp,avg_dp_1,avg_dp_2,avg_dp_3])
        
        throughput = af.tubeloss(
            dp,
            data["neg_sampleflow"],
            pipe_length,
            data["neg_temperature"],
            data["neg_pressure"])
        
        data["neg_cluster_conc"] = data["neg_cluster_conc"] / throughput.values[:,0]
        data["neg_conc_1"] = data["neg_conc_1"] / throughput.values[:,1]
        data["neg_conc_2"] = data["neg_conc_2"] / throughput.values[:,2]
        data["neg_conc_3"] = data["neg_conc_3"] / throughput.values[:,3]


    if ((data["pos_cluster_conc"] is None) or
        (data["pos_conc_1"] is None) or
        (data["pos_conc_2"] is None) or
        (data["pos_conc_3"] is None) or
        (data["pos_temperature"] is None) or 
        (data["pos_pressure"] is None) or 
        (data["pos_sampleflow"] is None) or
        (data["pos_cluster_mob"] is None) or
        (data["pos_mob_1"] is None) or
        (data["pos_mob_2"] is None) or
        (data["pos_mob_3"] is None)
        ):
        data["pos_cluster_conc"] = None
        data["pos_conc_1"] = None
        data["pos_conc_2"] = None
        data["pos_conc_3"] = None
    else:
        avg_cluster_dp = af.mob2diam(data["pos_cluster_mob"].mean(),data["pos_temperature"].mean(),data["pos_pressure"].mean())
        avg_dp_1 = af.mob2diam(data["pos_mob_1"].mean(),data["pos_temperature"].mean(),data["pos_pressure"].mean())
        avg_dp_2 = af.mob2diam(data["pos_mob_2"].mean(),data["pos_temperature"].mean(),data["pos_pressure"].mean())
        avg_dp_3 = af.mob2diam(data["pos_mob_3"].mean(),data["pos_temperature"].mean(),data["pos_pressure"].mean())

        dp = pd.Series([avg_cluster_dp,avg_dp_1,avg_dp_2,avg_dp_3])
        
        throughput = af.tubeloss(
            dp,
            data["pos_sampleflow"],
            pipe_length,
            data["pos_temperature"],
            data["pos_pressure"])
        
        data["pos_cluster_conc"] = data["pos_cluster_conc"] / throughput.values[:,0]
        data["pos_conc_1"] = data["pos_conc_1"] / throughput.values[:,1]
        data["pos_conc_2"] = data["pos_conc_2"] / throughput.values[:,2]
        data["pos_conc_3"] = data["pos_conc_3"] / throughput.values[:,3]

    return data

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
                                    
            # Message/flag only concerns negative polarity
            elif ("−" in message):   
                flags_neg_spectra.loc[combined_idx,message] = 1
            
            # Message/flag concerns both polarities
            else:
                flags_pos_spectra.loc[combined_idx,message] = 1
                flags_neg_spectra.loc[combined_idx,message] = 1
            
        return flags_neg_spectra, flags_pos_spectra


def save_as_netcdf(
    netcdf_save_path,
    data,
    negion_flags,
    posion_flags,
    flag_explanations,
    measurement_info):


    if (data is None):
        return False    
    elif ((data["neg_cluster_conc"] is None) and
        (data["neg_conc_1"] is None) and
        (data["neg_conc_2"] is None) and
        (data["neg_conc_3"] is None) and
        (data["pos_cluster_conc"] is None) and
        (data["pos_conc_1"] is None) and
        (data["pos_conc_2"] is None) and
        (data["pos_conc_3"] is None)
        ):
        return False
    else:    
        for spectra in [
            data["neg_cluster_conc"],
            data["neg_conc_1"],
            data["neg_conc_2"],
            data["neg_conc_3"], 
            data["pos_cluster_conc"],
            data["pos_conc_1"],
            data["pos_conc_2"],
            data["pos_conc_3"]
            ]:
            if spectra is not None:
                time = spectra.index.values
                LEN_TIME = len(time)
                nan_data = np.nan*np.ones(LEN_TIME)
                break

        ds = xr.Dataset()
        ds = ds.assign_coords(
            coords = {
                "time": time,
            }
        )
        
        ds.time.attrs["timezone"] = "utc"
                
        if data["neg_cluster_conc"] is not None:
            ds = ds.assign(neg_cluster_conc=(("time",),data["neg_cluster_conc"].values))
        else:
            ds = ds.assign(neg_cluster_conc=(("time",),nan_data))

        if data["neg_conc_1"] is not None:
            ds = ds.assign(neg_conc_1=(("time",),data["neg_conc_1"].values))
        else:
            ds = ds.assign(neg_conc_1=(("time",),nan_data))

        if data["neg_conc_2"] is not None:
            ds = ds.assign(neg_conc_2=(("time",),data["neg_conc_2"].values))
        else:
            ds = ds.assign(neg_conc_2=(("time",),nan_data))

        if data["neg_conc_3"] is not None:
            ds = ds.assign(neg_conc_3=(("time",),data["neg_conc_3"].values))
        else:
            ds = ds.assign(neg_conc_3=(("time",),nan_data))

        if data["pos_cluster_conc"] is not None:
            ds = ds.assign(pos_cluster_conc=(("time",),data["pos_cluster_conc"].values))
        else:
            ds = ds.assign(pos_cluster_conc=(("time",),nan_data))

        if data["pos_conc_1"] is not None:
            ds = ds.assign(pos_conc_1=(("time",),data["pos_conc_1"].values))
        else:
            ds = ds.assign(pos_conc_1=(("time",),nan_data))

        if data["pos_conc_2"] is not None:
            ds = ds.assign(pos_conc_2=(("time",),data["pos_conc_2"].values))
        else:
            ds = ds.assign(pos_conc_2=(("time",),nan_data))

        if data["pos_conc_3"] is not None:
            ds = ds.assign(pos_conc_3=(("time",),data["pos_conc_3"].values))
        else:
            ds = ds.assign(pos_conc_3=(("time",),nan_data))


        if data["neg_cluster_mob"] is not None:
            ds = ds.assign(neg_cluster_mob=(("time",),data["neg_cluster_mob"].values))
        else:
            ds = ds.assign(neg_cluster_mob=(("time",),nan_data))

        if data["neg_mob_1"] is not None:
            ds = ds.assign(neg_mob_1=(("time",),data["neg_mob_1"].values))
        else:
            ds = ds.assign(neg_mob_1=(("time",),nan_data))

        if data["neg_mob_2"] is not None:
            ds = ds.assign(neg_mob_2=(("time",),data["neg_mob_2"].values))
        else:
            ds = ds.assign(neg_mob_2=(("time",),nan_data))

        if data["neg_mob_3"] is not None:
            ds = ds.assign(neg_mob_3=(("time",),data["neg_mob_3"].values))
        else:
            ds = ds.assign(neg_mob_3=(("time",),nan_data))

        if data["pos_cluster_mob"] is not None:
            ds = ds.assign(pos_cluster_mob=(("time",),data["pos_cluster_mob"].values))
        else:
            ds = ds.assign(pos_cluster_mob=(("time",),nan_data))

        if data["pos_mob_1"] is not None:
            ds = ds.assign(pos_mob_1=(("time",),data["pos_mob_1"].values))
        else:
            ds = ds.assign(pos_mob_1=(("time",),nan_data))

        if data["pos_mob_2"] is not None:
            ds = ds.assign(pos_mob_2=(("time",),data["pos_mob_2"].values))
        else:
            ds = ds.assign(pos_mob_2=(("time",),nan_data))

        if data["pos_mob_3"] is not None:
            ds = ds.assign(pos_mob_3=(("time",),data["pos_mob_3"].values))
        else:
            ds = ds.assign(pos_mob_3=(("time",),nan_data))



        if data["neg_temperature"] is not None:
            ds = ds.assign(neg_temperature=(("time",),data["neg_temperature"].values))
        else:
            ds = ds.assign(neg_temperature=(("time",),nan_data))
            
        if data["pos_temperature"] is not None:
            ds = ds.assign(pos_temperature=(("time",),data["pos_temperature"].values))
        else:
            ds = ds.assign(pos_temperature=(("time",),nan_data))

        if data["neg_pressure"] is not None:
            ds = ds.assign(neg_pressure=(("time",),data["neg_pressure"].values))
        else:
            ds = ds.assign(neg_pressure=(("time",),nan_data))
            
        if data["pos_pressure"] is not None:
            ds = ds.assign(pos_pressure=(("time",),data["pos_pressure"].values))
        else:
            ds = ds.assign(pos_pressure=(("time",),nan_data))

        if data["neg_sampleflow"] is not None:
            ds = ds.assign(neg_sampleflow=(("time",),data["neg_sampleflow"].values))
        else:
            ds = ds.assign(neg_sampleflow=(("time",),nan_data))
            
        if data["pos_sampleflow"] is not None:
            ds = ds.assign(pos_sampleflow=(("time",),data["pos_sampleflow"].values))
        else:
            ds = ds.assign(pos_sampleflow=(("time",),nan_data))


        ds.neg_cluster_conc.attrs["units"] = "cm-3"
        ds.neg_cluster_conc.attrs["description"] = "Negative cluster ion number concentration"
        ds.neg_conc_1.attrs["units"] = "cm-3"
        ds.neg_conc_1.attrs["description"] = "Negative ion number concentration in channel 1"
        ds.neg_conc_2.attrs["units"] = "cm-3"
        ds.neg_conc_2.attrs["description"] = "Negative ion number concentration in channel 2"
        ds.neg_conc_3.attrs["units"] = "cm-3"
        ds.neg_conc_3.attrs["description"] = "Negative ion number concentration in channel 3"

        ds.pos_cluster_conc.attrs["units"] = "cm-3"
        ds.pos_cluster_conc.attrs["description"] = "Positive cluster ion number concentration"
        ds.pos_conc_1.attrs["units"] = "cm-3"
        ds.pos_conc_1.attrs["description"] = "Positive ion number concentration in channel 1"
        ds.pos_conc_2.attrs["units"] = "cm-3"
        ds.pos_conc_2.attrs["description"] = "Positive ion number concentration in channel 2"
        ds.pos_conc_3.attrs["units"] = "cm-3"
        ds.pos_conc_3.attrs["description"] = "Positive ion number concentration in channel 3"

        ds.neg_cluster_mob.attrs["units"] = "cm2s-1V-1"
        ds.neg_cluster_mob.attrs["description"] = "Negative cluster ion mobility"
        ds.neg_mob_1.attrs["units"] = "cm2s-1V-1"
        ds.neg_mob_1.attrs["description"] = "Negative ion mobility in channel 1"
        ds.neg_mob_2.attrs["units"] = "cm2s-1V-1"
        ds.neg_mob_2.attrs["description"] = "Negative ion mobility in channel 2"
        ds.neg_mob_3.attrs["units"] = "cm2s-1V-1"
        ds.neg_mob_3.attrs["description"] = "Negative ion mobility in channel 3"

        ds.pos_cluster_mob.attrs["units"] = "cm2s-1V-1"
        ds.pos_cluster_mob.attrs["description"] = "Positive cluster ion mobility"
        ds.pos_mob_1.attrs["units"] = "cm2s-1V-1"
        ds.pos_mob_1.attrs["description"] = "Positive ion mobility in channel 1"
        ds.pos_mob_2.attrs["units"] = "cm2s-1V-1"
        ds.pos_mob_2.attrs["description"] = "Positive ion mobility in channel 2"
        ds.pos_mob_3.attrs["units"] = "cm2s-1V-1"
        ds.pos_mob_3.attrs["description"] = "Positive ion mobility in channel 3"

        ds.neg_temperature.attrs["units"] = "K"
        ds.neg_temperature.attrs["description"] = "Negative polarity sample air temperature"

        ds.pos_temperature.attrs["units"] = "K"
        ds.pos_temperature.attrs["description"] = "Positive polarity sample air temperature"

        ds.neg_pressure.attrs["units"] = "Pa"
        ds.neg_pressure.attrs["description"] = "Negative polarity sample air pressure"

        ds.pos_pressure.attrs["units"] = "Pa"
        ds.pos_pressure.attrs["description"] = "Positive polarity sample air pressure"

        ds.neg_sampleflow.attrs["units"] = "lpm"
        ds.neg_sampleflow.attrs["description"] = "Flow rate at the negative polarity outlet"

        ds.pos_sampleflow.attrs["units"] = "lpm"
        ds.pos_sampleflow.attrs["description"] = "Flow rate at the positive polarity outlet"

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
            
        # Add measurement info
        ds = ds.assign_attrs(measurement_info)
        
        ds.to_netcdf(netcdf_save_path)
        
        return True
       
def cic_processor(config_file):
    """ Process CIC data
    
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
        instrument_model = config['instrument_model']
        longitude = config["longitude"]
        latitude = config["latitude"]
        end_date = config['end_date']
        allow_reprocess = config["allow_reprocess"]
        redo_database = config["redo_database"]
        pipelength = config['inlet_length']
        do_inlet_loss_correction = config['do_inlet_loss_correction']
        convert_to_standard_conditions = config['convert_to_standard_conditions']
        file_format = config["file_format"]
        resolution = config["resolution"]
        
    # Check the config file
    assert isinstance(start_date,date)
    assert (end_date=='' or isinstance(end_date,date))
    assert os.path.exists(save_path)
    assert all([os.path.exists(x) for x in load_path])
    assert isinstance(allow_reprocess,bool)
    assert isinstance(convert_to_standard_conditions,bool)
    assert isinstance(do_inlet_loss_correction,bool)
    assert isinstance(pipelength,float)
    assert isinstance(longitude,float)
    assert isinstance(latitude,float)
    assert isinstance(redo_database,bool)
    pd.tseries.frequencies.to_offset(resolution)
    
    # Extract relevant info for metadata from the config
    measurement_info = {
        'measurement_location':location,
        'description':description,
        'instrument_model':instrument_model,
        'longitude':longitude,
        'latitude':latitude,
        'inlet_length':pipelength,
        'do_inlet_loss_correction':str(do_inlet_loss_correction),
        'convert_to_standard_conditions':str(convert_to_standard_conditions),
        "resolution":resolution
    }    

    end_date = date.today() if end_date=='' else end_date

    if redo_database:
        try:
            os.remove(database)
        except:
            pass

    db = TinyDB(database)
    check = Query()

    filename_formats = [[f"%Y%m%d-{file_format}.records"]]

    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)

    start_date_str = start_dt.strftime("%Y%m%d")
    end_date_str = end_dt.strftime("%Y%m%d")

    # List existing dates based on if diagnostic file was found
    list_of_existing_dates = [x["timestamp"] for x in db.search(check.records.exists())]

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

                    records_filename = os.path.join(load_path[i],x.strftime(filename_formats[j][0]))

                    if os.path.exists(records_filename):

                        date_str = x.strftime("%Y%m%d")

                        db.insert(
                            {"timestamp":date_str,
                            "records":records_filename}
                            )

                        files_found=True

                        if (load_path[i]!=preferred_path):
                            preferred_path = load_path.pop(i)
                            load_path.insert(0,preferred_path)

                        if (filename_formats[i]!=preferred_format):
                            preferred_format = filename_formats.pop(j)
                            filename_formats.insert(0,preferred_format)

                        print(f"Added {x.strftime('%Y%m%d')} to database for {location}")
                        
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
            ((check.records.exists()) &
            (~check.processed_file.exists()) &
            (check.timestamp>=start_date_str) &
            (check.timestamp<=end_date_str))))
    else:
        database_iterator = iter(db.search(
            (check.recrods.exists() &
            (check.timestamp>=last_day) &
            (check.timestamp>=start_date_str) &
            (check.timestamp<=end_date_str))))

    for x in database_iterator:

        print("Processing %s (%s)" % (x["timestamp"], location))

        (ion_records, 
            ion_flags, 
            offset_flags, 
            flag_explanations) = read_raw(x["records"],x["timestamp"],resolution)

        negion_flags, posion_flags = flags2polarity(ion_flags, offset_flags, flag_explanations)
           
        data = get_data(ion_records)

        if convert_to_standard_conditions:
            data = bring_to_sealevel(data)
         
        if do_inlet_loss_correction:
            data = correct_inlet_losses(data,pipelength)

        my_save_path = os.path.join(save_path,"CIC_"+x["timestamp"]+".nc")
        
        saved = save_as_netcdf(
            my_save_path,
            data,
            negion_flags,
            posion_flags,
            flag_explanations,
            measurement_info
        )
        
        if saved:
            db.update({"processed_file": my_save_path},check.timestamp==x["timestamp"])
            
    print("Done!")
