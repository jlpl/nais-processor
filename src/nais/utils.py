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

def remove_flagged_rows(ds,flag):
    """
    
    Parameters
    ----------
    
    ds : xarray.Dataset
        NAIS dataset
    flag : str
	    Flag to be removed
    
    Returns
    -------
    
    xarray.Dataset
        NAIS dataset with flag rows set to NaN
    
    """
    
    ds2 = ds.copy(deep=True)
        
    idx0 = ds2.neg_ion_flags.sel(flag=flag)==1
    idx1 = ds2.pos_ion_flags.sel(flag=flag)==1
    idx2 = ds2.neg_particle_flags.sel(flag=flag)==1
    idx3 = ds2.pos_particle_flags.sel(flag=flag)==1
    
    ds2.neg_ions[idx0,:] = np.nan
    ds2.pos_ions[idx1,:] = np.nan
    ds2.neg_particles[idx2,:] = np.nan
    ds2.pos_particles[idx3,:] = np.nan     

    return ds2
    
def combine_databases(database_list, combined_database):
    """Combine JSON databases

    Parameters
    ----------

    database_list : list of str
        List of full paths to databases that should be combined

        First database should have the earliest data, second database
        the second earliest and so on
    combined_database : str
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
               
def combine_data(
    source_dir, 
    date_range, 
    time_reso, 
    flag_sensitivity=0.5):
    """
    
    Parameters
    ----------
    
    source_dir : str
        Directory for NAIS datafiles    
    date_range : pandas.DatetimeIndex
        Range of dates for combining data        
    time_reso : str
        A pandas date frequency string
        
        See for all options here: 
        https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases
    flag_sensitivity : float
        fraction of time flag needs to be present 
        in resampling
    
    Returns
    -------
    
    xarray.Dataset or None
        Combined dataset, None if no data in the 
        date range
    
    """
    
    assert pd.Timedelta(time_reso)>pd.Timedelta("5min")
    
    data_read = False
    for date in date_range:
        
        filename_date = os.path.join(source_dir,"NAIS_"+date.strftime("%Y%m%d")+".nc")
        
        if os.path.isfile(filename_date):       
            ds = xr.open_dataset(filename_date)
            ds_data = ds[["neg_ions","pos_ions","neg_particles","pos_particles"]]
            ds_flags = ds[["neg_ion_flags","pos_ion_flags","neg_particle_flags","pos_particle_flags"]]
            ds_data_resampled = ds_data.resample({"time":time_reso}).median()
            ds_flags_resampled = xr.where(
                ds_flags.resample({"time":time_reso}).mean()>flag_sensitivity,1,0
            )
            
            #ds_flags_resampled[]
            ds.close()

            if data_read==False:
                ds_data_combined = ds_data_resampled
                ds_flags_combined = ds_flags_resampled
                data_read = True
            else:
                ds_data_combined = xr.concat(
                    (ds_data_combined,ds_data_resampled),
                    dim="time"
                )
                ds_flags_combined = xr.concat(
                    (ds_flags_combined,ds_flags_resampled),
                    dim="time",
                    fill_value=0
                )
                       
    if data_read:
        ds_data_combined_resampled = ds_data_combined.resample({"time":time_reso}).median()
        ds_flags_combined_resampled = xr.where(
            ds_flags_combined.resample({"time":time_reso}).mean()>flag_sensitivity,1,0
        )
        
        # Combine the two data frames
        ds_data_and_flags =  xr.merge((ds_data_combined_resampled,ds_flags_combined_resampled))
        
        idx_final=pd.date_range(date_range[0],date_range[-1],freq=time_reso)
        ds_final=ds_data_and_flags.reindex(
            indexers={"time":idx_final.values}, 
            method="nearest",
            tolerance=time_reso
        )
        
        return ds_final
    
    else:
        return None

		

