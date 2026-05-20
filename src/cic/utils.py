import numpy as np
import pandas as pd
import os
import aerosol.functions as af
import xarray as xr

def remove_flagged_rows(ds,flag):
    """
    Set data associated with given flag as NaN
    
    Parameters
    ----------
    
    ds : xarray.Dataset
        CIC dataset
    flag : str
        Flag to be removed
    
    Returns
    -------
    
    xarray.Dataset
        CIC dataset with flag rows set to NaN
    
    """
    
    ds2 = ds.copy(deep=True)

    idx0 = ds2.neg_ion_flags.sel(flag=flag)==1
    idx1 = ds2.pos_ion_flags.sel(flag=flag)==1
    
    ds2.neg_ions[idx0,:] = np.nan
    ds2.pos_ions[idx1,:] = np.nan

    return ds2
               
def combine_data(
    files, 
    start,
    end, 
    time_reso, 
    flag_sensitivity=0.5):
    """
    Combine netcdf datafiles and resample to new resolution
    on continuous time index  
    
    Parameters
    ----------
    
    files : list
        List of CIC data file paths    
    start : str
        start time
    end : str
        end time   
    time_reso : str
        A pandas date frequency string. See for all options here: 
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
    
    data_read = False
    for f in files:
        
        filename_date = f
        
        if os.path.isfile(filename_date):       
            ds = xr.open_dataset(filename_date)
            
            data_variables = []
            flag_variables = []
            for v in list(ds):
                if ("_flags" in v):
                    flag_variables.append(v)
                else:
                    data_variables.append(v)
            ds_data = ds[data_variables]
            ds_flags = ds[flag_variables]

            # Sort and uniqueify
            ds_data = ds_data.sortby("time")
            unique_mask = ~pd.Index(ds_data.time.values).duplicated()
            ds_data = ds_data.isel(time=unique_mask)
    
            ds_flags = ds_flags.sortby("time")
            unique_mask = ~pd.Index(ds_flags.time.values).duplicated()
            ds_flags = ds_flags.isel(time=unique_mask)
            
            ds_data_resampled = ds_data.resample({"time":time_reso}).median()
            ds_flags_resampled = xr.where(
                ds_flags.resample({"time":time_reso}).mean()>flag_sensitivity,1,0
            )
            
            ds.close()

            if data_read==False:
                ds_data_combined = ds_data_resampled
                ds_flags_combined = ds_flags_resampled
                place = ds.attrs["measurement_location"]
                place_id = ds.attrs["id"]
                lon = ds.attrs["longitude"]
                lat = ds.attrs["latitude"]
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
        # Sort and uniqueify
        ds_data_combined = ds_data_combined.sortby("time")
        unique_mask = ~pd.Index(ds_data_combined.time.values).duplicated()
        ds_data_combined = ds_data_combined.isel(time=unique_mask)

        ds_flags_combined = ds_flags_combined.sortby("time")
        unique_mask = ~pd.Index(ds_flags_combined.time.values).duplicated()
        ds_flags_combined = ds_flags_combined.isel(time=unique_mask)

        ds_data_combined_resampled = ds_data_combined.resample({"time":time_reso}).median()
        ds_flags_combined_resampled = xr.where(
            ds_flags_combined.resample({"time":time_reso}).mean()>flag_sensitivity,1,0
        )
        
        # Combine the two data frames
        ds_data_and_flags =  xr.merge((ds_data_combined_resampled,ds_flags_combined_resampled))
        
        idx_final = pd.date_range(start,end,freq=time_reso)
        ds_final = ds_data_and_flags.reindex(
            indexers={"time":idx_final.values}, 
            method="nearest",
            tolerance=time_reso
        )

        ds_final.attrs = {
            "measurement_location": place,
            "id": place_id,
            "longitude": lon,
            "latitude": lat
        }
        
        return ds_final
    
    else:
        return None
