![data](./small_img.png)

# NAIS processor
Use this code package to process [NAIS](https://www.airel.ee/products/nais/) (Neutral cluster and Air Ion Spectrometer, Airel Ltd.) data files.

The code corrects for diffusion losses in the inlet line (Gromley and Kennedy, 1948) and applies an ion mode calibration (Wagner et al. 2016). Optionally the data can be corrected to standard conditions (273.15 K, 101325 Pa), which can be useful when comparing aerosol particle and ion data from various locations at different altitudes.

Optionally one can also apply a cleaning procedure to the data where the corona ion band is removed from the particle data and instances of electrometer noise are removed from ion and particle data.

[Documentation](https://jlpl.github.io/nais-processor/)

## Installation
```shell
pip install nais-processor
```
## Example usage
Open the python prompt and load methods from the `nais_processor` module.
Then use the [`make_config()`](https://jlpl.github.io/nais-processor/#nais_processor.make_config) method to create a configuration file that is used at processing the data files.
```
$ python
>>> from nais_processor import *
>>> make_config()

Enter absolute path to configuration file.
For example: /home/user/campaign.yml
> /home/user/viikki.yml

Enter absolute path(s) to raw data folder(s). Separate multiple paths with comma.
For example: /home/user/data/2021,/home/user/data/2022
> /home/user/data/2021,/home/user/data/2022

Enter absolute path to processed data folder.
For example: /home/user/campaign
> /home/user/viikki

Enter start of measurement (YYYY-MM-DD)
Leave empty if you want to start from the earliest date
> 2022-09-28   

Enter end of measurement (YYYY-MM-DD)
If empty processor assumes current day
> 2022-09-30

Enter absolute path to database file
For example: /home/user/campaign.json
> /home/user/viikki.json 

Allow reprocessing (True/False)
Overwrites already processed data when re-running the processor
> True

Enter measurement location
For example: Helsinki, Kumpula
> Viikki, Helsinki, Finland 

Apply data cleaning procedures (True/False)
Attempt to remove corona ions and electrometer noise from data
> True

Apply corrections (True/False)
Requires a NAIS with temperature and pressure sensors.
> True 

Length of the inlet in meters
> 1.0 

Correct concentrations to sealevel conditions (True/False)
> True

Configuration saved to: /home/user/viikki.yml
```
The resulting configuration file looks like this:
```yaml
allow_reprocess: true
apply_cleaning: true
apply_corrections: true
data_folder:
- /home/user/data/2021 
- /home/user/data/2022
database_file: /home/user/viikki.json
end_date: 2022-09-30
inlet_length: 1.0
location: Viikki, Helsinki, Finland
processed_folder: /home/user/viikki
sealevel_correction: true
start_date: 2022-09-28
```
Then process the data files by running [`nais_processor()`](https://jlpl.github.io/nais-processor/#nais_processor.nais_processor) method with the config file as the input argument.
```
>>> nais_processor("/home/user/viikki.yml")
/home/user/viikki.yml
building database...
processing 20220928 (Viikki, Helsinki, Finland)
processing 20220929 (Viikki, Helsinki, Finland)
processing 20220930 (Viikki, Helsinki, Finland)
Done!
```
The code produces daily processed data files for ion and particle data. These files are saved in the destinations given in the configuration file.

The processed data files are named

`NAIS[n|p][yyyymmdd][np|nds].sum`

where `n` and `p` refer to negative and positive polarity respectively. `yyyymmdd` tells the date in the year-month-day format. `np` and `nds` refer to particle and ion data respectively. The cleaned files have an additional `_cleaned` at the end of the filename before the suffix.

The data files have the following structure (sum matrix)
```
[0,0]  = UTC offset in hours
[1:,0] = Time (MATLAB datenum) 
[0,2:] = Geometric mean diameter of size-channel (m)
[1:,1] = Integrated total number concentration (cm-3)
[1:,2:] = Normalized number concentrations, dN/dlogDp (cm-3)
```
The locations of raw files, processed files and cleaned processed files are written in the `database_file`, which is in JSON format.

### Combining sumfiles
Once you have processed your NAIS data you can extract any time range in a sum matrix format using the [`combine_spectra()`](https://jlpl.github.io/nais-processor/#nais_processor.combine_spectra) function.

Example:
```python
import nais_processor as nais

start_time="2022-09-29 02:00:00"
end_time="2022-09-30 14:00:00"

combined_data = nais.combine_spectra(
    config_file,start_time,end_time,spectra_type="negion")
```
You can also plot sum matrices using [`plot_sumfile()`](https://jlpl.github.io/nais-processor/#nais_processor.plot_sumfile) or you can create daily plots using [`plot_nais()`](https://jlpl.github.io/nais-processor/#nais_processor.plot_nais) without playing around with sum matrices.

## License
This project is licensed under the terms of the GNU GPLv3.

## References
Gormley P. G. and Kennedy M., Diffusion from a Stream Flowing through a Cylindrical Tube, Proceedings of the Royal Irish Academy. Section A: Mathematical and Physical Sciences, 52, (1948-1950), pp. 163-169.

Wagner R., Manninen H.E., Franchin A., Lehtipalo K., Mirme S., Steiner G., Petäjä T. and Kulmala M., On the accuracy of ion measurements using a Neutral cluster and Air Ion Spectrometer, Boreal Environment Research, 21, (2016), pp. 230–241.



