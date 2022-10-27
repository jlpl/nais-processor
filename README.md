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
Use the [`make_config()`](https://jlpl.github.io/nais-processor/#nais_processor.make_config) method to create a configuration file that is used at processing the data files. You will be asked questions about the measurement and the data. The information is written into a yaml file.
```
$ python
>>> from nais_processor import *
>>> make_config()
```
The resulting configuration file may look for example like this:
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
remove_corona_ions: true
remove_noisy_electrometers: true
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

where `n` and `p` refer to negative and positive polarity respectively. `yyyymmdd` tells the date in the year-month-day format. `np` and `nds` refer to particle and ion data respectively.

In the processed data files the header contains the geometric mean diameters of the size bins, the first column is the time and the rest of the data is the number-size distribution matrix with normalized number concentrations (dN/dlogDp). 

The locations of raw files, processed files and cleaned processed files are written in the JSON formatted `database_file`.

### Combining sumfiles
Once you have processed your NAIS data you can extract any time range in a sum matrix format using the [`combine_spectra()`](https://jlpl.github.io/nais-processor/#nais_processor.combine_spectra) function.

Example:
```python
import nais_processor as nais

start_time="2022-09-29 02:00+02:00"
end_time="2022-09-30 14:00+02:00"

combined_data = nais.combine_spectra(
    config_file,start_time,end_time,spectra_type="negion")
```

Use the [`aerosol-functions`](https://github.com/jlpl/aerosol-functions) package to further analyze and plot the data.

## License
This project is licensed under the terms of the GNU GPLv3.

## References
Gormley P. G. and Kennedy M., Diffusion from a Stream Flowing through a Cylindrical Tube, Proceedings of the Royal Irish Academy. Section A: Mathematical and Physical Sciences, 52, (1948-1950), pp. 163-169.

Wagner R., Manninen H.E., Franchin A., Lehtipalo K., Mirme S., Steiner G., Petäjä T. and Kulmala M., On the accuracy of ion measurements using a Neutral cluster and Air Ion Spectrometer, Boreal Environment Research, 21, (2016), pp. 230–241.
