# NAIS processor

Use this code package to process [NAIS](https://www.airel.ee/products/nais/) (Neutral cluster and Air Ion Spectrometer, Airel Ltd.) data files.

The code corrects for diffusion losses in the inlet line (Gromley and Kennedy, 1948) and applies an ion mode calibration (Wagner et al. 2016). Optionally the data can be corrected to standard conditions (273.15 K, 101325 Pa), which can be useful when comparing aerosol particle and ion data from various locations at different altitudes.

## Installation
```
pip install nais-processor
```

## Example usage

Open the python prompt and load methods from the `nais_processor` module.
Then use the `make_config()` method to create a configuration file that
is used at processing the data files.

```
$ python
>>> from nais_processor import *
>>> make_config()

Enter name of configuration file.
E.g. ./configs/campaign.yml
> ./configs/nyc.yml

Give path(s) to raw data. If multiple paths give them as comma separated list.
E.g. /data/nais/2021,/data/nais/2022
> /campaigns/nyc/nais/2015,/campaigns/nyc/nais/2016

Path to where processed data is saved.
E.g. ./data/campaign/processed
> ./data/nyc/processed

Path to where figures are saved. Leave empty if no figures.
E.g. ./data/campaign/figures
> ./data/nyc/figs

Start of measurement (YYYY-MM-DD)
> 2015-01-05

End of measurement (YYYY-MM-DD)
If empty processor assumes current day, use for continuous measurements.
> 2016-06-08

Enter name of database file
E.g. ./logs/campaign.json
> ./logs/nyc.json

Measurement location
E.g. Helsinki, Kumpula
> New York City

Apply corrections to data? (True/False)
Requires a NAIS with temperature and pressure sensors.
> True

Length of the inlet in meters
> 1.0

Correct concentrations to sealevel conditions? (True/False)
> False

Configuration saved to: ./configs/nyc.yml
```
Then process the data files by running `nais_processor()` method with the config file as the input argument.

```
>>> nais_processor("./configs/nyc.yml")
Configuration file: ./configs/nyc.yml
processing 20150105
processing 20150106
processing 20150107
...
Done!
```
Run `do_figs()` with the config file in order to create plots of the processed data if needed.

```
>>> do_figs("./configs/nyc.yml")
plotting 20150105
plotting 20150106
plotting 20150107
...
Done!
```
The code produces daily processed data files and optionally figures for ion and particle data. These files are saved in the destinations given in the configuration file.

The data files are named

`NAIS[n|p][yyyymmdd][np|nds].sum`

where `n` and `p` refer to negative and positive polarity respectively. `yyyymmdd` tells the date in the year-month-day format. `np` and `nds` refer to particle and ion data respectively.

The data files have the following structure:

```
[0,0]  = UTC offset in hours
[1:,0] = Time (MATLAB datenum) 
[0,2:] = Geometric mean diameter of size-channel (m)
[1:,1] = Integrated total number concentration (cm-3)
[1:,2:] = Normalized number concentrations, dN/dlogDp (cm-3)
```

The locations of raw files, processed files and possible figures are written in the `database_file`, which is in JSON format.

## License

This project is licensed under the terms of the GNU GPLv3.

## References

Gormley P. G. and Kennedy M., Diffusion from a Stream Flowing through a Cylindrical Tube, Proceedings of the Royal Irish Academy. Section A: Mathematical and Physical Sciences, 52, (1948-1950), pp. 163-169.

Wagner R., Manninen H.E., Franchin A., Lehtipalo K., Mirme S., Steiner G., Petäjä T. and Kulmala M., On the accuracy of ion measurements using a Neutral cluster and Air Ion Spectrometer, Boreal Environment Research, 21, (2016), pp. 230–241.



