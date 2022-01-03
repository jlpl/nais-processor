# NAIS processor

Use this code package to process [NAIS](https://www.airel.ee/products/nais/) (Neutral cluster and Air Ion Spectrometer, Airel Ltd.) data files.

The code corrects for diffusion losses in the inlet line (Gromley and Kennedy, 1948) and applies an ion mode calibration (Wagner et al. 2016). Optionally the data can be corrected to standard conditions (273.15 K, 101325 Pa), which can be useful when comparing aerosol particle data from various locations at different altitudes.

## Example using conda 

```
$ git clone https://github.com/jlpl/nais-processor.git
$ cd nais-processor
$ conda env create -f nais_processor.yml
$ conda activate nais-processor
$ (nais-processor) python
>>> from nais_processor import *
>>> nais_processor('config.yml')
```

Where `config.yml` is a configuration file that contains information for processing. A configuration file template is provided.

The code produces daily figures and data files for ion and particle data. These files are saved in the destinations given in the configuration file. 

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

The locations of raw files, processed files and figures as well as possible errors during processing are written in the `database_file`, which is in JSON format. Also a YAML formatted database is written which is easier to read.

## Hints

- For continuous measurements: 
    * Leave the `end_date` out in the configuration file, it will default to current date. 
    * The NAIS creates a new file once a day, so run your processing script also once a day.

## License

This project is licensed under the terms of the GNU GPLv3.

## References

Gormley P. G. and Kennedy M., Diffusion from a Stream Flowing through a Cylindrical Tube, Proceedings of the Royal Irish Academy. Section A: Mathematical and Physical Sciences, 52, (1948-1950), pp. 163-169.

Wagner R., Manninen H.E., Franchin A., Lehtipalo K., Mirme S., Steiner G., Petäjä T. and Kulmala M., On the accuracy of ion measurements using a Neutral cluster and Air Ion Spectrometer, Boreal Environment Research, 21, (2016), pp. 230–241.



