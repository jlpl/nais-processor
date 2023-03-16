# NAIS processor
Code package to process [NAIS](https://www.airel.ee/products/nais/) (Neutral cluster and Air Ion Spectrometer, Airel Ltd.) data files.

The below options are possible:

* Inlet loss correction (Gromley and Kennedy, 1948)
* Ion mode correction (Wagner et al. 2016)
* Conversion to standard conditions (273.15 K, 101325 Pa)
* Remove charger ion band from total particle data
* Use fill values in case of missing environmental sensor data

Some additional utility methods are also included. See [documentation](https://jlpl.github.io/nais-processor/).

## Installation
```shell
pip install nais-processor
```

## Example usage
Use the [`make_config_template()`](https://jlpl.github.io/nais-processor/#nais_processor.make_config_template) method to create a configuration file template and fill it with necessary information. The configuration file is used at processing the data files.

For example:
```
$ python
>>> from nais_processor import *
>>> make_config_template("/home/user/viikki.yml")
```
This will create a configuration file template called `/home/user/viikki.yml`. After filling in the information for our example measurement
the file may look like this:
```yaml
measurement_location: Viikki, Helsinki, Finland
longitude: 25.02
latitude: 60.23
data_folder:
- /home/user/data/2021
- /home/user/data/2022
processed_folder: /home/user/viikki
database_file: /home/user/viikki.json 
start_date: 2022-09-28
end_date: 2022-09-30
inlet_length: 1.0
do_inlet_loss_correction: true
convert_to_standard_conditions: true
do_wagner_ion_mode_correction: true
remove_corona_ions: true
allow_reprocess: false
use_fill_values: true
fill_temperature: 273.15
fill_pressure: 101325.0
fill_flowrate: 54.0
```
Then process the data files by running [`nais_processor()`](https://jlpl.github.io/nais-processor/#nais_processor.nais_processor) method with the config file as the input argument.

In our example case:
```
>>> nais_processor("/home/user/viikki.yml")
Building database...
Processing 20220928 (Viikki, Helsinki, Finland)
Processing 20220929 (Viikki, Helsinki, Finland)
Processing 20220930 (Viikki, Helsinki, Finland)
Done!
```
The code produces daily processed data files `NAIS_yyyymmdd.nc` (netCDF format). These files are saved in the destination given in the configuration file.

The locations of raw and processed files for each day are written in the JSON formatted `database_file`.

## License
This project is licensed under the terms of the GNU GPLv3.

## References
Gormley P. G. and Kennedy M., Diffusion from a Stream Flowing through a Cylindrical Tube, Proceedings of the Royal Irish Academy. Section A: Mathematical and Physical Sciences, 52, (1948-1950), pp. 163-169.

Wagner R., Manninen H.E., Franchin A., Lehtipalo K., Mirme S., Steiner G., Petäjä T. and Kulmala M., On the accuracy of ion measurements using a Neutral cluster and Air Ion Spectrometer, Boreal Environment Research, 21, (2016), pp. 230–241.
