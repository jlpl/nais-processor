# NAIS processor

Use this code package to process [NAIS](https://www.airel.ee/products/nais/) (Neutral cluster and Air Ion Spectrometer, Airel Ltd.) measurements.

The code corrects for diffusion losses in the inlet line (Gromley and Kennedy, 1948) and applies an ion mode calibration (Wagner et al. 2016). Optionally the data can be corrected to standard conditions (273.15 K, 101325 Pa).

## Example usage


```
$ git clone https://github.com/jlpl/nais-processor.git
$ cd nais-processor
$ conda env create -f nais_processor.yml
$ mkdir processed figures
$ touch config.yml process.py
$ conda activate nais-processor
$ (nais-processor) python process.py
```

Files:

```YAML
# config.yml

database_file: './db.json'
inlet_length: 1.0 # in meters
location: 'My measurement site'
data_folder: '/path/to/raw/nais/data/'
processed_folder: './processed/'
figure_folder: './figures/'
start_date: '2018-01-01'
end_date: '2018-12-31'
sealevel_correction: False
```

```python
# process.py

from nais_processor import *
nais_processor('config.yml')
```

## Hints

- For continuous measurements: 
    * Leave `end_date` out, it will default to current date. 
    * The NAIS creates a new file once a day, so run `process.py` once a day.
- The locations of raw files, processed files and figures as well as possible errors during processing are written in the `database_file`. Inspect the database with `jq` and query with `tinydb`.


## License

This project is licensed under the terms of the GNU GPLv3.

## References

Gormley P. G. and Kennedy M., Diffusion from a Stream Flowing through a Cylindrical Tube, Proceedings of the Royal Irish Academy. Section A: Mathematical and Physical Sciences, 52, (1948-1950), pp. 163-169.

Wagner R., Manninen H.E., Franchin A., Lehtipalo K., Mirme S., Steiner G., Petäjä T. and Kulmala M., On the accuracy of ion measurements using a Neutral cluster and Air Ion Spectrometer, Boreal Environment Research, 21, (2016), pp. 230–241.



