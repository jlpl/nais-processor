# NAIS processor

Use this code package to process NAIS (Neutral cluster and Air Ion Spectrometer) measurements.

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
nais_plotter('config.yml')
```

## Hints

- For continuous measurements choose `end_date` in distant future. The NAIS creates a new file once a day, so run `process.py` once a day.
- The locations of raw files, processed files and figures as well as possible errors during processing are written in the `database_file`. Inspect the database with `jq` and query with `tinydb`.

## License

This project is licensed under the terms of the GNU GPLv3.


