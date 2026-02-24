# `nais-processor`
This code project is designed to facilitate [NAIS](https://www.airel.ee/products/nais/) (Neutral cluster and Air Ion Spectrometer, Airel Ltd.) and [CIC](https://www.airel.ee/products/cic/) (Cluster Ion Counter, Airel Ltd.) data analysis.

Latest version: 0.2.0

## Installation

Install from GitHub using `pip`
```shell
pip install git+https://github.com/jlpl/nais-processor.git
```

## Documentation
See [here](https://jlpl.github.io/nais-processor/)

## `nais` package

See [here](https://jlpl.github.io/nais-processor/tutorial.html#) for a tutorial on data analysis work flow.

### `nais.processor` module
Generate easy-to-handle netcdf data files from any NAIS instrument generation. Apply optional corrections and transforms.

### `nais.utils` module
Utility tools to work with the processed files.

### GUI app: `nais-data-checker`
GUI tool to visually inspect the data together with error flags and define regions of bad data.

## `cic` package

#### `cic.processor` module
Generate similar netcdf data files from CIC data. Apply optional corrections and transforms.

#### `cic.utils` module
Utility tools to work with the processed files from CIC.

## License
This project is licensed under the terms of the GNU GPLv3.
