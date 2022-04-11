# Algorithm for Tracking Convective Systems (ATRACKCS)

## Introduction

ATRACKCS (Algorithm for Tracking Convective Systems) is a Python package for automated mesoscale convective systems (MCS) detection and tracking. The operation of the algorithm is by using brightness temperature and precipitation, whose data coming from satellite data (MERGE IR and IMERG V06B). 

The detection of the MCS is from the cold top of the clouds, based on magnitude threshold, and the generation of an approximate horizontal area from the convex hull. The tracking in time and space is done by overlapping areas. The algorithm allows parameterization and can be adapted to the specific needs of each geographic environment and/or SCM detection need.

ATRACKCS is intended for researchers and students who are interested in the characterization of MCS both in the meteorological (short term) and climatological (long term) fields.

## Documentation



## Example use case



## Main Dependencies

* Python3.7.
* netCDF4 (developed 1.5.6 in py3)
* numpy (tdeveloped 1.20.1 in py3)
* pandas (developed 1.3.5 in py3)
* scipy (developed 1.7.3 in py3)
* geopandas (developed in 0.9.0 in py3)
* xarray (developed in 0.18.0 in py3)
* gdal (developed in 3.1.4 in py3)
* geopy (developed in 2.2.0 in py3)
* shapely (developed in 1.8.0 in py3)
* rioxarray (developed in 0.9.1 in py3)
* rasterio (developed in 1.2.1 in py3)
* folium (optional, only used for plotting in 0.12.1.post1 in py3)
* OS: Linux or Windows, may work in MAC.

## Installation

Recommend building the Python environment using [Anaconda](https://www.anaconda.com/distribution/).

### Create conda environment using environment file

This way will install the optional `cartopy` package and allow you to run
the notebook examples.

After Anaconda installation, git clone this repository:

```
git clone https://github.com/alramirezca/atrackcs
```

Then build a new conda environment using the environment file provided. For example:

```
cd IPART
conda env create -f env_py3.yml
```

This creates a new environment named `atrackcs_py3`. Activate the environment using

```
conda activate atrackcs_py3
```

After that, you can check the list of packages installed by

```
conda list
```

Finally install ATRACKCS using:

```
pip install -e .
```


## tests

To validate installation, issue a new Python session and run

```
import atrackcs
```

If nothing prints out, installation is successful.



## Inventory

* docs: readthedocs documentation.
* atrackcs: core module functions.
* notebooks: a series of jupyter notebooks illustrating the major functionalities of the package.
* scripts: example computation scripts. Can be used as templates to quickly develop your own working scripts.


## Changelog

### v1

Minor fixes:

* fix a bug in latitudinal range filtering when data cover both of the Northern and Southern Hemispheres.
* more robust handling of zonally cyclic data.
* (related to a change in v3.3.0) a better way to prevent potential [matplotlib memory leaking](https://github.com/matplotlib/matplotlib/issues/20490).


### v1.0

* initial upload. Can perform MSC detection and tracKing through time and space.

## Contributing

We welcome contributions from the community. Please create a fork of the project on GitHub
and use a pull request to propose your changes. 

## Citation

If you use `ATRACKCS` in published research, please cite it by referencing the

## Getting help

Please post issues on the project GitHub page.
