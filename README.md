# Algorithm for Tracking Convective Systems (ATRACKCS)

## Introduction

ATRACKCS (Algorithm for Tracking Convective Systems) is a Python package for the automated detection and tracking of mesoscale convective systems (MCS). The algorithm uses brightness temperature (Tb) and precipitation (P) coming from satellite data.

### The detection of the MCS (regions) is performed using these steps:

1. At a given time step, the algorithm finds all pixels where `Tb <= Tb_threshold [200 k - 240 k]` and defines approximate regions with the convex hull,  using a binary structure where the pixels that satisfy the described condition are equal to 1 and the remaining pixels are equal to 0.
2. Transform from geographic to plane coordinates and compute an approximate area of the defined regions. 
3. Discard all regions whose area is `<= area_threshold [>= 1000 km**2]`
4. Estimate Tb attributes of those regions.
5. Estimate P attributes of those regions. This is optional as the algorithm can operate only with Tb as input variable.

![](joss/resume_atrackcs_1.png)

### The tracks are performed using these steps:

1. **overlapping priority** principle: for any MCS at time `t`, the MCS with the highest olverlap percentage at time `t+1` "wins" and is associated with it. 
2. The MCS (with lower or no overlap percentages) at time `t+1` could form a track on their own, and are left to be associated in the next iteration between `t+1` and `t+2`.
3. No merging or splitting is allowed, any MCS at time `t` can only be linked to one MCS at time `t+1`. Similarly, any MCS at time `t+1` can only be linked to one MCS at time `t`.
4. All tracks that do not get updated during the `t` - `t+1` process terminate. This assumes that no gap in the track is allowed. 

![](joss/resume_atrackcs_2.png)

 The algorithm can be adapted to the needs of the MCS detection, as the user is allowed to define the thresholds of Tb and P (no required). ATRACKCS is intended for researchers and students who are interested in the characterization of MCS both in the meteorological and climatological fields.

## Input Data 

**Brightness Temperature**: NCEP/CPC L3 (Merge IR V1): Spatial and temporal resolution is 4 km and 30 minutes, data availability from February 7, 2000 to present. The interest variable of this dataset is `Tb` and the file format must be netCDF4. [Access link.](https://doi.org/10.5067/P4HZB9N27EKU/)

**Precipitation**: GPM (IMERG V06B): Spatial and temporal resolution is 0.1 degree (∼10 km) and 30 minutes, data availability from June 1, 2000 to present. The interest variable of this dataset is `PrecipitationCal` and the file format must be netCDF4. [Access link.](https://doi.org/10.5067/GPM/IMERG/3B-HH/06)

We suggest the option `subset/get data` and use `OpenDAP` method for downloading and refining the date range and interest region.

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
* folium (developed in in 0.12.1.post1 in py3)
* OS: Linux or Windows, may work in MAC.

## Installation

It is recommended to build the Python environment using [Anaconda](https://www.anaconda.com/distribution/). Installation may take several minutes.

### Create conda environment using environment file

After Anaconda installation, git clone this repository:

```
git clone https://github.com/alramirezca/ATRACKCS
```

Then build a new conda environment using the environment file provided. For example:

```
cd ATRACKCS
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

* atrackcs: core module functions.
* notebooks: a series of jupyter notebooks illustrating the major functionalities of the package.
* scripts: example computation scripts. 
* test: integration tests for detecting schemes.

## Changelog

### v1.0

* initial upload. Can perform MCS detection and tracking through time and space.

## Use in researchs

* Spatio-temporal Characterization of Mesoscale Convective Systems over Northern South America. In AGU Fall Meeting 2021.
[AGU fall meeting 2021](https://agu.confex.com/agu/fm21/meetingapp.cgi/Paper/874852).

* Cloud-resolving Simulations of Mesoscale Convective Systems in Colombia. In AGU Fall Meeting 2021.
[ResearchGate](https://www.researchgate.net/publication/357975142_Cloud-resolving_Simulations_of_Mesoscale_Convective_Systems_in_Colombia).
[AGU fall meeting 2021](https://agu.confex.com/agu/fm21/meetingapp.cgi/Paper/875417).

