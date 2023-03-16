[![DOI](https://zenodo.org/badge/479858216.svg)](https://zenodo.org/badge/latestdoi/479858216)

![](repo/ATRACKCS_LOGO_.png)


## Introduction

ATRACKCS is a Python package for the automated detection and tracking of Mesoscale convective systems (MCS). MCS are organized cloud clusters that produce regional rainfall and feature vertical development penetrating the mid-upper troposphere. The spatio-temporal characterization of MCS contributed to reducing the vulnerability to severe precipitation events, as well as understanding weather and regional climate. This Python package is a potential tool for characterizing spatio-temporal distribution and evolution of MCS and is intended for researchers and students interested in exploring MCS dynamics. 

ATRACKCS provides a set of functions designed for a workflow analysis that includes the detection and characterization of MCS, as well as the integration in tracks, allowing detailed monitoring of the MCS life cycle both in space and time. The algorithm uses brightness temperature (Tb) and precipitation (P) coming from public satellite data and can operate with Tb as the only input variable or associating precipitation features. The algorithm parameterization can be adapted to the needs of the MCS detection, as the user is allowed to define the thresholds of Tb and P. 

### Technical details and Usage 

The MCS (regions) detection and characterization are performed using these steps: 

1. The input data: Tb and P area processed. The algorithm can operate with Tb as the only input variable. At a given time step, the algorithm finds all pixels where `Tb <= Tb_threshold [200 k, 240 k]` and defines approximate regions with the convex hull, using a binary structure where the pixels that satisfy the described condition are equal to 1 and the remaining pixels are equal to 0. 
2. Transform from geographic to plane coordinates and compute an approximate area of the defined region, discard all regions where area is `<= area_threshold [>= 1000 km^2]` and estimate Tb attributes of those regions.
3. Estimate P attributes of those regions and discard regions based on P rate and P area containing. This is optional.

![](repo/resume_atrackcs_1.png)

The tracks are performed using these steps:

4. overlapping priority principle: for any MCS at time t, the MCS with the highest overlap percentage at time `t+1` "wins" and is associated with it. If there are MCS with lower overlap percentages, they are left to be associated in the next iteration and form a track on their own.
No merging or splitting is allowed, any MCS at time t can only be linked to one MCS at time `t+1`. Similarly, any MCS at time `t+1` can only be linked to one MCS at time `t`. All tracks that do not get updated during the `t` - `t+1` process terminate. 
5. Estimate MCS and tracks attributes.

![](repo/resume_atrackcs_2.png)

The algorithm can be adapted to the needs of the MCS detection, as the user is allowed to define the thresholds of Tb and P (no required). Notebooks and scripts presenting the algorithm functionalities based on use cases are available on the repository page. 

## Input Data 

**Brightness Temperature**: NCEP/CPC L3 (Merge IR V1): Spatial and temporal resolution is 4 km and 30 minutes, data availability from February 7, 2000 to present. The interest variable of this dataset is `Tb` and the file format must be netCDF4. [Access link.](https://doi.org/10.5067/P4HZB9N27EKU)

**Precipitation**: GPM (IMERG V06B): Spatial and temporal resolution is 0.1 degree (∼10 km) and 30 minutes, data availability from June 1, 2000 to present. The interest variable of this dataset is `PrecipitationCal` and the file format must be netCDF4. [Access link.](https://doi.org/10.5067/GPM/IMERG/3B-HH/06)

We suggest the option `subset/get data` and use `OpenDAP` method for downloading and refining the date range and interest region.

## Example use case
Different applications can be made with ATRACKCS. This example tracks a single MCS, but the analyses can be extended both spatially and temporally depending on the interest.

![](repo/example_3.png)

Cloud top Tb obtained from the Merge IR V1 dataset from 08:00-16:00 UTC-5 8th jul 2019 are shown in (a–c). The pixels contained within red contours have Tb less than 225 K. The MCS trajectory that form on 8th jul 2019 was determined using ATRACKCS is shown in (d). The blue and red dots displays the location of the geometric centroid of MCS iniciation and decay.

## Runtime and efficiency

Accurately measuring the runtime and efficiency of ATRACKCS is not a trivial task, since the results may vary significantly due to factors associated with the hardware, software, input region and algorithm parameterization. With the intention of giving a base time estimate, the algorithm was executed for the northern region of South America, composed of 300 pixels (latitude) by 300 pixels (longitude), for the dates between 2019-12-24 19:00 and 2019-12-31 18:00 (168 hours), with a computer equipped with Windows 10 operating system, with an Intel Core i5-8250U processor - 1. 80 GHz and 12 GB RAM; for more details of the algorithm parameterization, please refer the example hosted in the repository page [1_detect_and_track_MCS_scheme_Tb_and_P](https://github.com/alramirezca/ATRACKCS/blob/main/notebooks/1_detect_and_track_MCS_scheme_Tb_and_P.ipynb). This submitted example took 14.5 minutes to be executed and identified 3200 MCS organized in 1207 tracks.

## Main Dependencies

* Python 3.7
* netCDF4 (developed 1.5.6 in py3)
* numpy (developed 1.20.1 in py3)
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
* OS: Linux or Windows, may work in MAC

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

This creates a new environment named `atrackcs_py3`. Activate the environment using:

```
conda activate atrackcs_py3
```

After that, you can check the list of packages installed by:

```
conda list
```

Finally install ATRACKCS using:

```
pip install -e .
```
## tests

To validate installation, issue a new Python session and run:

```
import atrackcs
```

If nothing prints out, installation is successful.

## Inventory

* atrackcs: core module functions.
* notebooks: a series of jupyter notebooks illustrating the major functionalities of the package.
* repo: figures - repository content.
* scripts: example computation scripts. 
* test: integration tests for detecting schemes.


## Changelog

### v1.0

* Can perform MCS detection and tracking through time and space.

## Use in researchs

* Climatological characteristics of deep convection in Northwestern South America using a persistence tracking technique. In American Geophysical Union Fall Meeting 2022.
[ResearchGate](https://www.researchgate.net/publication/366461770_Climatological_characteristics_of_deep_convection_in_Northwestern_South_America_using_a_persistence_tracking_technique).
[AGU fall meeting 2022](https://agu.confex.com/agu/fm22/meetingapp.cgi/Paper/1175671).

* Spatio-temporal Characterization of Mesoscale Convective Systems over Northern South America. In American Geophysical Union Fall Meeting 2021.
[AGU fall meeting 2021](https://agu.confex.com/agu/fm21/meetingapp.cgi/Paper/874852).

* Cloud-resolving Simulations of Mesoscale Convective Systems in Colombia. In American Geophysical Union Fall Meeting 2021.
[ResearchGate](https://www.researchgate.net/publication/357975142_Cloud-resolving_Simulations_of_Mesoscale_Convective_Systems_in_Colombia).
[AGU fall meeting 2021](https://agu.confex.com/agu/fm21/meetingapp.cgi/Paper/875417).

## Cite as
Ramírez-Cardona Álvaro, Robledo Vanessa, Rendón A. Angela M., Henao Juan. J, Hernández, K. Santiago, Gómez-Ríos Sebastián,  & Mejía John F. (2022). Algorithm for Tracking Convective Systems (ATRACKCS) (v1.0). Zenodo. https://doi.org/10.5281/zenodo.7025990

## Acknowledgements

This research was funded by the Colombian Ministry of Science, Technology and Innovation (MINCIENCIAS) through the project: “Implementación de un sistema de investigación y pronóstico meteorológico de corto plazo con el modelo WRF, para apoyo a sistemas de comando y control de la Fuerza Aérea Colombiana” (code 1115-852-70955) with funds from “Patrimonio Autónomo Fondo Nacional de Financiamiento para la Ciencia, la Tecnología y la Innovación, Fondo Francisco José de Caldas”. 