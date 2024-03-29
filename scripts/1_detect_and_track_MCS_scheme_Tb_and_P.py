# -*- coding: utf-8 -*-
"""
Track Mesoescale Convective Systems (MCS) at consecutive time steps to
form tracks with a scheme based on brightness temperature and precipitation rates

@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
          
This script tracks detected MCS at individual time steps to form tracks.

ATRACKCS is a Python package for the automated detection and tracking of MCS. 
It is a potential tool for characterizing their spatio-temporal distribution 
and evolution. ATRACKCS provides a set of Python functions designed for a 
workflow analysis that includes the detection and characterization of MCS, 
as well as the integration in tracks, allowing detailed monitoring of the MCS
life cycle both in space and time.

ATRACKCS uses brightness temperature (Tb) and precipitation (P) coming from
satellite data and can operate with Tb as the only input variable or associating
precipitation features. Although the magnitude of Tb does not represent a
direct measurement of P, it is an indirect representation of cloud cover 
height associated with an MCS event. Some methodologies for MCS detection use 
other satellite spectral bands as a proxy for P. The algorithm parameterization 
can be adapted to the needs of the MCS detection, as the user is allowed to 
define the thresholds of Tb and P. The parameterization in this notebook was 
established from an extensive literature review, which can be consulted below.        

## The detection of the MCS (regions) is performed using these steps:

1. At a given time step, the algorithm finds all pixels where Tb  ≤225𝐾  and defines
approximate regions with the convex hull, using a binary structure where the pixels
that satisfy the described condition are equal to  1  and the remaining pixels are
equal to  0.
2. Transform from geographic to plane coordinates and compute an approximate area 
of the defined regions.
3. Discard all regions whose area is  <2000𝑘𝑚2 .
4. Estimate the average, minimum and maximum brightness temperature of those regions.
5. At a given time step, the algorithm finds all pixels where P  ≥2𝑚𝑚×ℎ−1  and discard
pixels that do not match with this condition.
6. Estimate an approximate area of the precipitation pixels that are contained in each region.
7. Estimate the average and maximum precipitation for each region whose area is  ≥500𝑘𝑚2 .
The algorithm has the option of discard those regions whose precipitaion area is  <500𝑘𝑚2 , 
but in this case those regions are going to be part of the possible tracks.

## The tracks are performed using these steps:

Specifically, assume we have detected  𝑛  MSC at time  𝑡 , and  𝑚  MSC at time  𝑡+1.
There are theoretically  𝑛×𝑚  possible associations to link these two groups of MCS.
Of course not all of them are meaningful. The rules that are applied in the association
process are:

1. Overlapping priority principle: for any MCS at time t, the MCS with the highest overlap
percentage at time t+1 is associated with it.
2. The MCS (with lower or no overlap percentages) at time  𝑡+1  could form a track on their own, 
and are left to be associated in the next iteration between  𝑡+1  and  𝑡+2 .
3. No merging or splitting is allowed, any MCS at time  𝑡  can only be linked to one MCS at 
time 𝑡+1 , similarly, any MCS at time  𝑡+1  can only be linked to one MCS at time 𝑡 .
4. All tracks that do not get updated during the  𝑡  -  𝑡+1  process terminate. This assumes 
that no gap in the track is allowed.
5. In this parameterization, no tracks are discarded based on their total duration. 
The algorithm has the option for filtering the tracks with a specific minimun duration.

## Input data

**Brightness Temperature:**
NCEP/CPC L3 (Merge IR V1): Spatial and temporal resolution is 4 km and 30 
minutes, data availability from February 7, 2000 to present. The interest 
variable of this dataset is `Tb` and the files format must be `netCDF4`.
https://doi.org/10.5067/P4HZB9N27EKU

**Precipitation:**
GPM (IMERG V06B): Spatial and temporal resolution is 10 km and 30 minutes, 
data availability from June 1, 2000 to present. The interest variable of 
this dataset is `PrecipitationCal` and the files format must be `netCDF4`.
https://doi.org/10.5067/GPM/IMERG/3B-HH/06

We suggest the option subset/get data and use OpenDAP method for downloading
and refining the date range and interest region.

In this case the algorithm is run for 5 days 
(2019/12/25/ 00 (UTC) - 2019/12/31/ 00 (UTC)) 
for northern South America. The input data are in the notebooks folder in the
repository. The raw data can be downloaded with the links at the top.

## Steps

1. Make sure you have successfully installed the ATRACKCS library.
2. Execute the following code blocks in sequence.

## Results

* `resume_Tb_P_2019_12_24_19_2019_12_31_17.csv`: a csv table listing various
 attributes for the tracks and MCS associated.
* `map_Tb_P_2019_12_24_19_2019_12_31_10.html` (folium): plot of the 
geographical locations of the MSC with informations that links to the associaded tracks and features of the MSC.

## Bibliography
* Feng, Z., Leung, L. R., Liu, N., Wang, J., Houze, R. A., Li, J., Hardin, J. C., Chen, D., & Guo, J. (2021). A Global High‐resolution Mesoscale Convective System Database using Satellite‐derived Cloud Tops, Surface Precipitation, and Tracking. Journal of Geophysical Research: Atmospheres. https://doi.org/10.1029/2020jd034202
* Li, J., Feng, Z., Qian, Y., & Leung, L. R. (2020). A high-resolution unified observational data product of mesoscale convective systems and isolated deep convection in the United States for 2004–2017. Earth System Science Data Discussions, October, 1–48. https://doi.org/10.5194/essd-2020-151
* Liu, W., Cook, K. H., & Vizy, E. K. (2019). The role of mesoscale convective systems in the diurnal cycle of rainfall and its seasonality over sub-Saharan Northern Africa. Climate Dynamics, 52(1–2), 729–745. https://doi.org/10.1007/s00382-018-4162-y
* Vizy, E. K., & Cook, K. H. (2018). Mesoscale convective systems and nocturnal rainfall over the West African Sahel: role of the Inter-tropical front. Climate Dynamics, 50(1–2), 587–614. https://doi.org/10.1007/s00382-017-3628-7    

"""
#________________________________#Set Paths_____________________________________________
#AFirst of all we assign the locations to the input (raw data) using `TBDIR` and `PDIR` 
#and output data using `OUTDIR`.

TBDIR=  r"C:/Users/ASUS/Documents/GitHub/atrackcs/notebooks/1_input_data/tb/"
PDIR= r"C:/Users/ASUS/Documents/GitHub/atrackcs/notebooks/1_input_data/p/"

OUTDIR=r"C:/Users/ASUS/Documents/GitHub/atrackcs/scripts/1_output_data/"

#_____________#Parameters used in the MCS detection and tracking process_____________________________________________
"""
* UTM_LOCAL_ZONE: int, is needed for converting the WGS geodetic coordinate 
  system to plane coordinate system. This is a constant that must be asociated
  with the interest region.

* UTC_LOCAL_HOUR: int, is needed for converting the raw data hour (UTC) to a
  local hour (interest region).

* UTC_LOCAL_SIGN: str (minus, plus, local), is needed for converting the raw
  data hour (UTC) to a local hour (interest region).

* DETECT_SCHEME: str (Tb, Both), association scheme: Tb or Both (Tb and P).

* TB: int, MCS based on limited maximun threshold cold cloud top.

* AREA_TB: int, MCS with a minimun largest area.

* MIN_P: int, MCS with a minimun precipitation pixel.

* AREA_P: int, MCS with an mimimum area precipitation.

* DROP_EMPTY_PRECIPITATION: boolean, if True eliminates MCS that do not contain
  precipitation with the MIN_P and AREA_P selected.

* THRESHOLD_OVERLAPPING_P: int, percentage overlap limit between MCS to be considered part
  of the track.

* LOCATION_FOLIUM: list (lat, lon), location for center the map_folium.

* MIN_DURATION: int, for filtering tracks based on minimun duration.
"""

#32718 is the UTM zone 18S plane coordinate system. 
#It was Used for tracking MCS in South America - Colombia.
UTM_LOCAL_ZONE = 32718 

#UTC-5 is the local hour for Colombia
UTC_LOCAL_HOUR = 5
UTC_LOCAL_SIGN = "minus"

#Scheme of association
DETECT_SCHEME = "Both"

TB = 225 #[Feng et al.,(2021); Li et al.,(2020)]

AREA_TB = 2000 # [Lui et al., (2019); Vizy & Cook,(2018)] 

MIN_P = 2 #[Feng et al.,(2021);Li et al.,(2020)]  

AREA_P = 500 #[Schumacher et al.,(2021)]

DROP_EMPTY_PRECIPITATION = False

THRESHOLD_OVERLAPPING_P = 0

LOCATION_FOLIUM = [5, -73.94]

MIN_DURATION = 0
#________________________________#Import modules_____________________________________________

from atrackcs.utils import funcs
from atrackcs.detect import detect_mcs
from atrackcs.features import features_Tb, features_P, features_Tracks
from atrackcs.track import track_mcs
import os

#_________________________#Reading raw Tb and P data._____________________________________________
"""
The read process is handled with the function funcs.readNC(). Also make sure
the output folder exists.

* ds is the xarray.Dataset object we just read in.
"""
ds = funcs.readNC(pathTb = TBDIR, pathP =PDIR, utc_local_hour = UTC_LOCAL_HOUR,
                  utc_local_sign = UTC_LOCAL_SIGN)

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
#___________________________________Detecting MCS______________________________________________________
"""
The detection process is handled with the function `detect_mcs()`.
The MCS Tb features is handled with the function features_Tb().
The MCS P features is handled with the function `features_P()`.

* mcs_l is the Geopandas.GeoDataFrame object we just created in the detection
 process.
 
Each `row` object in `mcs_l` stores a MCS record. We can have a peak into what
`msc_l` columns contains so far. 

* `time`: datetime64, hour (UTC-5) for this case.
* `Tb`: float, polygon index after being filtered by the established
  parameterization.
* `geometry`: geometry, MCS polygon (convex hull). To view the coordinate
  reference system of the geometry column, access the crs attribute. 
  For this case is `EPSG:32718`.
* `area_Tb`: float, area polygon [km^2]
* `centroid_`: geometry, geometric centroid polygon (convex hull). 
  The crs in  this case is `EPSG:32718`.
* `mean_tb`: float, brightness temperature average of the pixels composing
  the polygon. [K]
* `min_tb`: float, brightness temperature min value of the pixels composing
  the polygon. [K]
* `max_tb`: float, brightness temperature max value of the pixels composing
   the polygon. [K] 
"""
mcs_l = detect_mcs(ds, detect_scheme = DETECT_SCHEME, Tb = TB, 
                   area_Tb = AREA_TB, 
                   utm_local_zone = UTM_LOCAL_ZONE, path_save = None)

mcs_l = features_Tb(mcs_l, ds)

mcs_l = features_P(mcs_l, ds, min_precipitation = MIN_P, area_P = AREA_P,
                   drop_empty_precipitation = DROP_EMPTY_PRECIPITATION)

#_______________________________Tracking MCS______________________________________________________
"""
The tracking process is handled with the function track_mcs(). 
However, these have not been characterized in this phase. 
The tracks features is handled with the function features_Tracks()

* tracks_l is a geopandas.GeoDataFrame object we just created in the tracking
  process and stores the MCS polygons and the tracks.
* This GeoDataFrame contains exactly the information previously referenced 
  except for the dropping of some features associated to the brightness
  temperature and the addition of some characteristics associated to the 
  MCS precipitation and some characteristics generated for the tracks.
* The indexing of this object is a pandas.MultiIndex since it hierarchically
  associates an id for each track and an id for each MCS.
* The indexing of this object is encrypted generated for each track using the
  uuid library. To disable encryption for both MCS and tracks, use the parameter
  encrypt_index = False in the features_Tracks() function. The index generated
  when not using encryption is a int iteration result of the tracking process.
* The indexes encryption is useful when processing a long period of time and
  must iterate for smaller periods of time, for example every 6 months. 
  This is a limitation imposed by the hardware used to run the algorithm.

We can have a peak into what tracks_l new columns contains.

* belong: str, encrypted index generated for each track.
* id_gdf: str, encrypted index generated for each MCS.
* geometry: geometry, polygon. The crs is EPSG:4326-WGS84
* centroid_: geometry, geometric centroid polygon. The crs is EPSG:4326.
* mean_pp: float, precipitation average of the pixels composing the polygon
  [mm/h].
* max_pp: float, precipitation max value of the pixels composing the polygon
  [mm/h].
* intersection_percentage: float, percentage overlap between MCS [%].
* distance_c: float, distance between the overlapping geometric centroids of 
  the MCS [km].
* direction: float, direction between the overlapping geometric centroids of
  the MCS [°].
* total_duration: float, total duration of the event or the track. This value 
   is associated to each MCS of the corresponding track [h].
* total_distance: float, total distance of the event or the track. This value
  is associated to each MCS of the corresponding track [km].
* mean_velocity: float, velocity average of the event or the track. This value 
  is associated to each MCS of the corresponding track [km/h]. 
"""
tracks_l = track_mcs(mcs_l, threshold_overlapping_percentage = 
            THRESHOLD_OVERLAPPING_P, utm_local_zone = UTM_LOCAL_ZONE,
                     path_save = None)

tracks_l = features_Tracks(tracks_l, initial_time_hour = MIN_DURATION,
                         path_save = OUTDIR)

print ("total tracks detected: " + str(len(tracks_l.index.levels[0])) + 
       " and total MCS detected: " + str(len(tracks_l.index.levels[1])))

#_______________________________Plotting MCS______________________________________________________
"""
As the last step, let is reload the tracks results from the local storage 
and plot the MCS with the help of the folium library

The read process tracks is handled with the function `funcs.readTRACKS()`.
The function `funcs.plot_folium()` saves the `.html` result and return the path where was saved.
"""

#-------------------Load results-------------------
pathload = OUTDIR + "resume_Tb_P_2019_12_24_19_2019_12_31_17.csv"

tracks_l = funcs.readTRACKS(pathload)

#---------Discarding tracks that last less than 10 hours--------
tracks_l = tracks_l.reset_index()
tracks_l = tracks_l[tracks_l.total_duration >= 10].sort_index()
tracks_l = tracks_l.set_index(["belong", "id_gdf"])
tracks_l = tracks_l.sort_index()

print ("total tracks that last less than 10 hours : " + 
       str(len(tracks_l.index.levels[0])) + " and total MCS detected: " 
       + str(len(tracks_l.index.levels[1])))

#---------Generating plot map folium--------
path_html_folium = funcs.plot_folium(tracks_l, location = LOCATION_FOLIUM, 
                                     path_save = OUTDIR)


