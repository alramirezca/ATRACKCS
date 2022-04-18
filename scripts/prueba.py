# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 18:45:15 2022

@author: ASUS
"""

from atrackcs.utils import funcs
from atrackcs.detect import detect_mcs
from atrackcs.features import features_Tb, features_P, features_Tracks
from atrackcs.track import track_mcs


#________________________________#Global_____________________________________________
"""
UTM_LOCAL_ZONE: Is needed for converting the WGS geodetic coordinate system to plane coordinate system.
This is a constant that must be asociated with the interest region.
For example: 
32718 is the UTM zone 18S plane coordinate system. It was Used for tracking spots in
South America - Colombia

UTC_LOCAL_HOUR: Is needed for converting the raw data hour (UTC) to a local hour (interest region)
For example: 
UTC-5 is the local hour for Colombia. Change the sign (+ or -)   based on the interest local hour in the line: 
datedt = datex.to_pydatetime() (+ or -) timedelta(hours=UTC_LOCAL_HOUR) in function "read_resample_data"
"""

UTM_LOCAL_ZONE = 32718 
UTC_LOCAL_HOUR = 5
UTC_LOCAL_SIGN = "minus"

#___________________________________Functions______________________________________________________


phTb= "C:/Users/ASUS/Desktop/udea/Prueba_GPM/prueba_test_2/tb/"
phP= "C:/Users/ASUS/Desktop/udea/Prueba_GPM/prueba_test_2/p/"

#Tb and P

ds = funcs.readNC(pathTb = phTb, pathP =phP, utc_local_hour = UTC_LOCAL_HOUR, 
                  utc_local_sign = UTC_LOCAL_SIGN)



sup = detect_mcs(ds, detect_scheme = "Both", Tb = 225, area_Tb = 1000, 
                   utm_local_zone = UTM_LOCAL_ZONE, path_save = None)


sup = features_Tb(sup, ds)

sup = features_P(sup, ds, min_precipitation = 2, area_P = 500, drop_empty_precipitation = False)

sup = track_mcs(sup, threshold_overlapping_percentage = None, utm_local_zone = UTM_LOCAL_ZONE,
                path_save = None)

resume = features_Tracks(sup, encrypt_index = True,
                         path_save = None)

funcs.plot_folium(resume, location = [5, -73.94], path_save = r"C:\Users\ASUS\Desktop\udea/")

#Tb 
ds = funcs.readNC(pathTb = phTb, pathP =None, utc_local_hour = UTC_LOCAL_HOUR)



sup = detect_mcs(ds, detect_scheme  = "Tb", Tb = 225, area_Tb = 1000, 
                   utm_local_zone = UTM_LOCAL_ZONE, path_save = None)


sup = features_Tb(sup, ds)

#sup = features_P(sup, ds, min_precipitation = 2, area_P = 500, drop_empty_precipitation = False)

sup = track_mcs(sup, threshold_overlapping_percentage = None, utm_local_zone = UTM_LOCAL_ZONE,
                path_save = None)

resume = features_Tracks(sup, encrypt_index = True,
                         path_save = None)

funcs.plot_folium(resume, location = [5, -73.94], path_save = r"C:\Users\ASUS\Desktop\udea/")

