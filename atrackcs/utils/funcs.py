# -*- coding: utf-8 -*-
"""
Utility functions
@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
"""
from os import path
import xarray as xr  
import numpy as np
import pandas as pd
from scipy import ndimage
import geopandas as gpd
from datetime import timedelta
import rioxarray
import rasterio
from geopy.distance import geodesic
import math
import sys
import glob
import numbers as checknumbers
from shapely.geometry import MultiPolygon, Polygon, shape, Point, MultiPoint, mapping 
from shapely.wkt import loads
import uuid
import tqdm
import time
import warnings
import folium
import webbrowser
warnings.filterwarnings("ignore")


#___________________________________Functions______________________________________________________

def readNC(pathTb = None, pathP = None, utc_local_hour = 0, utc_local_sign = "minus"):
    """
    Function for reading and resampling the Tb and P DataArrays. 
    The spatial resampling is 0.1° - lineal interpolation.
    The temporal resampling is 1 h - nearest original coordinate to up-sampled frequency coordinates.
    
    Inputs:
    * pathTb: str, path where the Tb raw data are located.
    * pathP: str, path where the P raw data are located.
    The path must have the next structure:
    linux: r"/home....../"
    windows: r"C:/....../"
    
    * utc_local_hour: int, allows transform the raw data hour (UTC) 
    to a interest time zone (interest region). 
    * utc_local_sign: str (minus, plus, local), sets whether to add or subtract 
    in the conversion for the time zone of interest. If set "local" no conversion will be done
    and the time zone will be GMT or UTC.
    (ex: UTC-5 timing is determined by utc_local_hour = 5 and utc_local_sign = "minus".
    
    Outputs:
    * xarray.Dataset with the brightness temperature (Tb) and P (Precipitation) data.
    """
      
    if isinstance(pathTb, str) and isinstance(pathP, str): 
        try:
            #Globbing the Tb and P files
            filenamestb = glob.glob(pathTb+'*.nc4')
            filenamesp =glob.glob(pathP+'*.nc4')
            
            #Reading P data
            ds_p = xr.open_mfdataset(filenamesp)
            ds_p['P'] = ds_p['precipitationCal']; del ds_p['precipitationCal']
            ds_p = ds_p['P'].T
            #Temporal resampling precipitation data
            ds_p  = ds_p.resample(time ='1H').mean()  #promedio Horario de precipitacion
            
            #Reading Tb data
            ds_t = xr.open_mfdataset(filenamestb); ds_t = ds_t['Tb']
            #Temporal resampling Tb data
            ds_t  =  ds_t.resample(time="1H").nearest(tolerance="1H") 
            
            #Spatial resampling Tb DataArray. This is based on P coordinates (lat and lon).
            ds_t =  ds_t.interp(lat=ds_p.lat.values, lon=ds_p.lon.values)
            #Reorder levels from Tb DataArray for merging with P DataArray.
            ds_t = ds_t.transpose("lat", "lon", "time")
            
            #Merging DataArrays
            try:
                ds = xr.merge([ds_p, ds_t]) 
            except:
            #The raw P data from some years does not have a valid datetime index 
            #These lines convert CFTimeIndex to DatetimeIndex for merging.
                datetimeindex = ds_p.indexes['time'].to_datetimeindex()
                ds_p['time'] = datetimeindex
                ds = xr.merge([ds_p, ds_t]) 
        
            #Closing the raw DataArrays
            ds_p.close()
            ds_t.close()
            
            #Converting the UTC to local hour
            datex = ds.time.coords.to_index()
            
            #Replacing the datetimeindex based on UTC_LOCAL_HOUR
            if utc_local_sign == "minus":
                datedt = datex.to_pydatetime() - timedelta(hours=utc_local_hour)
            elif utc_local_sign == "plus":
                datedt = datex.to_pydatetime() + timedelta(hours=utc_local_hour)
            elif utc_local_sign == "local":
                datedt = datex.to_pydatetime() 
            else:
                raise TypeError("You must type a valid parameter for utc_local_sign: minus, plus or local. If you use local please enter utc_local_hour = 0")

            dates_64 = [np.datetime64(row) for row in datedt]
            ds =  ds.assign_coords({"time": dates_64})
                
            #Attaching Atributes to DataArray merged.
            ds.Tb.attrs["units"] = "K"
            ds.P.attrs["units"] = "mm/h"
            ds.Tb.attrs["_FillValue"] = 'NaN'
            ds.P.attrs["_FillValue"] = 'NaN'
            ds.lon.attrs['units'] = "degrees_east"
            ds.lat.attrs['units'] = "degrees_north"
            
            #Extracting dimensiones: time, lat and lon
            dates = ds.time.values; 
            lon, lat = np.float32(np.meshgrid(ds.lon,ds.lat))
         
            #Establishing EPSG:4326
            ds = ds.rio.set_crs(4326)
            ds.attrs['crs'] = ds.rio.crs
            
            initial_date_lecture = str(ds.time[0].values)[:16]
            final_date_lecture = str(ds.time[-1].values)[:16]
            
            print('Complete Tb and P data reading ' + initial_date_lecture + " - " + final_date_lecture)    
        except:
            raise FileNotFoundError("Make sure you are complying with the Tb and P paths parameters: /home../")
            
    
    elif isinstance(pathTb, str) and pathP is None: 
        try:        
            #Globbing the Tb files
            filenamestb = glob.glob(pathTb+'*.nc4')
                   
            #Reading Tb data
            ds_t = xr.open_mfdataset(filenamestb); 
            #Temporal resampling Tb data
            ds_t  =  ds_t.resample(time="1H").nearest(tolerance="1H") 
            
            #Reorder levels from Tb DataArray 
            ds = ds_t.transpose("lat", "lon", "time")
                       
            #Converting the UTC to local hour
            datex = ds.time.coords.to_index()
            
            #Replacing the datetimeindex based on UTC_LOCAL_HOUR
            if utc_local_sign == "minus":
                datedt = datex.to_pydatetime() - timedelta(hours=utc_local_hour)
            elif utc_local_sign == "plus":
                datedt = datex.to_pydatetime() + timedelta(hours=utc_local_hour)
            elif utc_local_sign == "local":
                datedt = datex.to_pydatetime() + timedelta(hours=utc_local_hour)
            else:
                raise TypeError("You must type a valid parameter for utc_local_sign: minus, plus or local. If you use local please enter tc_local_hour = 0")
            dates_64 = [np.datetime64(row) for row in datedt]
            ds =  ds.assign_coords({"time": dates_64})
                
            #Attaching Atributes to DataArray merged.
            ds.Tb.attrs["units"] = "K"
            ds.Tb.attrs["_FillValue"] = 'NaN'
            ds.lon.attrs['units'] = "degrees_east"
            ds.lat.attrs['units'] = "degrees_north"
            
            #Extracting dimensiones: time, lat and lon
            lon, lat = np.float32(np.meshgrid(ds.lon,ds.lat))
         
            #Establishing EPSG:4326
            ds = ds.rio.set_crs(4326)
            ds.attrs['crs'] = ds.rio.crs
            
            initial_date_lecture = str(ds.time[0].values)[:16]
            final_date_lecture = str(ds.time[-1].values)[:16]
            
            print('Complete Tb data reading ' + initial_date_lecture + " - " + final_date_lecture) 
        except:
            raise FileNotFoundError("Make sure you are complying with the Tb path parameters: /home../")
          
    else:
        raise FileNotFoundError("There must be at least a valid path for Tb data.")
    return ds


def readTRACKS(path):
    """
    function for reading tracks results. 

    Inputs:
    * path: str, path where the tracks and MCS results is located.
        
    Outputs:
    * Geopandas.GeoDataFrame with the tracks and MCS associated.
    """
    
    df = pd.read_csv(path, index_col = ["belong", "id_gdf"], parse_dates = ["time"])
    df['geometry'] = gpd.GeoSeries.from_wkt(df['geometry'])
    df['centroid_'] = gpd.GeoSeries.from_wkt(df['centroid_'])
    df = gpd.GeoDataFrame(df, geometry='geometry', crs = 4326)
    return df

def plot_folium(resume, location, path_save):
    """
    function for plotting tracks results in folium map. 

    Inputs:
    * resume: GeoDataFrame, data related with the tracks and MCS's.
    * location list (lat, lon), location for center the map_folium.
    * path_save: str, path where the .html folium map will be saved   
    
    Outputs:
    * the .html folium map will be open with the librarie "webbrowser"
    * path_saved: str, path where the .html was saved.
    """
    m = folium.Map(location=location, zoom_start=5, tiles='CartoDB positron')
    
    df = resume.reset_index()
    
    for i in df.belong.unique():
        #Sorting index by time
        tracks = df.loc[df.belong == i].reset_index()
        tracks  = tracks.set_index("time").sort_index()
        tracks = tracks.reset_index()
        for idn, r in tracks.iterrows():

            sim_geo = gpd.GeoSeries(r['geometry']).simplify(tolerance=0.001)
            geo_j = sim_geo.to_json()
            geo_j = folium.GeoJson(data=geo_j,
                                  style_function=lambda x: {'fillColor': 'orange'})
            folium.Popup(r.index).add_to(geo_j)
            try: #Tb and P methodlogy
                folium.Marker(location=[r['centroid_'].y, r['centroid_'].x], popup='id_track: {} <br> id_mcs: {} <br> hour_mcs: {} <br> time: {} <br> area[km2]: {} <br> distance_traveled[km]: {} <br> direction[°]: {} <br> intersection_percentage[%]: {} <br> mean_tb[K]: {} <br> mean_p[mm/h]: {} <br> total_distance_traveled[km]: {} <br> total_duration[h]: {} <br>'.format(r['belong'], r["id_gdf"], idn, r["time"], round(r['area_tb'],1), round(r["distance_c"],1), r["direction"],  r["intersection_percentage"], round(r["mean_tb"],1), round(r["mean_pp"],1), round(r["total_distance"],1), r["total_duration"])).add_to(m)
                extra_name = "Tb_P_"
            except: #Tb methodlogy
                folium.Marker(location=[r['centroid_'].y, r['centroid_'].x], popup='id_track: {} <br> id_mcs: {} <br> hour_mcs: {} <br> time: {} <br> area[km2]: {} <br> distance_traveled[km]: {} <br> direction[°]: {} <br> intersection_percentage[%]: {} <br> mean_tb[K]: {} <br> total_distance_traveled[km]: {} <br> total_duration[h]: {} <br>'.format(r['belong'], r["id_gdf"], idn, r["time"], round(r['area_tb'],1), round(r["distance_c"],1), r["direction"],  r["intersection_percentage"], round(r["mean_tb"],1), round(r["total_distance"],1), r["total_duration"])).add_to(m)            
                extra_name = "Tb_"
            geo_j.add_to(m)
        
    min_time = str(resume.time.min())[:-6].replace("-","_").replace(" ","_")
    max_time = str(resume.time.max())[:-6].replace("-","_").replace(" ","_")
    path_result = path_save+'map_'+extra_name+min_time+"_"+max_time+".html"
    m.save(path_result)
    try:
        webbrowser.open(path_result)
    except:
        pass
    
    return path_result
