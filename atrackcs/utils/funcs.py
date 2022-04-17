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
warnings.filterwarnings("ignore")


#___________________________________Functions______________________________________________________

def readNC(pathTb = None, pathP = None, utc_local_hour = None):
    """
    Function for reading and resampling the Tb and P DataArrays. 
    The spatial resampling is 0.1° - lineal interpolation.
    The temporal resampling is 1 h - nearest original coordinate to up-sampled frequency coordinates.
    
    Inputs:
    pathTb: Path where the Tb raw data are located.
    pathP: Path where the P raw data are located.
    The path must have the next structure:
    linux: r"/home....../"
    windows: r"C:/....../"
    
    UTC_LOCAL_HOUR: Is needed for converting the raw data hour (UTC) to a local hour (interest region)
    
    Outputs:
    Dataset with the  P and Tb data.
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
            datedt = datex.to_pydatetime() - timedelta(hours=utc_local_hour)
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
            datedt = datex.to_pydatetime() - timedelta(hours=utc_local_hour)
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

def plot_folium(resume, location):
    m = folium.Map(location=[5, -73.94], zoom_start=5, tiles='CartoDB positron')
    
    df = resume.reset_index()
    
    for i in df.belong.unique():
        xx = sup.loc[df.belong == i].reset_index()
        for idn, r in xx.iterrows():
            #without simplifying the representation of each borough, the map might not be displayed
            #sim_geo = gpd.GeoSeries(r['geometry'])
            sim_geo = gpd.GeoSeries(r['geometry']).simplify(tolerance=0.001)
            geo_j = sim_geo.to_json()
            geo_j = folium.GeoJson(data=geo_j,
                                  style_function=lambda x: {'fillColor': 'orange'})
            folium.Popup(r.index).add_to(geo_j)
            folium.Marker(location=[r['centroid'].y, r['centroid'].x], popup='track: {} <br> id_track: {} <br> id_gdf: {} <br> area[km2]: {} dirección: {} <br> distancia[km]: {}  <br> tiempo {} <br>'.format(r['belong'], idn ,r["id_gdf"] , r['area_tb'], r["direction"], round(r["total_distance"],1) ,r['time'])).add_to(m)
            geo_j.add_to(m)
        

            m.save(pathResultados+'map_tb225_dur6hr_MAM2019.html')