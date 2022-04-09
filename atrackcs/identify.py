# -*- coding: utf-8 -*-
"""
identifying mesoscale convective systems functions
@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
"""
#Libraries
import xarray as xr  
import numpy as np
import pandas as pd
from scipy import ndimage
import geopandas as gpd
from datetime import timedelta
import rioxarray
import rasterio
from geopy.distance import geodesic
import numbers as checknumbers
from shapely.geometry import MultiPolygon, Polygon, shape, Point, MultiPoint, mapping 
from shapely.wkt import loads
import warnings
warnings.filterwarnings("ignore")


#___________________________________Functions______________________________________________________


def polygon_identify(data_id, variable):
    """
    Function to draw polygons associated to each P or Tb spot; the polygon delimitation 
    is from the Convex Hull. This function execute inside the function identify_msc2. 
    
    Inputs:
    data_id = DataArray associated with each variable: P or Tb
    variable = string name variable: "Tb" or "P"
    
    Outputs:
    Geodataframe with the polygons associated to the P and tb spots.
    """
    #structure that determines which elements are neighbors to each pixel.
    structure = ndimage.generate_binary_structure(2,2)

    for step in data_id.time:
        blobs, number_of_blobs = ndimage.label(data_id.sel(time = step).values, structure=structure)
        distance = ndimage.distance_transform_edt(blobs)
    
        #only the outline of the polygons (spots) based on binary structure.
        distance[distance != 1] = 0
        blobs[distance == 0]= 0
        blobs = blobs.astype('float')
        blobs[blobs  == 0] = 'nan'
        
        #replacing the values generated in the original Datarray
        data_id.loc[dict(time = step.values)] = blobs
            
    #___________________Polygon generation of the spots______________________________
        
    #Datarray to Dataframe
    df = data_id.to_dataframe().reset_index()
    df.dropna(inplace = True)
    #Dataframe to Geodataframe
    gdf = gpd.GeoDataFrame(
            df[["time", variable]], geometry=gpd.points_from_xy(df.lon,df.lat))
     
    #Extracting the coordinates of the georeferenced points in a separate column
    gdf['geo'] = gdf['geometry'].apply(lambda x: x.coords[0])
    
    #Converting the points-groups to polygons based on the Convex Hull
    df = gdf.groupby(by=["time" , variable])['geo'].apply(lambda x: MultiPoint(sorted(x.tolist()))).to_frame().reset_index()
    df["geo"] = df["geo"].apply(lambda x: x.convex_hull.wkt).to_frame()
    df["geo"] = df["geo"].apply(lambda x: loads(x)).to_frame()
    
    #Making the Geodataframe final 
    gdf = gpd.GeoDataFrame(df, geometry=df.geo, crs=4326); del df; del gdf["geo"]
    #Dropping lines and solitary points from the Geodataframe
    gdf = gdf.loc[gdf['geometry'].geom_type=='Polygon']
         
    return gdf      

def thresholds(func):
    def valid_thresholds(data, variables, Tb, area_Tb, utm_local_zone):
        if Tb not in range(200,241):
            raise ValueError("You must enter a value of Tb between 200 K - 240 K")
        if variables not in ["Both", "Tb"]:
            raise TypeError("You must type a valid parameter for variables: Tb or Both")
        if not area_Tb >999:
            raise ValueError("You must enter an area_Tb greater or equal than 1000 km2")
        if not isinstance(utm_local_zone, int):
            raise ValueError("You must enter a valid UTM local zone")
        return func(data, variables, Tb, area_Tb, utm_local_zone)
    return valid_thresholds

def confirm_P_data(func):
    def valid_variables(data, variables, Tb, area_Tb, utm_local_zone):
        if variables == "Both":  
            try:
                data["P"]
            except:
                raise TypeError("There is not P data, please change the parameter variables for Tb")
        return func(data, variables, Tb, area_Tb, utm_local_zone)
    return valid_variables
           

@confirm_P_data
@thresholds
def identify_mcs(data, variables = "Both", Tb = 225, area_Tb = 2000, utm_local_zone = None):
    """
    Function for identifying the spots, in a time step, based on the methodology
    of brightness temperature and precipitation association.
 
    Inputs
    data = DataArray of the variables
    variables: Methodology of association "Tb" or "Both":Tb and P. 
    If "Tb" the spots area are delimited by brightness temperature cold cloud top
    Elif "Both" the spots area are delimited by brightness temperature cold cloud top and
    precipitation threshold in each spot


    #default parameters based on literature
    Tb (Brightness Temperature): spots based on limited maximun threshold cold cloud top (ex. <225 K)[Feng et al.,(2021); Li et al.,(2020)]
    area_Tb: spots with a minimun largest area polygon (ex. > 2000 km2)[Lui et al., (2019); Vizy & Cook,(2018)] 
    P (Precipitation): By default is the spots greater than 2 mm/h with length[Feng et al.,(2021)].         

    Outputs:
    Geodataframe with the polygons associated to the methodology choosed
    """ 
    
    global counter, storms_counter
    
    #Converting DataArray to Dataframe
    dataf = data.to_dataframe()

    #Methodology Tb
    if variables == "Tb":
        #Masking the dataframe except for the outline of the spots based
        #on limited threshold cold cloud top 
        dataf.loc[dataf["Tb"]<=Tb] = 1
        dataf.loc[dataf["Tb"]!=1] = 0
        
        #Converting Dataframe to DataArray
        datax = dataf.to_xarray()["Tb"]; del dataf

    #Methodology Tb and P        
    elif variables =='Both':  
        #Masking the dataframe except for the outline of the spots based
        #on limited threshold cold cloud top         
        dataf.loc[dataf["Tb"]<=Tb] = 1
        dataf.loc[dataf["Tb"]!=1] = 0
       
        #Converting Dataframe to DataArray
        datax = dataf.to_xarray()['Tb']
        
        #Establishing EPSG:4326
        datax = datax.rio.write_crs(4326)

        #Masking the P dataframe except for the outline of the spots based on 2mm/hr       
        dataf.loc[dataf["P"]>=2] = 1
        dataf.loc[dataf["P"]!=1] = 0
        
        #Converting Dataframe to DataArray
        datay = dataf.to_xarray()["P"]
        
        #Establishing EPSG:4326 - WGS geodetic coordinate system 
        datay = datay.rio.write_crs(4326)
        del dataf 
        
    #Identifying polygons of Tb and P of the original dataArray  
    gdf_tb = polygon_identify(datax, variable = 'Tb')

    #Establishing EPSG:32718 - UTM zone 18S plane coordinate system 
    gdf = gdf_tb.to_crs(utm_local_zone) 
    
    #Calculating polygon area   
    gdf["area_tb"] = gdf['geometry'].area/ 10**6 #meters to kilometers
    gdf["area_tb"] = gdf["area_tb"].round(1)
    
    #Dropping polygons with an area less than area_Tb (default: 2000 km2)
    gdf = gdf.loc[gdf['area_tb'] >= area_Tb]    
    
    #Calculating centroids
    gdf["centroides"] = gdf.geometry.centroid
    gdf.reset_index(inplace = True, drop = True)
        
    print("Spots identification completed")

    return gdf

    
