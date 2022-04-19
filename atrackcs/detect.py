# -*- coding: utf-8 -*-
"""
identifying mesoscale convective systems functions
@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
          
Bibliography
* Feng, Z., Leung, L. R., Liu, N., Wang, J., Houze, R. A., Li, J., Hardin, J. C., Chen, D., & Guo, J. (2021). A Global High‐resolution Mesoscale Convective System Database using Satellite‐derived Cloud Tops, Surface Precipitation, and Tracking. Journal of Geophysical Research: Atmospheres. https://doi.org/10.1029/2020jd034202
* Li, J., Feng, Z., Qian, Y., & Leung, L. R. (2020). A high-resolution unified observational data product of mesoscale convective systems and isolated deep convection in the United States for 2004–2017. Earth System Science Data Discussions, October, 1–48. https://doi.org/10.5194/essd-2020-151
* Liu, W., Cook, K. H., & Vizy, E. K. (2019). The role of mesoscale convective systems in the diurnal cycle of rainfall and its seasonality over sub-Saharan Northern Africa. Climate Dynamics, 52(1–2), 729–745. https://doi.org/10.1007/s00382-018-4162-y
* Vizy, E. K., & Cook, K. H. (2018). Mesoscale convective systems and nocturnal rainfall over the West African Sahel: role of the Inter-tropical front. Climate Dynamics, 50(1–2), 587–614. https://doi.org/10.1007/s00382-017-3628-7          
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

class MCS():

    def __init__(self, geodataframe):
        self.data = geodataframe
        
    def save(self, path):
        self.save = self.data.to_csv(path)
    

def polygon_identify(data_id, variable):
    """
    Function to estimate polygons (convexhull) associated to the
    interest variable (Tb or P). This function execute inside the 
    function detect_mcs().
    
    Inputs:
    * data_id: DataFrame, with data associated to the detect scheme. 
    * variable: str (Tb or P), variable's name for estimating the polygons.
    
    Outputs:
    * Geodataframe, with the estimated polygons.
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
            
    #_________________Spots's polygon generation______________________________
        
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
    """
    Decorator for checking the thresholds and the type of enter data in detect_mcs()
    """
    def valid_thresholds(data, detect_scheme, Tb, area_Tb, utm_local_zone, path_save):
        if Tb not in range(200,241):
            raise ValueError("You must enter a value of Tb between 200 K - 240 K")
        if detect_scheme not in ["Both", "Tb"]:
            raise TypeError("You must type a valid parameter for detect_scheme: Tb or Both")
        if not area_Tb >999:
            raise ValueError("You must enter an area_Tb greater or equal than 1000 km2")
        if not isinstance(utm_local_zone, int):
            raise ValueError("You must enter a valid UTM local zone")
        if path_save == None:
            pass
        return func(data, detect_scheme, Tb, area_Tb, utm_local_zone, path_save)
    return valid_thresholds

def confirm_P_data(func):
    """
    Decorator for checking and blocking the use of P data when the detect_scheme is Tb
    """
    def valid_variables(data, detect_scheme, Tb, area_Tb, utm_local_zone, path_save):
        if detect_scheme == "Both":  
            try:
                data["P"]
            except:
                raise TypeError("There is not P data, please change the parameter detect_scheme for Tb")
        return func(data, detect_scheme, Tb, area_Tb, utm_local_zone, path_save)
    return valid_variables
           
@confirm_P_data
@thresholds
def detect_mcs(data, detect_scheme = "Both", Tb = 225, area_Tb = 2000, utm_local_zone = None, path_save = None):
    """
    Function for detecting MCS, in a time step, based on the detect scheme 
    and Tb criteria selection.

    Inputs
    * data = Dataset, with the data resulting from the process funcs.readNC().
    * detect_scheme: str(Tb, Both), association scheme. If "Tb" the MCS area are 
    delimited by Tb cold cloud top. Elif "Both" the MCS area are delimited
    by Tb cold cloud top and precipitation criteria selection.
    * utm_local_zone: int, is needed for converting the WGS geodetic coordinate system 
    * path_save: str, path where the .csv will be saved.
    
    #default parameters based on literature
    * Tb, int: MCS based on limited maximun threshold cold cloud top 
     (ex. <225 K)[Feng et al.,(2021); Li et al.,(2020)]
    * area_Tb: MCS with a minimun largest area polygon 
     (ex. > 2000 km2)[Lui et al., (2019); Vizy & Cook,(2018)] 

    Outputs:
    * GeoDataFrame, with the MCS records.
    """ 
    
    global counter, storms_counter
    
    #Converting DataArray to Dataframe
    dataf = data.to_dataframe()

    #Methodology Tb
    if detect_scheme == "Tb":
        #Masking the dataframe except for the outline of the spots based
        #on limited threshold cold cloud top 
        dataf.loc[dataf["Tb"]<=Tb] = 1
        dataf.loc[dataf["Tb"]!=1] = 0
        
        #Converting Dataframe to DataArray
        datax = dataf.to_xarray()["Tb"]; del dataf

    #Methodology Tb and P        
    elif detect_scheme =='Both':  
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
    gdf["centroid_"] = gdf.geometry.centroid
    gdf.reset_index(inplace = True, drop = True)
        
    print("MCS detection completed")
    
    gdf = MCS(gdf)

    
    if path_save is None:
        pass
    else:
        try:
            gdf.save(path_save)
        except:
            raise FileNotFoundError("No such file or directory: " + str(path_save))

    return gdf.data

    
