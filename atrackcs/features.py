# -*- coding: utf-8 -*-
"""
characterizing brightness temperature, precipitation and tracks
@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
"""

#Libraries
import math
import xarray as xr  
import numpy as np
import pandas as pd
import geopandas as gpd
import tqdm
import time
import warnings
import uuid
from geopy.distance import geodesic
warnings.filterwarnings("ignore")

#___________________________________Functions______________________________________________________

class TRACKS():

    def __init__(self, geodataframe):
        self.data = geodataframe
        
    def save(self, path):
        self.save = self.data.to_csv(path)


#________________Brightness Temperature
def features_Tb(sup, ds):
    """
    Function for estimating max, mean and min brightness temperature of each MCS.
    
    Inputs
    * sup: GeodataFrame, with the Tb polygons generated in the process detect_mcs().
    * ds: Dataset, with the Tb and P data resulting from the process funcs.readNC().
    
    Outputs:
    * GeoDataFrame, which identifies the MCS's polygons with some Tb attributtes
    * Attributes:
        * time: datetime64, hour (UTC-5) for this case.
        * Tb: float, polygon index after being filtered by the established parameterization.
        * geometry: geometry, MCS polygon (convex hull). 
          The crs in this case corresponds to utm_local_zone.
        * area_Tb: float, area polygon [km**2]
        * centroid_: geometry, geometric centroid polygon (convex hull). 
          The crs in this case corresponds to utm_local_zone.
        * mean_tb: float, Tb average of the pixels composing the polygon [K].
        * min_tb: float, Tb min value of the pixels composing the polygon [K].
        * max_tb: float, Tb max value of the pixels composing the polygon [K]. 
    """
    #Preparing the DataArray
    data_magnitud = ds.copy()
    data_magnitud = data_magnitud.to_dataframe()
    data_magnitud_ = data_magnitud.to_xarray()["Tb"]
    
    #Establishing EPSG:4326 - WGS geodetic coordinate system 
    data_magnitud_ = data_magnitud_.rio.write_crs(4326)  
    
    #Renaming dimensions 
    data_magnitud_= data_magnitud_.rename({"lat":"y", "lon":"x"})
    gdf_tb = sup['geometry']

    #Creating new columns for the mean, max and min Tb
    sup['mean_tb'] = None; sup['min_tb']= None; sup['max_tb']= None

    print ("Estimating MCS's brightness temperature attributes: ")    
    #Estimating Tb features
    for index_t,_dates, progress in zip(gdf_tb.index,sup.time, tqdm.tqdm(range(len(sup.time)-1))):
        _polygon = gdf_tb.geometry.loc[index_t]
        coordinates = np.dstack((_polygon.exterior.coords.xy[0], _polygon.exterior.xy[1]))
        geometries = [{'type': 'Polygon', 'coordinates': [coordinates[0]]}]    
        blob_clipped = data_magnitud_.loc[:,:,str(_dates)].rio.clip(geometries, gdf_tb.crs, drop=False, invert=False)
        
        temp_minima = round(float(blob_clipped.min(skipna=True).values),4)
        temp_promedio = round(float(blob_clipped.mean(skipna=True).values),4)
        temp_max = round(float(blob_clipped.max(skipna=True).values),4)
        
        sup.loc[index_t,"mean_tb"] = temp_promedio
        sup.loc[index_t,"min_tb"] = temp_minima
        sup.loc[index_t,"max_tb"] = temp_max
        
        sup['mean_tb'] = sup["mean_tb"].astype(float)
        sup['min_tb'] = sup["min_tb"].astype(float)
        sup['max_tb'] = sup["max_tb"].astype(float)

        #Progress bar
        time.sleep(0.01)
        
    return sup

#________________Precipitation

def confirm_P_data(func):
    """
    Decorator for checking and blocking the use of P data when the detect_scheme is only by Tb
    """
    def valid_P(sup, ds, min_precipitation, area_P, drop_empty_precipitation):
        try:
            ds["P"]
        except:
            raise FileNotFoundError("Is not possible to estime P features because there is no P data")
        return func(sup, ds, min_precipitation, area_P, drop_empty_precipitation)
    return valid_P
        
@confirm_P_data
def features_P(sup, ds, min_precipitation = 2, area_P = 500, drop_empty_precipitation = False):
    """
    Function to estimate the existing precipitation (mean P and max P) inside the MCS. 
    These characteristics will only be calculated for MCS that satisfy the conditions set 
    in "min_precipitation" and "area_P".
    
    Inputs:
    * sup: GeoDataframe, result data generated in the process features_Tb().
    * ds: Dataset, with the Tb and P data resulting from the process funcs.readNC().
    * min_precipitation: minimum precipitation limit in each pixel (ex. 2 mm/h)
    * P_area: minimun area P polygon in Tb polygon (ex. 500 km2)
    * drop_empty_precipitation: boolean, if True remove MCS that do not contain precipitation 
      values within the established criteria. This is useful if you want to establish tracks
      only with MCS that satisfy the precipitation criteria. 

    Outputs:
    * GeoDataFrame, which identifies the MCS's polygons with some Tb and P attributtes
    * Attributes: this GeoDataFrame contains exactly the information referenced in features_Tb()
      except for the dropping of some features associated to the Tb and the addition of 
      some characteristics associated to the MCS precipitation.
        * mean_pp: float, P average of the pixels composing the polygon [mm/h].
        * max_pp: float, P max value of the pixels composing the polygon [mm/h].
    """
    
    if ds["P"].time.size != 0:
    
        #Preparing the DataArray
        data_magnitud = ds.copy()
        data_magnitud = data_magnitud.to_dataframe()
        data_magnitud_ = data_magnitud.to_xarray()["P"]
        
        #Establishing EPSG:4326 - WGS geodetic coordinate system 
        data_magnitud_ = data_magnitud_.rio.write_crs(4326)  
    
        #Renaming dimensions 
        data_magnitud_= data_magnitud_.rename({"lat":"y", "lon":"x"}) 
        gdf_tb = sup['geometry']
        
        #Creating new columns for the mean and max P
        sup['mean_pp'] = None; sup['max_pp']= None
    
        #Estimating P features    
        print ("Estimating MCS's precipitation attributes: ")    
        
        #Translating "area_P" in pixels (~10 km x ~10km)
        valid_pixels_precipitation = int((area_P/100.))
        
        
        for index_t,_dates, progress in zip(gdf_tb.index,sup.time,tqdm.tqdm(range(len(sup.time)-1))):
            _polygon = gdf_tb.geometry.loc[index_t]
            coordinates = np.dstack((_polygon.exterior.coords.xy[0], _polygon.exterior.xy[1]))
            geometries = [{'type': 'Polygon', 'coordinates': [coordinates[0]]}]
    
            #Applying criteria 2 mm/hr
            blob_clipped = data_magnitud_.loc[:,:,str(_dates)].rio.clip(geometries, gdf_tb.crs, drop=False, invert=False)
            blob_clipped = blob_clipped.where(blob_clipped>=min_precipitation) 
            
            #Checking if the DataArray has at least the valid pixels with precipitation.
            if blob_clipped.notnull().sum() >= valid_pixels_precipitation:
                max_magnitud = round(float(blob_clipped.max(skipna=True).values),4)
                mean_magnitud = round(float(blob_clipped.mean(skipna=True).values),4)
            else:
                max_magnitud = np.nan
                mean_magnitud = np.nan
            
            sup.loc[index_t,"mean_pp"] = mean_magnitud
            sup.loc[index_t,"max_pp"] = max_magnitud
            
            sup['mean_pp'] = sup["mean_pp"].astype(float)
            sup['max_pp'] = sup["max_pp"].astype(float)
            
            #Progress bar
            time.sleep(0.01)       
            
            #Dropping polygons that do not meet the precipitation parameters
            #in terms of area and amount.
        if drop_empty_precipitation == True:
            sup = sup[sup['mean_pp'].notna()].reset_index(drop=True)
        else: 
            pass
    else:
        raise TypeError("There is no P data")
    return sup

#________________Trajectories
    
def distance_centroids(row):
    """
    Function to estimate the distance between two georreferenced points.
    This function execute inside the function distance_direction_Tracks().
    Inputs:
    * row: GeoDataFrame, containing the shift of the MCS at time t (centroid_)
      and the MCS at time t+1 (centroids_distance)
    
    Outputs:
    * distance: float, between the geometric centroids.
    """
    try:
        d = geodesic((row.centroid_.y,row.centroid_.x),(row.centroids_distance.y,row.centroids_distance.x)).km
    except:
        d = np.nan
    return d

def direction_points(row, mode = "u"):
    """
    Function to estimate the direction between two georreferenced points.
    
    Inputs:
    * row: GeoDataFrame, containing the shift of the MCS at time t (centroid_)
      and the MCS at time t+1 (centroid_2)
    * mode: str(u, v, deg), select the result of based on the output

        
    Outputs: The output is in function of the mode selected
    * deg: direction in degrees between the two points
    * u: normalized vector component
    * v: normalized vector component
    """
    
    try:
        #distance vector between the centroids of two points - in x: lon and y: lat
        distance = [row.centroid_.x - row.centroid_2.x, row.centroid_.y - row.centroid_2.y]
        
        #Normal Vector
        norm = math.sqrt(distance[0] ** 2 + distance[1] ** 2) #Based on Pitagoras Teorem
        direction_normalized = [distance[0] / norm, distance[1] / norm] #Divided by normalized vector
        u = direction_normalized[0]; v = direction_normalized[1]
        
        #Conversion of vector components into direction degrees     
        direccion_deg = uv_to_degrees(u,v)
        direccion_deg = direccion_deg.round(1)
    except:
        direccion_deg = np.nan
        u = np.nan; v = np.nan
        
    if mode == "u":
        var = u
    elif mode == "v":
        var = v
    elif mode == "deg":
        var = direccion_deg
    return var

def uv_to_degrees(U,V):
    """
    Calculates the direction in degrees from the u and v components.
    Takes into account the direction coordinates is different than the 
    trig unit circle coordinate. If the direction is 360 then returns zero
    (by %360). This function execute inside the function distance_direction_Tracks().
    
    Inputs:
    * U: float, west/east direction (from the west is positive, from the east is negative)
    * V: float, south/noth direction (from the south is positive, from the north is negative)
      
    Outputs:
    * direction in degrees: float  
    """
    WDIR= (270-np.rad2deg(np.arctan2(V,U)))%360
    return WDIR


def distance_direction_Tracks(sup):
    """
    Function to estimate the distance and direction of the MCS's that compose 
    the tracks. This function execute inside the function feature_Tracks().
     
    Inputs:
    * sup: GeoDataframe, result data generated in the process features_Tb() or features_P().

    Outputs:
    * Geodataframe, with the distances and directions of the MCS's estimated.
    """
    
    #Making new columns
    sup["distance_c"] = None; contador = 0
    sup["direction"] = None; sup["u"] = None; sup["v"] = None
    sup["estado"] = None;
    sup["estado_porcentaje"] = None

    tracks = sup.belong.unique(); len_track = len(tracks)

    #Estimating distance and direction between two points   
    
    print("Estimating distance and direction between geometrics centroids: ")
    for track, progress in zip(tracks, tqdm.tqdm(range(len_track-1))):
        
        #____________________distance___________________________
        track_df = sup.loc[sup["belong"] == track, sup.columns]        
        #Generating centroid to compare (time i + 1) in another column               
        track_df["centroids_distance"] = track_df["centroid_"].shift()
        
        #Calculating the distance between centroids (time i and time i + 1)        
        track_df["distance_c"] = track_df.apply(lambda row: distance_centroids(row), axis = 1)
        
        #Obtaining the ids of the centroids (spots) that were obtained.
        index_track_distancia = track_df.index    
        
        #Replacing the obtained distance values in the original Dataframe
        sup.loc[index_track_distancia, "distance_c"] = track_df["distance_c"]
        
        #____________________direction___________________________

        #Generating centroid to compare (time i - 1) in another column               
        track_df["centroid_2"] = track_df["centroid_"].shift(-1)
        
        #Calculating the direction between centroids ((time i and time i - 1))
        track_df["direction_fake"] = track_df.apply(lambda row: direction_points(row, mode = "deg"), axis = 1)
        track_df["u_fake"] = track_df.apply(lambda row: direction_points(row, mode = "u"), axis = 1)
        track_df["v_fake"] = track_df.apply(lambda row: direction_points(row, mode = "v"), axis = 1)

        #Aligning the Dataframe so that the calculated direction corresponds to time i + 1
        track_df["direction"] = track_df["direction_fake"].shift(1)
        track_df["u"] = track_df["u_fake"].shift(1)
        track_df["v"] = track_df["v_fake"].shift(1)

        del track_df["direction_fake"]; del track_df["v_fake"]; del track_df["u_fake"]
               
        #Generating qualitative and quantitative information on the growth or decrease of the spot based on last state
        track_df["estado"] = np.where(track_df["area_tb"].diff() > 0, 'CRECIENDO', 'DECRECIENDO')
        track_df["diferencia_area"] = track_df["area_tb"].diff()
        track_df["area_tb_fake"] = track_df["area_tb"].shift(1)
        track_df["estado_porcentaje"] = track_df.apply(lambda track_df: 
            track_df["diferencia_area"]*100/track_df["area_tb_fake"], axis = 1)
        track_df.loc[(track_df["estado_porcentaje" ] >= - 10.) & (track_df["estado_porcentaje" ] <= 10.), "estado"] = "ESTABLE"
        
        first_index_main = track_df.iloc[0].name
        
        #Assigning nan to the first record of the track
        track_df.loc[first_index_main, "estado"] = np.nan

        del track_df["area_tb_fake"]; del track_df["diferencia_area"] 
        
        #Replacing the obtained values of direction (degrees) in the original dataframe
        index_track_direccion = track_df.index    
        sup.loc[index_track_direccion, "direction"] = track_df["direction"]   
        sup.loc[index_track_direccion, "u"] = track_df["u"]   
        sup.loc[index_track_direccion, "v"] = track_df["v"]   
        
        #Replacing status and percentage status values
        sup.loc[index_track_direccion, "estado"] = track_df["estado"]   
        sup.loc[index_track_direccion, "estado_porcentaje"] = track_df["estado_porcentaje"]   
                               
        contador +=1
        #porcentaje = round((contador*100./len_track),2)
        
        #Progress bar
        time.sleep(0.01)     

    sup.distance_c = sup.distance_c.astype(float)
    sup.direction = sup.direction.astype(float)
    sup.u = sup.u.astype(float)
    sup.v = sup.v.astype(float)

    #print ("Porcentaje de distancias-direcciones-estados calculadas " + str(porcentaje) + " %")
    #print ("distances-directions-states of tracks completed")
    return sup
    
def features_Tracks(sup, initial_time_hour = 0, encrypt_index = True,
                    path_save = None):
    """
    Function for calculating attributes associated with the tracks
    
    Inputs:
    * sup: GeoDataFrame, result data generated in the process track_mcs() 
    * initial_time_hour: int, default is 0 but could change based for filtering
      tracks based on total duration.
    * encrypt_index: boolean, default is True. The index generated is encrypted 
      and generated using the uuid library for each MCS and each track.
    * path_save: str, path where the .csv results will be saved.   
    
    Outputs:
    * GeoDataFrame, which identifies the tracks and the MCS's polygons with some Tb, P 
      and tracks attributtes
    * Attributes: this GeoDataFrame contains exactly the information referenced in features_Tb()
      and features_P() except for the addition of some attributes associated to the tracks.
        * belong: str, encrypted index generated for each track.
        * id_gdf: str, encrypted index generated for each MCS.
        * geometry: geometry, polygon. The crs in this case corresponds to WGS84 - EPSG:4326.
        * centroid_: geometry, geometric centroid polygon. The crs in this case corresponds
          to WGS84 - EPSG:4326.
        * intersection_percentage: float, percentage overlap between MCS [%].
        * distance_c: float, distance between the overlapping geometric 
          centroids of the MCS [km].
        * direction: float, direction between the overlapping geometric 
          centroids of the MCS [°].
        * total_duration: float, total duration of the event or the track. 
          This value is associated to each MCS of the corresponding track [h].
        * total_distance: float, total distance of the event or the track. 
          This value is associated to each MCS of the corresponding track [km].
        * mean_velocity: float, velocity average of the event or the track. 
          This value is associated to each MCS of the corresponding track [km/h].
    """      
		
    #Estimating distance and direction between geometrics centroids
    sup = distance_direction_Tracks(sup)

    #Preparing Dataframe   
    NUEVODF = sup.copy()
    NUEVODF = NUEVODF.set_index(['belong', 'id_gdf']).sort_index()

    #Replacing old id for spots and tracks based on alphanumeric code 16 and 20
    #characteres respectively    
    new_belong = []
    new_track_id = []
    

    #Encriptying index track and index spot     
    if encrypt_index == True:
        for track_id in NUEVODF.index.levels[1]:
            new_track_id.append(str(uuid.uuid4())[-22:])
        dic_replace_track_id = dict(zip(NUEVODF.index.levels[1].values, new_track_id))
        
        for belong_id in NUEVODF.index.levels[0]:
            new_belong.append(str(uuid.uuid4())[:13])
        dic_replace_belong = dict(zip(NUEVODF.index.levels[0].values, new_belong))

        reg_sup_res =  NUEVODF.reset_index()
        
        reg_sup_res.belong = reg_sup_res.belong.replace(dic_replace_belong)
        reg_sup_res.id_gdf = reg_sup_res.id_gdf.replace(dic_replace_track_id)
    else:
        reg_sup_res =  NUEVODF.reset_index()
        pass
    

       
    reg_sup_res = reg_sup_res.drop(labels=['Tb', 'min_tb', 'max_tb','u','v', 'estado','estado_porcentaje'],axis=1)
       
    reg_sup = reg_sup_res.set_index(["belong" , "id_gdf"]).sort_index()
    
    #Attaching to the dataframe total duration and total distance for each spot. Each spot 
    #has the register of the features of his own track
    reg_sup["total_duration"] = None
    reg_sup["total_distance"] = None
    
    count_df = reg_sup_res.groupby(by = [reg_sup_res.belong]).count()
    sum_df = reg_sup_res.groupby(by = [reg_sup_res.belong]).sum()
    
    for _b in count_df.index: 
    
        count_value = count_df.loc[_b,"id_gdf"]
        sum_value = sum_df.loc[_b,"distance_c"]

        reg_sup.loc[_b, "total_duration"] = count_value
        reg_sup.loc[_b, "total_distance"] = sum_value

    #Estimating track mean velocity      
    reg_sup['mean_velocity'] = reg_sup["total_distance"]/reg_sup["total_duration"]
    
    #Leaving only the tracks based on initial_time_hour
    reg_sup = reg_sup.reset_index()
    reg_sup = reg_sup[reg_sup['total_duration'] > initial_time_hour].sort_index()
    reg_sup = reg_sup.set_index(["belong", "id_gdf"])
    reg_sup = reg_sup.sort_index()
    
    #Saving as .csv the results
    #reg_sup.to_csv(pathResultados + "resume_"+str(reg_sup.time.min())[:-7]+"_"+str(reg_sup.time.max())[:-7]+"_"+str(extra_name)+".csv")
    
    #reg_sup.columns = ['time', 'geometry', 'area_tb', 'centroid', 'mean_tb', 'mean_p',
    #   'max_p', 'intersection_percentage', 'distance_c', 'direction',
    #   'total_duration', 'total_distance', 'mean_velocity']
    
    
    reg_sup_ = TRACKS(reg_sup)
    
    if path_save is None:
        pass
    else:
        try:
            #Saving as .csv the results
            
            min_time = str(reg_sup.time.min())[:-6].replace("-","_").replace(" ","_")
            max_time = str(reg_sup.time.max())[:-6].replace("-","_").replace(" ","_")
            
            if 'max_pp' in reg_sup.columns:
                extra_name = "Tb_P_"
            else:
                extra_name = "Tb_"
            
            path_result = path_save+'resume_'+extra_name+min_time+"_"+max_time+".csv"
            #print(path_result)
            reg_sup.to_csv(path_result)
        except:
            raise FileNotFoundError("No such file or directory: " + str(path_save))
       
    return reg_sup



