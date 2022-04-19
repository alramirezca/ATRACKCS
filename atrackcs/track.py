# -*- coding: utf-8 -*-
"""
tracking mesoscale convective systems functions
@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
"""

import xarray as xr  
import numpy as np
import pandas as pd
import geopandas as gpd
from datetime import timedelta
import rioxarray
import rasterio
import numbers as checknumbers
import tqdm
import time
import warnings
from shapely.geometry import Polygon
warnings.filterwarnings("ignore")


#___________________________________Functions______________________________________________________

class TRACKS():

    def __init__(self, geodataframe):
        self.data = geodataframe
        
    def save(self, path):
        self.save = self.data.to_csv(path)

def min_dist(point, gpd2):
    """
    Function to find the nearest polygon based on a georeferenced point.
    
    Inputs:
    * point: geometry, georeferenced point.
    * gpd2: GeoDataFrame, data polygons in a specific time and region
    
    Outputs:
    * GeoSeries, the nearest polygon to the georeferenced point.
    """
    gpd = gpd2.copy()
    gpd['Dist'] = gpd.apply(lambda row:  point.distance(row.geometry),axis=1)
    geoseries = gpd.iloc[gpd['Dist'].argmin()]
    return geoseries

def track_mcs(sup, threshold_overlapping_percentage = None, utm_local_zone = None, path_save = None):
    global msc_counter
    """
    Function for tracking convective systems according to the threshold_overlapping_percentage.
    This functions works based on identified convective systems in a period of time.   

    Inputs:
    * sup: GeoDataFrame, result data generated in the process features_Tb() 
    and features_P (when is using this detect scheme).
    * threshold_overlapping_percentage: float, by default there is no percentage overlap limit 
    between polygons, however, this can be set. 
    (ex. 10: this means that the percentage of overlap that is not greater than or equal to 10%, 
    in reference to the largest area polygon, is not taken into account in the trajectory).
    * utm_local_zone: int, is needed for converting the WGS geodetic coordinate system 
    to plane coordinate system. This is a constant that must be asociated with the interest 
    region. (ex. 32718 : is the UTM zone 18S plane coordinate system and is used for tracking 
    MCS in Colombia, Brasil, Peru, Ecuador, Chile.).
    * path_save: str, path where the .csv will be saved.
    
    Outputs:
    * GeoDataFrame, which identifies the tracks and his associated MCS's.              
    """
    msc_counter = 1
    pd.options.mode.chained_assignment = None
    range_sites = sup.index
    #Making column belong
    sup.loc[sup.index, "belong"] = 0
    #Making column percentage intersection
    sup.loc[sup.index, "intersection_percentage"] = np.nan
        
    #Generating tracks
    
    print ("Estimating trajectories: " )
    for isites, progress in zip(range_sites, tqdm.tqdm(range(len(range_sites)-1))):
        # Checking when the record (spot) has not yet been processed or attaching to a track
        if sup.at[isites, "belong"] == 0:
            try:
                #Assigning the counter number to identify the track(range 1) 
                sup.at[isites, "belong"] = msc_counter #Index, col  
                
                #Generating start and end time date (time t and time t+1)
                date_initial = sup.loc[isites].time
                date_final = date_initial + timedelta(hours=1)
                
                #Reference spot
                poli_ti = sup.loc[isites].to_frame().T
                #Establishing EPSG:32718 - UTM zone 18S plane coordinate system 
                poli_ti = gpd.GeoDataFrame(poli_ti, geometry='geometry', crs=utm_local_zone)
                
                #Spots one hour later
                poli_ti1 = sup.loc[sup.time == date_final]
                
                #Intersection of reference spot with spots one hour later
                intersect = gpd.overlay(poli_ti, poli_ti1, how='intersection')
                                                      
                #Nearest intersection point (centroid)
                
                try: #When there is only one intersection spot
                    
                    #Estimating the distance between centroids
                    dist = min_dist(intersect.centroid, poli_ti1)
                                
                    #Extracting the spot with the maxima area
                    intersect_area_max = max(intersect.area_tb_1.values, intersect.area_tb_2.values)
    
                    #Estimating percentage of intersection
                    intersect_percentage = ((intersect.area/1000000)*100/intersect_area_max)[0].round(1)                
                    
                    #Condition if threshold overlapping percentage is activated
                    #If the spot fails with threshold overlapping percentage is dropped
                    if isinstance(threshold_overlapping_percentage, checknumbers.Real):
                        if intersect_percentage>= threshold_overlapping_percentage:
                            sup.at[isites, "intersection_percentage"] = intersect_percentage #Index, col 
                            #print ("There was a " + str(intersect_percentage) + "% intersection") 

                        else:
                            sup.drop(isites, inplace = True)
                            #print ("The convective system " + str(isites) + " was dropped to the data due intersection is below that the threshold: "+ str(intersect_percentage)) 

                    #Condition if threshold overlapping percentage is not activated
                    else:
                        sup.at[isites, "intersection_percentage"] = intersect_percentage #Index, col 
                        #print ("There was a " + str(intersect_percentage) + "% intersection") 

                except: #When there are more than one intersection spot is selected the bigger spot
                    
                    #dimensioning the Geodataframe columns only for the intersected spots
                    intersect = intersect[intersect.columns[10:]]
                    
                    #Extracting the spot with the maxima area
                    intersect_ = intersect.loc[:, intersect.columns.str.contains('area*')]                    
                    bigger_polygon = intersect_.idxmax().values[0]
                    
                    #Preparing Geodataframe - renaming columns
                    interest_polygon = intersect.loc[bigger_polygon].to_frame()
                    name_columns = ['intersection_percentage', 'time', 'Tb', 'area_tb', 
                                    'centroides', 'mean_tb', 'min_tb', 'max_tb', 'mean_pp','max_pp', 'belong', 
                                    'intersection_percentage', 'geometry']
                    interest_polygon_t = interest_polygon.T
                    interest_polygon_t.columns = name_columns

                    #Area of the bigger spot                
                    area_bigger_polygon = interest_polygon_t.area_tb.values[0]
                    
                    covex =  interest_polygon
                    covex = gpd.GeoDataFrame(interest_polygon.T, geometry=interest_polygon.T.geometry)
                    interest_polygon =  covex.geometry.apply(lambda g: Polygon(g))
                    interest_polygon = interest_polygon.T
                    
                    #Estimating the distance between centroids                 
                    dist = min_dist(interest_polygon, poli_ti1)
                                        
                    #Estimating percentage of intersection                
                    intersect_percentage = ((intersect.area/1000000)*100/area_bigger_polygon)[0].round(1)                

                    #Condition if threshold overlapping percentage is activated
                    #If the spot fails with threshold overlapping percentage is dropped                    
                    if isinstance(threshold_overlapping_percentage, checknumbers.Real):
                        if intersect_percentage>= threshold_overlapping_percentage:
                            sup.at[isites, "intersection_percentage"] = intersect_percentage 
                            #print ("There was a " + str(intersect_percentage) + "% intersection") 
                        else:
                            sup.drop(isites, inplace = True)
                            #print ("The convective system " + str(isites) + " was dropped to the data due intersection is below that the threshold: "+ str(intersect_percentage)) 
                    #Condition if threshold overlapping percentage is not activated
                    else:
                        sup.at[isites, "intersection_percentage"] = intersect_percentage 
                        #print ("There was a " + str(intersect_percentage) + "% intersection") 
                    #print ("Warning: in the tracking of this convective system were more than 1 intersecting polygon: "  + str(isites))

                #Attaching found spot to the track
                time_intersect = dist.to_frame().T["time"].values[0]
                dist = tb_intersect = dist.to_frame().T["Tb"].values[0]
                sup.loc[(sup.time == time_intersect) & (sup.Tb == tb_intersect), "belong"] = msc_counter
                msc_counter +=1
                #print ("The convective system: " + str(isites) + " - belongs to the track: " +str(int(msc_counter)))
            except:
               msc_counter +=1

        # Checking when the record (spot) has been processed or attaching to a track            
        elif sup.at[isites, "belong"] in range(1, msc_counter):
            try:
                #Generating start and end time date (time t and time t+1)                
                date_initial = sup.iloc[isites].time
                date_final = date_initial + timedelta(hours=1)
            
                #Reference spot
                poli_ti = sup.loc[isites].to_frame().T
                #Establishing EPSG:32718 - UTM zone 18S plane coordinate system                 
                poli_ti = gpd.GeoDataFrame(poli_ti, geometry='geometry', crs=utm_local_zone)
                #Getting the id track previously estimate
                index_track = poli_ti.belong
                
                #Spots one hour later
                poli_ti1 = sup.loc[(sup.time == date_final) & (sup.belong == 0)]
    
                #Intersection of reference spot with spots one hour later
                intersect = gpd.overlay(poli_ti, poli_ti1, how='intersection')
                
                try:#When there is only one intersection spot
                    
                    #Estimating the distance between centroids
                    dist = min_dist(intersect.centroid, poli_ti1)
                    
                    #Extracting the spot with the maxima area
                    intersect_area_max = max(intersect.area_tb_1.values, intersect.area_tb_2.values)
    
                    #Estimating percentage of intersection                    
                    intersect_percentage = ((intersect.area/1000000)*100/intersect_area_max)[0].round(1)                
                    #print ("There was a " + str(intersect_percentage) + "% intersection") 
                    
                    #Condition if threshold overlapping percentage is activated
                    #If the spot fails with threshold overlapping percentage is dropped
                    if isinstance(threshold_overlapping_percentage, checknumbers.Real):
                        if intersect_percentage>= threshold_overlapping_percentage:
                            sup.at[isites, "intersection_percentage"] = intersect_percentage 
                            #print ("There was a " + str(intersect_percentage) + "% intersection") 
                        else:
                            sup.drop(isites, inplace = True)
                            #print ("The convective system " + str(isites) + " was dropped to the data due intersection is below that the threshold: "+ str(intersect_percentage)) 
                    
                    #Condition if threshold overlapping percentage is not activated
                    else:
                        sup.at[isites, "intersection_percentage"] = intersect_percentage 
                        #print ("There was a " + str(intersect_percentage) + "% intersection") 

                except:  #When there are more than one intersection spot is selected the bigger spot
                   
                    #dimensioning the Geodataframe columns only for the intersected spots
                    intersect = intersect[intersect.columns[10:]]
                    
                    #Extracting the spot with the maxima area
                    intersect_ = intersect.loc[:, intersect.columns.str.contains('area*')]                    
                    bigger_polygon = intersect_.idxmax().values[0]
                                        
                    #Preparing Geodataframe - renaming columns                    
                    interest_polygon = intersect.loc[bigger_polygon].to_frame()
                    name_columns = ['intersection_percentage', 'time', 'Tb', 'area_tb', 
                                    'centroides', 'mean_tb', 'min_tb', 'max_tb', 'mean_pp','max_pp', 'belong', 
                                    'intersection_percentage', 'geometry']
                    interest_polygon_t = interest_polygon.T
                    interest_polygon_t.columns = name_columns
                    
                    #Area of the bigger spot                
                    area_bigger_polygon = interest_polygon_t.area_tb.values[0]

                    #Estimating the distance between centroids                                                                                            
                    intersect = intersect.loc[bigger_polygon].to_frame().T
                    dist = min_dist(intersect.centroid, poli_ti1)
                
                    #Estimating percentage of intersection                    
                    intersect_percentage = ((intersect.area/1000000)*100/area_bigger_polygon)[0].round(1)                

                    #Condition if threshold overlapping percentage is activated
                    #If the spot fails with threshold overlapping percentage is dropped                    
                    if isinstance(threshold_overlapping_percentage, checknumbers.Real):
                        if intersect_percentage>= threshold_overlapping_percentage:
                            sup.at[isites, "intersection_percentage"] = intersect_percentage 
                            #print ("There was a " + str(intersect_percentage) + "% intersection") 
                        else:
                            sup.drop(isites, inplace = True)
                            #print ("The convective system " + str(isites) + " was dropped to the data due intersection is below that the threshold: "+ str(intersect_percentage)) 
                    #Condition if threshold overlapping percentage is not activated
                    else:
                        sup.at[isites, "intersection_percentage"] = intersect_percentage 
                        #print ("There was a " + str(intersect_percentage) + "% intersection") 
                #print ("Warning: in the tracking of this convective system were more than 1 intersecting polygon: "  + str(isites))

                #Attaching found spot to the track
                time_intersect = dist.to_frame().T["time"].values[0]
                tb_intersect = dist.to_frame().T["Tb"].values[0]
                sup.loc[(sup.time == time_intersect) & (sup.Tb == tb_intersect), "belong"] = index_track.values[0]
                #print ("The convective system: " + str(isites) + " - belongs to the track: " +str(int(index_track.values[0])))
            except:
                pass   
            
        #Progress bar
        time.sleep(0.01)    
       
    #Transforming plane coordinates to geodesics  coordinates   
    sup["centroid_"] = sup["centroid_"].to_crs(4326)
    sup = sup.to_crs(4326) 
    
    #Creating an original index
    sup["id_gdf"] = None
    sup["id_gdf"] = sup.index
    
    sup = TRACKS(sup)
    
    if path_save is None:
        pass
    else:
        try:
            sup.save(path_save)
        except:
            raise FileNotFoundError("No such file or directory: " + str(path_save))
    
    
    return sup.data
