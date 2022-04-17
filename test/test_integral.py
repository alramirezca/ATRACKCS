# -*- coding: utf-8 -*-
"""
test algorithm tb adn tb_p parameterization
@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
"""
#Libraries
from atrackcs.utils import funcs
from atrackcs.identify import identify_mcs
from atrackcs.features import features_Tb, features_P, features_Tracks
from atrackcs.track import track_mcs
import pandas as pnd
import geopandas as gpd
from shapely.wkt import loads
from pandas.util.testing import assert_frame_equal
import geopandas as gpd
from shapely.wkt import loads
import re



UTM_LOCAL_ZONE = 32718 
UTC_LOCAL_HOUR = 5

phTb= r"C:\Users\ASUS\Documents\GitHub\atrackcs\test\csv/tb/"
phP= r"C:\Users\ASUS\Documents\GitHub\atrackcs\test\csv/p/"

#___________________________________Functions______________________________________________________

#Significant decimals
simpledec = re.compile(r"\d*\.\d+")
def mround(match):
    return "{:.5f}".format(float(match.group()))


def run_tb(phTb, Tb, area_Tb):
    ds = funcs.readNC(pathTb = phTb, pathP =None, utc_local_hour = UTC_LOCAL_HOUR)

    sup = identify_mcs(ds, variables = "Tb", Tb = Tb, area_Tb = area_Tb, 
                       utm_local_zone = UTM_LOCAL_ZONE, path_save = None)
        
    sup = features_Tb(sup, ds)
    
    sup = track_mcs(sup, threshold_overlapping_percentage = None, 
                    utm_local_zone = UTM_LOCAL_ZONE, 
                    path_save = None)
    
    resume = features_Tracks(sup, encrypt_index = False, path_save = None)
    # resume.duracion_total =  resume.duracion_total.astype(int)
    # resume.distancia_total =  resume.distancia_total.astype(float)
    # resume.velocidad_promedio =  resume.velocidad_promedio.astype(float)
    
    return resume

def run_tb_p(phTb, phP, Tb, area_Tb, min_precipitation, area_P):
    ds = funcs.readNC(pathTb = phTb, pathP =phP, utc_local_hour = UTC_LOCAL_HOUR)

    sup = identify_mcs(ds, variables = "Both", Tb = Tb, area_Tb = area_Tb, 
                       utm_local_zone = UTM_LOCAL_ZONE, path_save = None)
        
    sup = features_Tb(sup, ds)
    sup = features_P(sup, ds, min_precipitation = min_precipitation, 
                     area_P = area_P, drop_empty_precipitation = False)

    sup = track_mcs(sup, threshold_overlapping_percentage = None, 
                    utm_local_zone = UTM_LOCAL_ZONE, 
                    path_save = None)
    
    resume = features_Tracks(sup, encrypt_index = False, path_save = None)
    # resume.duracion_total =  resume.duracion_total.astype(int)
    # resume.distancia_total =  resume.distancia_total.astype(float)
    # resume.velocidad_promedio =  resume.velocidad_promedio.astype(float)
    
    return resume


def test_run_tb():
    
    df = pnd.read_csv(r"C:/Users/ASUS/Documents/GitHub/atrackcs/test/csv/resultados/track_test.csv",
                                   index_col = ["belong", "id_gdf"], parse_dates = ["time"])
    df['geometry'] = gpd.GeoSeries.from_wkt(df['geometry'])
    df.geometry = df.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    df['centroid_'] = gpd.GeoSeries.from_wkt(df['centroid_'])
    df['centroid_'] = df['centroid_'].apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    expected_tracks = gpd.GeoDataFrame(df, geometry='geometry', crs = 4326)
    
    tracks = run_tb(phTb, 225, 1000)
    tracks.geometry = tracks.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))
    tracks.centroid_ = tracks.centroid_.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))
   
    
    assert_frame_equal(tracks.sort_index(),expected_tracks.sort_index(), check_dtype = False)   
    
def test_run_tb_p(): 
    df = pnd.read_csv(r"C:/Users/ASUS/Documents/GitHub/atrackcs/test/csv/resultados/track_test_with_p.csv",
                                   index_col = ["belong", "id_gdf"], parse_dates = ["time"])
    df['geometry'] = gpd.GeoSeries.from_wkt(df['geometry'])
    df.geometry = df.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    df['centroid_'] = gpd.GeoSeries.from_wkt(df['centroid_'])
    df.centroid_ = df.centroid_.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    expected_tracks = gpd.GeoDataFrame(df, geometry='geometry', crs = 4326)
    
    tracks = run_tb_p(phTb, phP, 225, 1000, 2, 500)
    tracks.geometry = tracks.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))
    tracks.centroid_ = tracks.centroid_.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))
   
    
    assert_frame_equal(tracks.sort_index(),expected_tracks.sort_index(), check_dtype = False)   
    

