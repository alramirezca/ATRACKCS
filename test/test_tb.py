# -*- coding: utf-8 -*-
"""
test algorithm tb parameterization
@authors: Álvaro Ramírez Cardona (alramirezca@unal.edu.co)
          Vanessa Robledo Delgado (vanessa.robledo@udea.edu.co)
"""

from atrackcs.utils import funcs
from atrackcs.identify import identify_mcs
from atrackcs.features import features_Tb, features_P, features_Tracks
from atrackcs.track import track_mcs
import pandas as pnd
import geopandas as gpd
from shapely.wkt import loads
from pandas.util.testing import assert_frame_equal



UTM_LOCAL_ZONE = 32718 
UTC_LOCAL_HOUR = 5

phTb= "C:/Users/ASUS/Desktop/udea/Prueba_GPM/prueba_test_2/tb/"
phP= "C:/Users/ASUS/Desktop/udea/Prueba_GPM/prueba_test_2/p/"


def run_tb(phTb, Tb, area_Tb):
    ds = funcs.readNC(pathTb = phTb, pathP =None, utc_local_hour = UTC_LOCAL_HOUR)

    sup = identify_mcs(ds, variables = "Tb", Tb = Tb, area_Tb = area_Tb, 
                       utm_local_zone = UTM_LOCAL_ZONE, path_save = None)
        
    sup = features_Tb(sup, ds)
    
    sup = track_mcs(sup, threshold_overlapping_percentage = None, 
                    utm_local_zone = UTM_LOCAL_ZONE, 
                    path_save = None)
    
    resume = features_Tracks(sup, encript_index = False, path_save = None)
    resume.duracion_total =  resume.duracion_total.astype(int)
    resume.distancia_total =  resume.distancia_total.astype(float)
    resume.velocidad_promedio =  resume.velocidad_promedio.astype(float)
    
    return resume

import geopandas as gpd
from shapely.wkt import loads
import re

simpledec = re.compile(r"\d*\.\d+")
def mround(match):
    return "{:.5f}".format(float(match.group()))


def test_run_tb():
    
    df = pnd.read_csv("C:/Users/ASUS/Desktop/udea/Prueba_GPM/prueba_test_2/resultados/track_test.csv",
                                   index_col = ["belong", "id_gdf"], parse_dates = ["time"])
    df['geometry'] = gpd.GeoSeries.from_wkt(df['geometry'])
    df.geometry = df.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    df['centroides'] = gpd.GeoSeries.from_wkt(df['centroides'])
    df.centroides = df.centroides.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))

    expected_tracks = gpd.GeoDataFrame(df, geometry='geometry', crs = 4326)
    
    tracks = run_tb(phTb, 225, 1000)
    tracks.geometry = tracks.geometry.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))
    tracks.centroides = tracks.centroides.apply(lambda x: loads(re.sub(simpledec, mround, x.wkt)))
   
    
    assert_frame_equal(tracks.sort_index(),expected_tracks.sort_index(), check_dtype = False)   
    
    

