# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 14:59:04 2021

@author: MialykO
"""
from psycopg2 import DatabaseError, connect
import pickle
import pandas as pd
import sys
import numpy as np
import modules.acea.acea_core as ac

#%% Database
# Database configurations, _ for internal variables
_db_name = 'ctw_wemwm_postgis'
_db_host = 'postgresql-wemwm.utwente.nl'
_db_user = 'ctw_wemwm_postgis'
_db_password = 'iR7Vookiexie'
_db_port = 5432

def Dev_GetCroplandCells30():
    "Select global croplands cells at 30 arc minutes"
    
    try:
        # Connent to database
        con = connect(database=_db_name, host=_db_host, user=_db_user, password=_db_password, port=_db_port)
        gridcells = pd.read_sql("SELECT c.* FROM ACEA.land_mask m INNER JOIN(SELECT id30, colx, rowy, st_y(centroid) lat FROM grid.grid30_world) c ON m.id30=c.id30", con)
        con.close()
    except DatabaseError as e: print(f'Error {e}'); sys.exit(1) # Check for errors
        
    with open(ac.GetPath(["data", "acea",'cropland_mask_30arc.pckl']), 'wb') as f:
        pickle.dump(gridcells, f) # (id30,x,y,lat)
    
# def Dev_GetCroplandCells5(save_path):
#     "Select global croplands cells at 5 arc minutes"
    
#     try:
#         # Connent to database
#         con = connect(database=_db_name, host=_db_host, user=_db_user, password=_db_password, port=_db_port); cur = con.cursor()
#         cur.execute("SELECT c.* FROM ACEA.land_mask m INNER JOIN(SELECT id5,id30, colx, rowy, st_y(centroid) lat FROM grid.grid5_world) c ON m.id30=c.id30")
#         cells = cur.fetchall(); con.close()
#     except DatabaseError as e: print(f'Error {e}'); sys.exit(1) # Check for errors
        
#     gridcells = np.zeros(len(cells), dtype=[('id5',np.int32),('x',np.int32),('y',np.int32),('lat',np.float32)])
#     for j in range(len(cells)): gridcells[j] = (cells[j][0], cells[j][1], cells[j][2], cells[j][3])
        
#     f = open(save_path, 'wb'); pickle.dump(gridcells, f); f.close()
def Dev_GetHistCountryCells5():
    "Get id5 by FAOSTAT country code "
    
    try: # Connent to database
        con = connect(database=_db_name, host=_db_host, user=_db_user, password=_db_password, port=_db_port)
        id5 = pd.read_sql("SELECT DISTINCT c.id5, c.faostat FROM grid.grid5_countries_lastknown_extent c", con)
        con.close()
    except DatabaseError as e: print(f'Error {e}'); sys.exit(1) # Check for errors
    
    id5.to_csv('faostat_country_id5.csv', sep='\t', index=False)
    
def Dev_GetCurrentCountryCells5():
    "Get 5 arc minute country cells by FAOSTAT country code (based on 2021 country map)"
    
    try: # Connent to database
        con = connect(database=_db_name, host=_db_host, user=_db_user, password=_db_password, port=_db_port)
        gridcells = pd.read_sql("SELECT DISTINCT c.id5, g5.id30, g5.colx-1 colx5, g5.rowy-1 rowy5, g30.colx-1 colx30, g30.rowy-1 rowy30, c.iso1al3 iso1al3, c.faostat \
                                    FROM grid.grid5_countries_lastknown_extent c \
                                    INNER JOIN(SELECT id5, id30, colx, rowy FROM grid.grid5_world) g5 ON g5.id5=c.id5 \
                                    INNER JOIN(SELECT id30, colx, rowy FROM grid.grid30_world) g30 ON g30.id30=g5.id30 \
                                    WHERE c.faostat IN (1,2,3,4,6,7,8,9,10,11,12,13,14,16,18,19,20,21,23,25,26,27,28,29,32,33,35,37,38,39,40,41,44,45,46,47,48,49,50,52,53,54,55,56,57,58,59,60,61,63,64,66,67,68,72,73,74,75,79,80,81,83,84,86,89,90,91,93,95,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,115,116,117,118,119,120,121,122,123,124,126,127,129,130,131,132,133,134,136,137,138,140,141,143,144,145,146,147,148,149,150,154,155,156,157,158,159,160,162,165,166,167,168,169,170,171,173,174,175,176,178,179,180,181,183,184,185,188,189,191,192,193,194,195,196,197,198,199,200,201,202,203,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,225,226,227,229,230,231,233,234,235,236,237,238,244,249,250,251,255,256,272,273,276,277)", con)
        con.close()
    except DatabaseError as e: print(f'Error {e}'); sys.exit(1) # Check for errors
    
    parameters = ['id5', 'id30', 'colx30', 'rowy30', 'faostat', 'colx5', 'rowy5']
    title = 'Country map'
    path = ac.GetPath(["data", "acea",'grid','current_countries_5arc.nc'])
    values5 = np.ones((7,2160,4320)) * np.ma.masked
    
    for index, cell in gridcells.iterrows():
        if index % 10000 == 0:
            print(f'Cells read: {index}')
        row5 = int(cell['rowy5']); col5 = int(cell['colx5'])
        for var in range(len(parameters)):
            values5[var, row5, col5] = cell[parameters[var]]

    ac.CreateLayeredRaster(values5, parameters, path, title)
    # with open(GetPath(["data", "acea",'country_mask_5arc.pckl']), 'wb') as f:
    #     pickle.dump(gridcells, f) # (id30,x,y,lat)
    
def GetFAOtoFreshFactors():
    "Get official FAO dry to fresh convertion factors by a FAOSTAT crop code"
    # Connent to database
    con = connect(database=_db_name, host=_db_host, user=_db_user, password=_db_password, port=_db_port)
    crop_ratios = pd.read_sql("SELECT faostat_code, faostat_to_dry_yield FROM crop.crop_groupings ", con)
    
    file_name = ac.GetPath(["data", "acea",'crop_parameters', 'crop_dry_to_fresh_ratios.npz'])
    with open(file_name, 'wb') as f:
        np.savez_compressed(f, crop_ratios=crop_ratios) 

def GetMesfinHoekstraAvg():
    "Get Mesfin& Hoekstra data from 47 Report"
    # Connent to database
    con = connect(database=_db_name, host=_db_host, user=_db_user, password=_db_password, port=_db_port)
    crop_cwu = pd.read_sql("SELECT DISTINCT r.crop_code, rg.cwu_m3ha/10 as green, rb.cwu_m3ha/10 as blue, rgr.cwu_m3ha/10 as grey FROM report47.cwu_allcrops_globalaverage as r \
                                INNER JOIN(SELECT * FROM report47.cwu_allcrops_globalaverage WHERE color = 'green') rg ON rg.crop_code=r.crop_code \
                                INNER JOIN(SELECT * FROM report47.cwu_allcrops_globalaverage WHERE color = 'blue') rb ON rb.crop_code=r.crop_code \
                                INNER JOIN(SELECT * FROM report47.cwu_allcrops_globalaverage WHERE color = 'grey') rgr ON rgr.crop_code=r.crop_code \
                                ORDER BY r.crop_code ASC", con)
                                    
    file_name = ac.GetPath(['mesfin_hoekstra_rep47_cwu_avg.npz'])
    with open(file_name, 'wb') as f:
        np.savez_compressed(f, crop_cwu=crop_cwu) 
#%% Call functions
# Dev_GetCroplandCells30()
# Dev_GetCurrentCountryCells5()
Dev_GetHistCountryCells5()
# GetFAOtoFreshFactors()
# GetMesfinHoekstraAvg()