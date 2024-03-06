# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 13:59:08 2022

@author: MialykO
"""
#%% Import packages
import numpy as np
import modules.acea.acea_core as ac
import pandas as pd

#%% Functions
def GetListOfHistCountries():
    'Lest all the current countries'
    
    with open(ac.GetPath(["data", "footprints",'regions','FAOSTAT_hist_countries.txt']), 'rb') as f:
        country_data = pd.read_csv(f, sep='\t', lineterminator='\r')

    return country_data

def GetListOfCurCountries():
    'Lest all the current countries'
    
    with open(ac.GetPath(["data", "acea",'grid','regions_countries.txt']), 'rb') as f:
        country_data = pd.read_csv(f, sep='\t', lineterminator='\r')

    return country_data

def GetHistCountryCells5(faostat_code):
    "Select grid cells of a country at 5 arc minutes"
    
    # For current countries
    m49_countries = ac.ReadNC(ac.GetPath(["data", "acea",'grid',"countries_m49_5arc.nc"]))
    id5 = ac.ReadNC(ac.GetPath(["data", "acea",'grid',"id5.nc"]))
    
    with open(ac.GetPath(["data", "footprints",'regions',"FAOSTAT_hist_countries.txt"]), 'rb') as f:
        _country_data = pd.read_csv(f, sep='\t', lineterminator='\r')
    m49_code = _country_data.loc[_country_data['Country Code'] == faostat_code]['M49 Code'].values[0]
    row5_list, col5_list = np.where(m49_countries==m49_code); id5_list = id5[row5_list,col5_list]
    country_cells = pd.DataFrame(np.array([id5_list,row5_list,col5_list]).astype(int).transpose(), columns=['id5', 'rowy5', 'colx5'])
    
    # For no longer existing countries (e.g. USSR) 
    if len(country_cells) == 0:
        with open(ac.GetPath(["data", "acea",'grid',"faostat_country_id5.csv"]), 'rb') as f: hist_id5_faostat = pd.read_csv(f, sep='\t', lineterminator='\r').dropna().astype(int)
        cells_id5 = hist_id5_faostat[hist_id5_faostat['faostat']==faostat_code]['id5'].to_list()
        row5_list, col5_list = np.where(np.isin(id5, cells_id5)); id5_list = id5[row5_list,col5_list]
        country_cells = pd.DataFrame(np.array([id5_list,row5_list,col5_list]).astype(int).transpose(), columns=['id5', 'rowy5', 'colx5'])
    
    return country_cells

def GetSimulatedCrops(crop_group='all'):
    "Get the list with simulated crops"
    
    path = ac.GetPath(["data", "footprints", 'Simulated_crops.txt'])
    data = pd.read_csv(path, sep='\t')
    
    if crop_group == 'all':
        crops = data
    else:
        crops = data[data['Group']==crop_group]
    
    return crops

def GetHistHarvestedAreas5(crop_code, irrigated):
    "Get historical harvested areas"
    
    if not irrigated:
        harv_area = ac.ReadNC(ac.GetPath(["data", 'footprints', "harvested_areas",crop_code, f'{crop_code}_area_rf_harv_global_annual_1990_2019.nc']), 'harvested_area')
        harv_area[np.isnan(harv_area)] = np.ma.masked
    else:
        harv_area = ac.ReadNC(ac.GetPath(["data", 'footprints', "harvested_areas",crop_code, f'{crop_code}_area_ir_harv_global_annual_1990_2019.nc']), 'harvested_area')
    
    harv_area[np.isnan(harv_area)] = np.ma.masked
    return harv_area

def GetIrrSystems5(consider_all=True):
    "Get the fraction of irrigation systems per irrigated cell in 5 acr minute arrays"
    
    ir_syst_data = ac.ReadNC(ac.GetPath(["data", 'footprints', 'cell_share_irrigation_systems_5arcmin.nc']), '3-dimensional Scientific Dataset')[:,:,:3] # Surface and Sprinkler fractions. Reference: https://doi.org/10.5194/hess-19-3073-2015
    ir_syst_data[ir_syst_data<0] = 0
    
    if consider_all:
        ir_syst_data[np.sum(ir_syst_data,2)==0] = [1, 0, 0]; sys_sum = np.sum(ir_syst_data,2)
        ir_syst_data_sur = ir_syst_data[:,:,0]/sys_sum
        ir_syst_data_spr = ir_syst_data[:,:,1]/sys_sum
        ir_syst_data_drp = ir_syst_data[:,:,2]/sys_sum
            	# Convert to system fraction, assume surface if nothing is reported
    else: # Only surface
        ir_syst_data_sur = np.ones((ir_syst_data.shape[0],ir_syst_data.shape[1]))
        ir_syst_data_spr = ir_syst_data_sur*0
        ir_syst_data_drp = ir_syst_data_sur*0
    
    return ir_syst_data_sur, ir_syst_data_spr, ir_syst_data_drp