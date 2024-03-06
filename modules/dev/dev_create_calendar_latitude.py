# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:23:42 2022

@author: MialykO
"""

import modules.acea.acea_core as ac
import numpy as np
import pickle
import pandas as pd
import os
import datetime
from netCDF4 import Dataset
from functools import partial

def Create2DRaster(file_path, plt_data, harv_data, source_data):
    "Create NETCDF rasters"
    
    cols30 = 720; rows30 = int(cols30/2) # Latitude and Longitude for 30 arcmin grid
    if os.path.exists(file_path):os.remove(file_path)
    
    raster = Dataset(file_path,'w', format='NETCDF4')
    raster.title = 'Crop calendar'
    raster.institution = 'University of Twente, Netherlands'
    raster.contact = 'Oleksandr Mialyk o.mialyk@utwente.nl'
    
    # Dimensions    
    raster.createDimension('lat', rows30);raster.createDimension('lon', cols30)

    # Variables
    lon = raster.createVariable('lon', 'f8', 'lon'); lat = raster.createVariable('lat', 'f8', 'lat')
    planting_day = raster.createVariable('planting_day', 'f4', ('lat','lon',), zlib = True, complevel = 6, fill_value = 1.e20)
    maturity_day = raster.createVariable('maturity_day', 'f4', ('lat','lon',), zlib = True, complevel = 6, fill_value = 1.e20)
    growing_season_length = raster.createVariable('growing_season_length', 'f4', ('lat','lon',), zlib = True, complevel = 6, fill_value = 1.e20)
    data_source_used = raster.createVariable('data_source_used', 'f4', ('lat','lon',), zlib = True, complevel = 6, fill_value = 1.e20)
    planting_day.units = 'day of year'; maturity_day.units = 'day of year'; growing_season_length.units = 'days'
    data_source_used.long_name = '[1] Rule-based calendar'; 
    
    # Adding values
    incrm = 0.5
    lat[:] = np.arange(90-incrm/2, -rows30/4, -incrm)
    lat.long_name = 'lat'; lat.units = 'degrees_north'
    lon[:] = np.arange(-180+incrm/2, cols30/4, incrm)
    lon.long_name = 'lon'; lon.units = 'degrees_east'
    
    gs_len = harv_data-plt_data
    gs_len[gs_len < 0] = gs_len[gs_len < 0]+365
    
    planting_day[:,:] = plt_data; maturity_day[:,:] = harv_data; growing_season_length[:,:] = gs_len; data_source_used[:,:] = source_data
    
    raster.close()
    
    
# Input data
crop_calendar_path = r'O:/Oleks/crop_calendar_latitude.nc'
growing_season_length = 360
row_switch = 180
planting_days = [1, 180] # [above row_switch, below]

gridcells30 = ac.GetAllCells30()
idxs = gridcells30.index.to_numpy()

plant_days = np.ones((360, 720))*np.ma.masked; harvest_days = plant_days*1

for idx in idxs:
    row = gridcells30.loc[idx]
    cell_id = row['id30']; row30 = row['rowy30']; col30 = row['colx30']
    if row30 <= row_switch:
        plant_days[row30, col30] = planting_days[0]
    else:
        plant_days[row30, col30] = planting_days[1]
    harvest_days[row30, col30] = plant_days[row30, col30] + growing_season_length
    
print(f'Number of run cells final: {harvest_days.count()}')
Create2DRaster(crop_calendar_path, plant_days, harvest_days, np.ones((harvest_days.shape)))