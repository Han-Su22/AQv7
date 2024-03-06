# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 17:41:24 2022

@author: MialykO
"""
import modules.acea.acea_core as ac
import numpy as np
from netCDF4 import Dataset
import os
from osgeo import gdal
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
    data_source_used.long_name = '[1] GGCMI crop calendar (doi: 10.5281/zenodo.5062513) and own rule-based map'; 
    
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

#%%
dominant_calendar = r"O:\Oleks\ACEA_v2\data\acea\crop_calendar\aff2_ir_crop_calendar.nc"
background_calendar = r"O:\Oleks\ACEA_v2\data\acea\crop_calendar\pea_ggcmi_ir_crop_calendar.nc"
output_folder = r"O:\Oleks\ACEA_v2\data\acea\crop_calendar"

# Dominant
dm_plant_day = ac.ReadNC(dominant_calendar, 'planting_day')
dm_harvest_day = ac.ReadNC(dominant_calendar, 'maturity_day')
# dm_harvest_day[dm_plant_day==1] = np.ma.masked; dm_plant_day[dm_plant_day==1] = np.ma.masked

# Background
bk_plant_day = ac.ReadNC(background_calendar, 'planting_day')
bk_harvest_day = ac.ReadNC(background_calendar, 'maturity_day')

# Final
dm_plant_day[dm_plant_day.mask==True] = bk_plant_day[dm_plant_day.mask==True]
dm_harvest_day[dm_harvest_day.mask==True] = bk_harvest_day[dm_harvest_day.mask==True]

source_data = np.ones((dm_plant_day.shape))
path = ac.GetPath([output_folder,'aff_ir_crop_calendar.nc'])
Create2DRaster(path, np.round(dm_plant_day), np.round(dm_harvest_day), source_data)