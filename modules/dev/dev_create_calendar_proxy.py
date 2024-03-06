# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 11:30:04 2022

@author: MialykO
"""
import modules.acea.acea_core as ac
import numpy as np
from netCDF4 import Dataset
import os
from osgeo import gdal

def Covert5to30arcmin(data):
    cols30 = 720; rows30 = int(cols30/2) # Latitude and Longitude for 30 arcmin grid
    output = np.ones((rows30, cols30))*np.ma.masked
    for col in range(cols30):
        for row in range(rows30):
            [x5, y5] = [int(col)*6, int(row)*6]
            x_range5 = [x5, x5+6]; y_range5 = [y5, y5+6]
            if (data.mask[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]==False).any():
                output[row,col] = 1
    return output

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
    data_source_used.long_name = '[1] GGCMI crop calendar (doi: 10.5281/zenodo.5062513)'; 
    
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

#%% Crop calendar&
proxy_crop_code_name = 'aff2'
final_crop_code_name = 'aff'; crop_fao = 641 #!!!!!!!!!!
f_growing_season = 1 # growing season factor relative to proxy
minimum_gs_duration = 180
maximum_gs_duration = 300
output_folder = r'O:/Oleks/ACEA_v2/data/acea/crop_calendar'
cutoff_cells = True

# Read rainfed
proxy_plant_day_rf = ac.ReadNC(f"O:/Oleks/ACEA_v2/data/acea/crop_calendar/{proxy_crop_code_name}_rf_crop_calendar.nc", 'planting_day')
proxy_gs_rf = ac.ReadNC(f"O:/Oleks/ACEA_v2/data/acea/crop_calendar/{proxy_crop_code_name}_rf_crop_calendar.nc", 'growing_season_length')
if cutoff_cells: # limit to actual growing areas
    cutoff_mask = gdal.Open(f"O:/Oleks/ACEA_v2/data/acea/harvested_areas/{crop_fao}/spam2010V2r0_global_H_{crop_fao}_R.tif")
    cutoff_mask = cutoff_mask.GetRasterBand(1).ReadAsArray()
    cutoff_mask = Covert5to30arcmin(np.ma.masked_where(cutoff_mask <= 0, cutoff_mask))
    proxy_plant_day_rf[cutoff_mask.mask] = np.ma.masked
    proxy_gs_rf[cutoff_mask.mask] = np.ma.masked
    
new_gs_rf = proxy_gs_rf * f_growing_season
new_gs_rf[new_gs_rf<minimum_gs_duration] = minimum_gs_duration; new_gs_rf[new_gs_rf>maximum_gs_duration] = maximum_gs_duration
new_harvest_day_rf = proxy_plant_day_rf + new_gs_rf
new_harvest_day_rf[new_harvest_day_rf > 365] = new_harvest_day_rf[new_harvest_day_rf > 365]-365

source_data = np.ones((proxy_plant_day_rf.shape))
path = ac.GetPath([output_folder,f'{final_crop_code_name}_rf_crop_calendar.nc'])
Create2DRaster(path, np.round(proxy_plant_day_rf), np.round(new_harvest_day_rf), source_data)


# Read irrigated
proxy_plant_day_ir = ac.ReadNC(f"O:/Oleks/ACEA_v2/data/acea/crop_calendar/{proxy_crop_code_name}_ir_crop_calendar.nc", 'planting_day')
proxy_gs_ir = ac.ReadNC(f"O:/Oleks/ACEA_v2/data/acea/crop_calendar/{proxy_crop_code_name}_ir_crop_calendar.nc", 'growing_season_length')
if cutoff_cells: # limit to actual growing areas
    cutoff_mask = gdal.Open(f"O:/Oleks/ACEA_v2/data/acea/harvested_areas/{crop_fao}/spam2010V2r0_global_H_{crop_fao}_I.tif")
    cutoff_mask = cutoff_mask.GetRasterBand(1).ReadAsArray()
    cutoff_mask = Covert5to30arcmin(np.ma.masked_where(cutoff_mask <= 0, cutoff_mask))
    proxy_plant_day_ir[cutoff_mask.mask] = np.ma.masked
    proxy_gs_ir[cutoff_mask.mask] = np.ma.masked
    
new_gs_ir = proxy_gs_ir * f_growing_season
new_gs_ir[new_gs_ir<minimum_gs_duration] = minimum_gs_duration; new_gs_ir[new_gs_ir>maximum_gs_duration] = maximum_gs_duration
new_harvest_day_ir = proxy_plant_day_ir + new_gs_ir
new_harvest_day_ir[new_harvest_day_ir > 365] = new_harvest_day_ir[new_harvest_day_ir > 365]-365

source_data = np.ones((proxy_plant_day_ir.shape))
path = ac.GetPath([output_folder,f'{final_crop_code_name}_ir_crop_calendar.nc'])
Create2DRaster(path, np.round(proxy_plant_day_ir), np.round(new_harvest_day_ir), source_data)