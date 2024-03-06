# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 14:57:58 2022

@author: MialykO
"""
from datetime import datetime
import modules.acea.acea_core as ac
import numpy as np
from netCDF4 import Dataset
import os

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
    data_source_used.long_name = '[1] MIRCA2000 (https://doi.org/10.1029/2008GB003435)'; 
    
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

#%% From MIRCA
input_file_rf = r"C:\Users\MialykO\surfdrive\Data\Crop\Masks\MIRCA\-1_sub1_rf_MIRCA_data.nc"
input_file_ir = r"C:\Users\MialykO\surfdrive\Data\Crop\Masks\MIRCA\-1_sub1_ir_MIRCA_data.nc"
output_folder = r'C:\Work\Oleks\2nd paper\MIRCA'
cols30 = 720; rows30 = int(cols30/2) # Latitude and Longitude for 30 arcmin grid

# Magic
mirca_rf = ac.ReadNC(input_file_rf, 'value')
mirca_rf_plant = mirca_rf[1,:,:]; mirca_rf_harvest = mirca_rf[2,:,:]
mirca_ir = ac.ReadNC(input_file_ir, 'value')
mirca_ir_plant = mirca_ir[1,:,:]; mirca_ir_harvest = mirca_ir[2,:,:]
dt1 = datetime(day=1, month=1, year=2000)

plant_day_rf = np.ones((rows30, cols30))*np.ma.masked; harvest_day_rf=plant_day_rf*1
plant_day_ir = plant_day_rf*1; harvest_day_ir=plant_day_rf*1

                    
                    
for col in range(cols30):
    for row in range(rows30):
        [x5, y5] = [int(col)*6, int(row)*6]
        x_range5 = [x5, x5+6]; y_range5 = [y5, y5+6]
        
        if (mirca_rf_plant.mask[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]==False).any():
            _month_plant = mirca_rf_plant[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
            _month_plant = int(np.median(_month_plant[_month_plant.mask==False].data))
            _month_harv = mirca_rf_harvest[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
            _month_harv = int(np.median(_month_harv[_month_harv.mask==False].data))
            
            plant_day_rf[row,col] = (datetime(day=15, month=_month_plant, year=2000)-dt1).days
            harvest_day_rf[row,col] = (datetime(day=15, month=_month_harv, year=2000)-dt1).days
        
        if (mirca_ir_plant.mask[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]==False).any():
            _month_plant = mirca_ir_plant[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
            _month_plant = int(np.median(_month_plant[_month_plant.mask==False].data))
            _month_harv = mirca_ir_harvest[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
            _month_harv = int(np.median(_month_harv[_month_harv.mask==False].data))
            
            plant_day_ir[row,col] = (datetime(day=15, month=_month_plant, year=2000)-dt1).days
            harvest_day_ir[row,col] = (datetime(day=15, month=_month_harv, year=2000)-dt1).days


source_data = np.ones((plant_day_rf.shape))
path = ac.GetPath([output_folder,'annual_rf.nc'])
Create2DRaster(path, plant_day_rf, harvest_day_rf, source_data)

path = ac.GetPath([output_folder,'annual_ir.nc'])
Create2DRaster(path, plant_day_ir, harvest_day_ir, source_data)



