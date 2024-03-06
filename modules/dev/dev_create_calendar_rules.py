# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 11:47:41 2022

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

def CropCalendar(calendar, idx):
    
    row = calendar.gridcells30.loc[idx]
    # row = calendar.gridcells30.loc[calendar.gridcells30['id30'] == 63023]
    # if cell_count % 500 == 0: print(f'Cells read: {cell_count} out of {len(gridcells30)} ({cell_count/len(gridcells30)*100:.2f} %)')
    cell_id = row['id30']; row30 = row['rowy30']; col30 = row['colx30']
    weather_file_path = f'{calendar.climate_folder}{int(cell_id)}.pckl'
    plant = False; harvest = False
    # if not row['id30'] == 82484: continue
    # if row30 < 100: continue
    if os.path.exists(weather_file_path): # Check if file exists
        with open(weather_file_path, 'rb') as f: tmax, tmin,_,_ = pickle.load(f) # Read climate
        tmean = pd.Series((tmax + tmin)/2).rolling(10).mean().to_numpy() # Get moving avg og mean temperature
        sow_days = []; harv_days = []; sim_day = np.where(calendar.day_idxs == datetime.datetime(1990,1,1))[0][0]; gs = False
        
        while sim_day < len(calendar.day_idxs):
            if (tmean[sim_day]>calendar.budbreak_temp)\
            and (np.median(tmean[sim_day:sim_day+20])>calendar.budbreak_temp+5)\
            and not gs:
                if len(harv_days)>0 and (harv_days[-1]+calendar.gs_break)>sim_day: sim_day+=1; continue
                else: sow_days.append(sim_day); gs = True
            if gs\
            and sim_day > (calendar.mig_gs_len + sow_days[-1])\
            and (tmean[sim_day]<calendar.harvest_temp)\
            and (np.median(tmean[sim_day:sim_day+10])<calendar.harvest_temp):
                harv_days.append(sim_day)
                gs = False
            
            sim_day +=1
        
        if len(sow_days)>0 and len(harv_days)>0:
            plant = int(np.median([d.timetuple().tm_yday for d in calendar.day_idxs[sow_days]]))
            harvest = int(np.median([d.timetuple().tm_yday for d in calendar.day_idxs[harv_days]]))
            
    return [row30,col30,plant,harvest]

def mycallback(res):
    for i in range(len(res)):
        if res[i][2] and res[i][3]:
            plant_days[res[i][0], res[i][1]] = res[i][2]
            harvest_days[res[i][0], res[i][1]] = res[i][3]

class calendar:
    # Input data
    climate_folder = r"O:/Oleks/ACEA_v2/data/acea/climate/gswp3-w5e5_obsclim/gswp3-w5e5_obsclim_01-01-1981_31-12-2019_"
    crop_calendar_path = r'O:/Oleks/crop_calendar_rules_7_15_1501.nc'
    
    # Set of rules
    budbreak_temp = 7 # to determine the sowing/budbreak
    harvest_temp = 15 # to determine the harvest
    mig_gs_len = 150 # minimum number of days in gs
    gs_break = 30 # break between growing seasons
    
    def __init__(self):
        self.gridcells30 = ac.GetAllCells30()
        self.day_idxs = pd.date_range(start='01/01/1981', end='31/12/2019').to_pydatetime()
        
        
    def Run(self):
        idxs = self.gridcells30.index.to_numpy()
        CropCalendar_partial = partial(CropCalendar, self)
        
        # a = CropCalendar(self, 1)
        
        
        from multiprocessing import get_context
        with get_context("spawn").Pool(14) as p:
            r = p.map_async(CropCalendar_partial, idxs, callback=mycallback)
            r.wait()

        print(f'Number of run cells final: {harvest_days.count()}')
        Create2DRaster(self.crop_calendar_path, plant_days, harvest_days, np.ones((harvest_days.shape)))

plant_days = np.ones((360, 720))*np.ma.masked; harvest_days = plant_days*1
if __name__ == '__main__':
    a = calendar()
    a.Run()
