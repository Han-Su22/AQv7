# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 14:59:04 2021

@author: MialykO
"""
from netCDF4 import Dataset
import datetime
import pickle
import numpy as np
import modules.acea.acea_core as ac

class Dev_CropParameters:
    "Get crop parameters"
    
    def __init__(self, conf):
        self.conf = conf
        self.data_folder = ac.GetPath(["data", "acea"])
        self.clim_end = conf.climate_end; self.num_years = 30
        self.clim_start = self.clim_end-self.num_years
        self.crop = self.CropDB(self.conf.crop_fao)
        self.gridcells30 = ac.GetAllCells30()
        
    
    def CreateGDDraster(self, irrigated):
        _gridcells = ac.GetSimulationCells30(self.conf.gridcells)
        final_data = np.ones((11,360,720)) *np.ma.masked
        first_year =  datetime.datetime(self.conf.climate_start,1,1)
        errors = []
        
        # Get crop calendar
        if irrigated:
            print(f'Generating crop parameters for irrigated {self.conf.crop_name}')
            self.crop_planting = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_ir_crop_calendar.nc']), 'planting_day')
            self.crop_harvest = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_ir_crop_calendar.nc']), 'maturity_day')
        else:
            print(f'Generating crop parameters for rainfed {self.conf.crop_name}')
            self.crop_planting  = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_rf_crop_calendar.nc']), 'planting_day')
            self.crop_harvest = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_rf_crop_calendar.nc']), 'maturity_day')
        
        cell_count = 0
        for cell_id in _gridcells: # Check which cells are valid
            cell_count += 1
            if cell_count % 1000 == 0: print(f'Cells read: {cell_count} out of {len(_gridcells)} ({cell_count/len(_gridcells)*100:.2f} %)')
                
            cell = self.gridcells30[self.gridcells30['id30']==cell_id]
            rowy = int(cell['rowy30']); colx = int(cell['colx30'])
            
            if self.crop_planting[rowy, colx] > 0:
                # Get climate data
                try:
                    with open(ac.GetPath([self.data_folder,'climate', f'{self.conf.climate_name}_{cell_id}.pckl']), 'rb') as f:
                        tmax, tmin,_,_ = pickle.load(f)
                except:
                    errors = np.append(errors, cell)
                    continue
            
                # Define sowing/harvest day
                sowing_day = int(np.ceil(self.crop_planting[rowy, colx]))
                if self.crop_harvest[rowy, colx] > 0: harvest_day = int(self.crop_harvest[rowy, colx]) # Check if the harvest day is provided
                else: harvest_day = 0 # No harvest date is provided
                
                # Calculate GDDs
                gs_len = np.zeros(self.num_years)
                gdd = gs_len*1;tmax_avg = gs_len*1; tmin_avg = gs_len*1
                GDD_accumulated = np.cumsum(self.ComputeGDD(tmax, tmin))
                
                for gs in range(self.num_years): # Iterate growing seasons
                    year_days = (datetime.datetime(self.clim_start+gs+1,1,1) - datetime.datetime(self.clim_start+gs,1,1)).days # Number of days in the current year
                    idx_start = (datetime.datetime(self.clim_start+gs,1,1)-first_year).days + sowing_day - 1 # Get sowing day index
                    if gs > 1 and idx_start + np.mean(gs_len) > len(tmax): continue # Skip if expected harvest date is beyond available temperature data
                    
                    if harvest_day == 0:
                        harvest_day = self.CalcHarvestDay(GDD_accumulated[idx_start:idx_start+year_days] - GDD_accumulated[idx_start-1], sowing_day, year_days)
                        reset_hd = True
                    else: reset_hd = False
                        
                    # Get harvest day index
                    if harvest_day > sowing_day: idx_end = idx_start + (harvest_day - sowing_day)
                    elif idx_start + harvest_day + year_days - sowing_day < len(tmax): idx_end = idx_start + (harvest_day + year_days - sowing_day)# For winter crops
                    else: continue # Skip if idx_end is beyond the available data
                    
                    # Get accumulated GDDs
                    if idx_start == 0: temp = 0 # Don't substract if sowing day index of the first growing season is 0
                    else: temp = GDD_accumulated[idx_start-1]
                        
                    gdd[gs] = GDD_accumulated[idx_end] - temp
                    tmax_avg[gs] = np.nanmean(tmax[idx_start-1:idx_end])
                    tmin_avg[gs] = np.nanmean(tmin[idx_start-1:idx_end])
                    
                    # Calculate growing season length (needed if harvest day is not provided)
                    if sowing_day > harvest_day:
                        gs_len[gs] = harvest_day + year_days - sowing_day
                    else:
                        gs_len[gs] = harvest_day - sowing_day
                    
                    # Recalculate harvest date for the cells where this date is not initially provided
                    if reset_hd: harvest_day = 0
                
                # Average growing season GGDs
                if any(gdd!=0): avg_GDD = np.nanmean(gdd[gdd!=0])
                else: avg_GDD = self.crop['gdd_maturity']/2 # If the environment is very cold, take a half of the original value
                coef = avg_GDD/self.crop['gdd_maturity'] # Adjustment coefficient for crop phenology
                if coef<.5: coef = .5
                if coef>1.5: coef = 1.5
                
                # Canopy growth coefficient
                CCini = self.crop['plants_per_ha'] * self.crop['seedling_cover'] * 10**-8
                tgrowth = -1/self.crop['cgc'] * np.log(0.08*CCini/self.crop['CCmax']) * coef
                cgc = -1/tgrowth * np.log(0.08*CCini/self.crop['CCmax'])
                
                # Canopy decline coefficient
                CCend = (np.exp(self.crop['cdc'] * (self.crop['gdd_maturity'] - self.crop['gdd_senescence'])/self.crop['CCmax'])-21)/-20 * self.crop['CCmax']
                cdc = self.crop['CCmax']/(self.crop['gdd_maturity'] * coef - self.crop['gdd_senescence'] * coef) * np.log(21-20*CCend/self.crop['CCmax'])
                               
                # Define avg harvest day
                # avg_gs = np.nanmean(gs_len[gs_len!=0]) # Avg growing season
                # harvest_day = sowing_day + avg_gs # (needed if harvest day is not provided)
                # if harvest_day > 365: harvest_day = np.ceil(harvest_day - 365) # Ceil to avoid 0
                
                # Define final crop parameters
                final_data[0, rowy, colx] = int(self.crop['gdd_emergence'] * coef) # Growing degree/Calendar days from sowing to emergence/transplant recovery
                final_data[1, rowy, colx] = int(self.crop['gdd_max_root'] * coef) # Growing degree/Calendar days from sowing to maximum rooting
                final_data[2, rowy, colx] = int(self.crop['gdd_senescence'] * coef) # Growing degree/Calendar days from sowing to senescence
                final_data[3, rowy, colx] = int(self.crop['gdd_maturity'] * coef) # Growing degree/Calendar days from sowing to maturity
                final_data[4, rowy, colx] = int(self.crop['gdd_yield_form'] * coef) # Growing degree/Calendar days from sowing to start of yield formation
                final_data[5, rowy, colx] = int(self.crop['gdd_duration_flowering'] * coef) # Duration of flowering in growing degree/calendar days (-999 for non-fruit/grain crops)
                final_data[6, rowy, colx] = int(self.crop['gdd_duration_yield_form'] * coef) # Duration of yield formation in growing degree/calendar days
                final_data[7, rowy, colx] = float(cdc) # Canopy decline coefficient (fraction per GDD/calendar day)
                final_data[8, rowy, colx] = float(cgc) # Canopy growth coefficient (fraction per GDD)
                final_data[9, rowy, colx] = np.nanmean(tmax_avg) # Max temperature used to adjust heat stress
                final_data[10, rowy, colx] = np.nanmean(tmin_avg) # Min temperature used to adjust heat stress
        # Save rasters
        if irrigated:
            _path = ac.GetPath([self.data_folder,'crop_parameters',  f'{self.conf.crop_name_short}_ir_crop_parameters.nc'])
            self.CreateParamRaster30(final_data, _path, 'Irrigated crop parameters')
        else:
            _path = ac.GetPath([self.data_folder,'crop_parameters',  f'{self.conf.crop_name_short}_rf_crop_parameters.nc'])
            self.CreateParamRaster30(final_data, _path, 'Rainfed crop parameters')
            
    def ComputeGDD(self, tmax, tmin):
        "Calculate Growing Degree Days (according to AquaCrop's method 3)"
        
        # Calculate the adjusted Tmin and Tmax based on Tbase and Tupper
        tmax[tmax > self.crop['tupper']] = self.crop['tupper']
        tmax[tmax < self.crop['tbase']] = self.crop['tbase']
        tmin[tmin > self.crop['tupper']] = self.crop['tupper']
        
        # Calculate Tavg and the growing degree days for each day
        tavg = np.round((tmax+tmin)/2, 1)
        tavg[tavg < self.crop['tbase']] = self.crop['tbase']
        
        return tavg - self.crop['tbase']
    
    def CalcHarvestDay(self, GDD_accumulated, sowing_day, year_days): 
        """ Determines the harvest day based on the accumulated Growing Degree Days 
        (according to AquaCrop's method 3) since the day of sowing. """
        
        harvest_day = sowing_day + np.argmax(GDD_accumulated > self.crop['gdd_maturity'])
        
        if harvest_day > year_days: harvest_day = harvest_day - year_days # If a winter crop
                              
        # If maturity can't be reached or harvest is close to sowing, harvest 2 weeks before the next sowing
        if sowing_day >= harvest_day >= sowing_day - 14:
            harvest_day = sowing_day - 14 
            # print('Harvest was triggered 14 days before the next sowing')
        
        if harvest_day == 366: harvest_day = 365
        
        return harvest_day
        
    def CropDB(self,crop_fao):
        'Get needed crop parameters'
        crop = dict()
        if crop_fao == 27: # Rice
            crop['tupper'] = 30; crop['tbase'] = 8
            crop['plants_per_ha'] = 1_000_000
            crop['seedling_cover'] = 6.
            crop['cgc'] = .007; crop['cdc'] = .005
            crop['CCmax'] = .95
            crop['gdd_emergence'] = 50
            crop['gdd_max_root'] = 370
            crop['gdd_senescence'] = 1300
            crop['gdd_maturity'] = 1900
            crop['gdd_yield_form']  = 1150
            crop['gdd_duration_flowering'] = 350
            crop['gdd_duration_yield_form'] = 680
        elif crop_fao == 44: # barley
            crop['tupper'] = 28; crop['tbase'] = 2
            crop['plants_per_ha'] =1_500_000
            crop['seedling_cover'] = 1.5
            crop['cgc'] = 0.008697; crop['cdc'] = 0.006
            crop['CCmax'] = 0.8
            crop['gdd_emergence'] = 98
            crop['gdd_max_root'] = 854
            crop['gdd_senescence'] = 924
            crop['gdd_maturity'] = 1296
            crop['gdd_yield_form']  = 867
            crop['gdd_duration_flowering'] = 160
            crop['gdd_duration_yield_form'] = 351
        elif crop_fao == 75: # Oats
            crop['tupper'] = 30; crop['tbase'] = 0
            crop['plants_per_ha'] =1_600_000
            crop['seedling_cover'] = 1.5
            crop['cgc'] = 0.00591; crop['cdc'] = 0.00291
            crop['CCmax'] = 0.98
            crop['gdd_emergence'] = 132
            crop['gdd_max_root'] = 890
            crop['gdd_senescence'] = 890
            crop['gdd_maturity'] = 1308
            crop['gdd_yield_form']  = 684
            crop['gdd_duration_flowering'] = 177
            crop['gdd_duration_yield_form'] = 559
        elif crop_fao == 79: # Millet
            crop['tupper'] = 33; crop['tbase'] = 10
            crop['plants_per_ha'] = 180_000
            crop['seedling_cover'] = 5
            crop['cgc'] = 0.02091; crop['cdc'] = 0.008
            crop['CCmax'] = 0.95
            crop['gdd_emergence'] = 61
            crop['gdd_max_root'] = 1583
            crop['gdd_senescence'] = 1579
            crop['gdd_maturity'] = 1760
            crop['gdd_yield_form']  = 930
            crop['gdd_duration_flowering'] = 306
            crop['gdd_duration_yield_form'] = 789
        elif crop_fao == 108: # Cereals
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 180_000
            crop['seedling_cover'] = 5
            crop['cgc'] = 0.02; crop['cdc'] = 0.008
            crop['CCmax'] = 0.9
            crop['gdd_emergence'] = 60
            crop['gdd_max_root'] = 1400
            crop['gdd_senescence'] = 1500
            crop['gdd_maturity'] = 1700
            crop['gdd_yield_form']  = 900
            crop['gdd_duration_flowering'] = 300
            crop['gdd_duration_yield_form'] = 800
        elif crop_fao == 116: # Potato
            crop['tupper'] = 26; crop['tbase'] = 2
            crop['plants_per_ha'] = 58_000
            crop['seedling_cover'] = 10
            crop['cgc'] = .0149; crop['cdc'] = .0038
            crop['CCmax'] = .94
            crop['gdd_emergence'] = 332
            crop['gdd_max_root'] = 967
            crop['gdd_senescence'] = 1468
            crop['gdd_maturity'] = 2324
            crop['gdd_yield_form']  = 553
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 1748
        elif crop_fao == 122: # Sweet potato
            crop['tupper'] = 35; crop['tbase'] = 15
            crop['plants_per_ha'] = 40_000
            crop['seedling_cover'] = 15
            crop['cgc'] = .00966; crop['cdc'] = .00798
            crop['CCmax'] = .94
            crop['gdd_emergence'] = 77
            crop['gdd_max_root'] = 538
            crop['gdd_senescence'] = 1091
            crop['gdd_maturity'] = 1294
            crop['gdd_yield_form']  = 415
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 875
        elif crop_fao == 125: # Cassava
            crop['tupper'] = 30; crop['tbase'] = 10
            crop['plants_per_ha'] = 12_500
            crop['seedling_cover'] = 10
            crop['cgc'] = .00948; crop['cdc'] = .0037
            crop['CCmax'] = .88
            crop['gdd_emergence'] = 110
            crop['gdd_max_root'] = 770
            crop['gdd_senescence'] = 3300
            crop['gdd_maturity'] = 3960
            crop['gdd_yield_form']  = 880
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 2750
        elif crop_fao == 136: # Taro
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 27_750
            crop['seedling_cover'] = 25
            crop['cgc'] = .007; crop['cdc'] = .0058
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 460
            crop['gdd_max_root'] = 1000
            crop['gdd_senescence'] = 2115
            crop['gdd_maturity'] = 2406
            crop['gdd_yield_form']  = 1512
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 861
        elif crop_fao == 137: # Yams
            crop['tupper'] = 30; crop['tbase'] = 15
            crop['plants_per_ha'] = 8_800
            crop['seedling_cover'] = 10
            crop['cgc'] = .00948; crop['cdc'] = .0037
            crop['CCmax'] = .91
            crop['gdd_emergence'] = 54
            crop['gdd_max_root'] = 378
            crop['gdd_senescence'] = 1621
            crop['gdd_maturity'] = 1945
            crop['gdd_yield_form']  = 432
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 1351
        elif crop_fao == 156: # Sugarcane
            crop['tupper'] = 32; crop['tbase'] = 12
            crop['plants_per_ha'] = 140_000
            crop['seedling_cover'] = 6.5
            crop['cgc'] = .01087; crop['cdc'] = .006
            crop['CCmax'] = .95
            crop['gdd_emergence'] = 82
            crop['gdd_max_root'] = 700
            crop['gdd_senescence'] = 3851
            crop['gdd_maturity'] = 4259
            crop['gdd_yield_form']  = 0
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 852
        elif crop_fao == 157: # Sugar beet
            crop['tupper'] = 25; crop['tbase'] = 3
            crop['plants_per_ha'] = 100_000
            crop['seedling_cover'] = 1
            crop['cgc'] = .010541; crop['cdc'] = .003857
            crop['CCmax'] = .98
            crop['gdd_emergence'] = 23
            crop['gdd_max_root'] = 408
            crop['gdd_senescence'] = 1704
            crop['gdd_maturity'] = 2203
            crop['gdd_yield_form']  = 865
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 852
        elif crop_fao == 176: # Beans
            crop['tupper'] = 30; crop['tbase'] = 9
            crop['plants_per_ha'] = 131_579
            crop['seedling_cover'] = 10
            crop['cgc'] = .009879; crop['cdc'] = .008813
            crop['CCmax'] = .99
            crop['gdd_emergence'] = 59
            crop['gdd_max_root'] = 888
            crop['gdd_senescence'] = 903
            crop['gdd_maturity'] = 1298
            crop['gdd_yield_form']  = 556
            crop['gdd_duration_flowering'] = 233
            crop['gdd_duration_yield_form'] = 668
        elif crop_fao == 187: # Pea
            crop['tupper'] = 27; crop['tbase'] = 5
            crop['plants_per_ha'] = 810_000
            crop['seedling_cover'] = 5
            crop['cgc'] = .00949; crop['cdc'] = .00209
            crop['CCmax'] = .95
            crop['gdd_emergence'] = 144
            crop['gdd_max_root'] = 581
            crop['gdd_senescence'] = 875
            crop['gdd_maturity'] = 942
            crop['gdd_yield_form']  = 685
            crop['gdd_duration_flowering'] = 162
            crop['gdd_duration_yield_form'] = 215
        elif crop_fao == 191: # Chickpea
            crop['tupper'] = 32; crop['tbase'] = 8
            crop['plants_per_ha'] = 330_000
            crop['seedling_cover'] = 5
            crop['cgc'] = .009879; crop['cdc'] = .0053
            crop['CCmax'] = .93
            crop['gdd_emergence'] = 144
            crop['gdd_max_root'] = 900
            crop['gdd_senescence'] = 1392
            crop['gdd_maturity'] = 1597
            crop['gdd_yield_form']  = 672
            crop['gdd_duration_flowering'] = 350
            crop['gdd_duration_yield_form'] = 750
        elif crop_fao == 195: # Cowpea
            crop['tupper'] = 36; crop['tbase'] = 10
            crop['plants_per_ha'] = 200_000
            crop['seedling_cover'] = 10
            crop['cgc'] = .0097; crop['cdc'] = .0053
            crop['CCmax'] = .93
            crop['gdd_emergence'] = 88
            crop['gdd_max_root'] = 831
            crop['gdd_senescence'] = 1101
            crop['gdd_maturity'] = 1259
            crop['gdd_yield_form']  = 916
            crop['gdd_duration_flowering'] = 217
            crop['gdd_duration_yield_form'] = 288
        elif crop_fao == 197: # Pigeon pea
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 810_000
            crop['seedling_cover'] = 5
            crop['cgc'] = .0095; crop['cdc'] = .0021
            crop['CCmax'] = .95
            crop['gdd_emergence'] = 150
            crop['gdd_max_root'] = 600
            crop['gdd_senescence'] = 900
            crop['gdd_maturity'] = 1000
            crop['gdd_yield_form']  = 700
            crop['gdd_duration_flowering'] = 150
            crop['gdd_duration_yield_form'] = 250
        elif crop_fao == 201: # Lentils
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 130_000
            crop['seedling_cover'] = 10
            crop['cgc'] = .01; crop['cdc'] = .0088
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 60
            crop['gdd_max_root'] = 900
            crop['gdd_senescence'] = 900
            crop['gdd_maturity'] = 1300
            crop['gdd_yield_form']  = 550
            crop['gdd_duration_flowering'] = 200
            crop['gdd_duration_yield_form'] = 700
        elif crop_fao == 211: # Pulses
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 130_000
            crop['seedling_cover'] = 10
            crop['cgc'] = .01; crop['cdc'] = .0088
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 60
            crop['gdd_max_root'] = 900
            crop['gdd_senescence'] = 900
            crop['gdd_maturity'] = 1300
            crop['gdd_yield_form']  = 550
            crop['gdd_duration_flowering'] = 200
            crop['gdd_duration_yield_form'] = 700
        elif crop_fao == 217: # Cashew
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 230_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .002; crop['cdc'] = .0015
            crop['CCmax'] = .65
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 4000
            crop['gdd_maturity'] = 4200
            crop['gdd_yield_form']  = 400
            crop['gdd_duration_flowering'] = 1500
            crop['gdd_duration_yield_form'] = 3200
        elif crop_fao == 242: # Groundnuts (peanut)
            crop['tupper'] = 30; crop['tbase'] = 9
            crop['plants_per_ha'] = 88_889
            crop['seedling_cover'] = 3
            crop['cgc'] = .00942; crop['cdc'] = .00683
            crop['CCmax'] = .7
            crop['gdd_emergence'] = 127
            crop['gdd_max_root'] = 1040
            crop['gdd_senescence'] = 1592
            crop['gdd_maturity'] = 1850
            crop['gdd_yield_form']  = 595
            crop['gdd_duration_flowering'] = 798
            crop['gdd_duration_yield_form'] = 943
        elif crop_fao == 249: # coconut
            crop['tupper'] = 40; crop['tbase'] = 15
            crop['plants_per_ha'] = 390_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .0018; crop['cdc'] = .0014
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 5844
            crop['gdd_maturity'] = 5966
            crop['gdd_yield_form']  = 1461
            crop['gdd_duration_flowering'] = 2192
            crop['gdd_duration_yield_form'] = 4444
        elif crop_fao == 254: # Oil palm
            crop['tupper'] = 40; crop['tbase'] = 15
            crop['plants_per_ha'] = 390_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .0018; crop['cdc'] = .0014
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 5844
            crop['gdd_maturity'] = 5966
            crop['gdd_yield_form']  = 1461
            crop['gdd_duration_flowering'] = 2192
            crop['gdd_duration_yield_form'] = 4444
        elif crop_fao == 260: # Olives
            crop['tupper'] = 31; crop['tbase'] = 6
            crop['plants_per_ha'] = 250_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .003234; crop['cdc'] = .00119
            crop['CCmax'] = .5
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 1426
            crop['gdd_maturity'] = 3185
            crop['gdd_yield_form']  = 856
            crop['gdd_duration_flowering'] = 285
            crop['gdd_duration_yield_form'] = 1720
        elif crop_fao == 263: # Sheanuts
            crop['tupper'] = 40; crop['tbase'] = 15
            crop['plants_per_ha'] = 390_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .0018; crop['cdc'] = .0014
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 5800
            crop['gdd_maturity'] = 5950
            crop['gdd_yield_form']  = 1450
            crop['gdd_duration_flowering'] = 2200
            crop['gdd_duration_yield_form'] = 4400
        elif crop_fao == 267: # Sunflower
            crop['tupper'] = 30; crop['tbase'] = 4
            crop['plants_per_ha'] = 57_000
            crop['seedling_cover'] = 5
            crop['cgc'] = 0.014993; crop['cdc'] = 0.006
            crop['CCmax'] = 0.98
            crop['gdd_emergence'] = 170
            crop['gdd_max_root'] = 1784
            crop['gdd_senescence'] = 1900
            crop['gdd_maturity'] = 2400
            crop['gdd_yield_form']  = 1266
            crop['gdd_duration_flowering'] = 350
            crop['gdd_duration_yield_form'] = 1087
        elif crop_fao == 270: # Rapeseed
            crop['tupper'] = 30; crop['tbase'] = 0
            crop['plants_per_ha'] = 440_000
            crop['seedling_cover'] = 5
            crop['cgc'] = .0074; crop['cdc'] = .00433
            crop['CCmax'] = .8
            crop['gdd_emergence'] = 240
            crop['gdd_max_root'] = 840
            crop['gdd_senescence'] = 1993
            crop['gdd_maturity'] = 2473
            crop['gdd_yield_form']  = 1308
            crop['gdd_duration_flowering'] = 528
            crop['gdd_duration_yield_form'] = 1050
        elif crop_fao == 289: # Sesame
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 300_000
            crop['seedling_cover'] = 4.5
            crop['cgc'] = .01025; crop['cdc'] = .00908
            crop['CCmax'] = .95
            crop['gdd_emergence'] = 154
            crop['gdd_max_root'] = 1148
            crop['gdd_senescence'] = 1246
            crop['gdd_maturity'] = 1470
            crop['gdd_yield_form']  = 574
            crop['gdd_duration_flowering'] = 224
            crop['gdd_duration_yield_form'] = 784
        elif crop_fao == 339: # Oilseeds
            crop['tupper'] = 30; crop['tbase'] = 0
            crop['plants_per_ha'] = 440_000
            crop['seedling_cover'] = 5
            crop['cgc'] = .0075; crop['cdc'] = .0043
            crop['CCmax'] = .8
            crop['gdd_emergence'] = 170
            crop['gdd_max_root'] = 900
            crop['gdd_senescence'] = 1900
            crop['gdd_maturity'] = 2400
            crop['gdd_yield_form']  = 1200
            crop['gdd_duration_flowering'] = 350
            crop['gdd_duration_yield_form'] = 1100
        elif crop_fao == 388: # Tomato
            crop['tupper'] = 28; crop['tbase'] = 7
            crop['plants_per_ha'] = 33_333
            crop['seedling_cover'] = 20
            crop['cgc'] = .007504; crop['cdc'] = .004
            crop['CCmax'] = .75
            crop['gdd_emergence'] = 43
            crop['gdd_max_root'] = 891
            crop['gdd_senescence'] = 1553
            crop['gdd_maturity'] = 1933
            crop['gdd_yield_form']  = 525
            crop['gdd_duration_flowering'] = 750
            crop['gdd_duration_yield_form'] = 1050
        elif crop_fao == 402: # Pepper
            crop['tupper'] = 30; crop['tbase'] = 7
            crop['plants_per_ha'] = 65_000
            crop['seedling_cover'] = 15
            crop['cgc'] = .00681; crop['cdc'] = .004
            crop['CCmax'] = .75
            crop['gdd_emergence'] = 60
            crop['gdd_max_root'] = 420
            crop['gdd_senescence'] = 1600
            crop['gdd_maturity'] = 1680
            crop['gdd_yield_form']  = 402
            crop['gdd_duration_flowering'] = 915
            crop['gdd_duration_yield_form'] = 1270
        elif crop_fao == 403: # Onion
            crop['tupper'] = 30; crop['tbase'] = 10
            crop['plants_per_ha'] = 330_000
            crop['seedling_cover'] = 5
            crop['cgc'] = .0075; crop['cdc'] = .0053
            crop['CCmax'] = .65
            crop['gdd_emergence'] = 60
            crop['gdd_max_root'] = 343
            crop['gdd_senescence'] = 1263
            crop['gdd_maturity'] = 1450
            crop['gdd_yield_form']  = 816
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 634
        elif crop_fao == 463: # Vegetables
            crop['tupper'] = 30; crop['tbase'] = 7
            crop['plants_per_ha'] = 50_000
            crop['seedling_cover'] = 20
            crop['cgc'] = .0075; crop['cdc'] = .004
            crop['CCmax'] = .8
            crop['gdd_emergence'] = 40
            crop['gdd_max_root'] = 900
            crop['gdd_senescence'] = 1500
            crop['gdd_maturity'] = 1900
            crop['gdd_yield_form']  = 500
            crop['gdd_duration_flowering'] = 700
            crop['gdd_duration_yield_form'] = 1000
        elif crop_fao == 10000: # Vegetables, leafy
            crop['tupper'] = 30; crop['tbase'] = 7
            crop['plants_per_ha'] = 80_000
            crop['seedling_cover'] = 7
            crop['cgc'] = .0075; crop['cdc'] = .004
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 40
            crop['gdd_max_root'] = 900
            crop['gdd_senescence'] = 1800
            crop['gdd_maturity'] = 1900
            crop['gdd_yield_form']  = 0
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 1800
        elif crop_fao == 10001: # Vegetables, tuber
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 100_000
            crop['seedling_cover'] = 1
            crop['cgc'] = .0105; crop['cdc'] = .0039
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 40
            crop['gdd_max_root'] = 400
            crop['gdd_senescence'] = 1300
            crop['gdd_maturity'] = 1900
            crop['gdd_yield_form']  = 700
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 1100
        elif crop_fao == 486: # Banana
            crop['tupper'] = 40; crop['tbase'] = 15
            crop['plants_per_ha'] = 410_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .002; crop['cdc'] = .0015
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 4100
            crop['gdd_maturity'] = 4200
            crop['gdd_yield_form']  = 2500
            crop['gdd_duration_flowering'] = 800
            crop['gdd_duration_yield_form'] = 1500
        elif crop_fao == 489: # Plantain
            crop['tupper'] = 40; crop['tbase'] = 15
            crop['plants_per_ha'] = 410_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .002; crop['cdc'] = .0015
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 4100
            crop['gdd_maturity'] = 4200
            crop['gdd_yield_form']  = 2500
            crop['gdd_duration_flowering'] = 800
            crop['gdd_duration_yield_form'] = 1500
        elif crop_fao == 490: # Oranges
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 390000
            crop['seedling_cover'] = 200
            crop['cgc'] = .002; crop['cdc'] = .0015
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 2959
            crop['gdd_maturity'] = 4000
            crop['gdd_yield_form']  = 822
            crop['gdd_duration_flowering'] = 300
            crop['gdd_duration_yield_form'] = 3178
        elif crop_fao == 515: # Apples
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 280_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .0045; crop['cdc'] = .0021
            crop['CCmax'] = .8
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 1550
            crop['gdd_maturity'] = 2650
            crop['gdd_yield_form']  = 400
            crop['gdd_duration_flowering'] = 300
            crop['gdd_duration_yield_form'] = 2000
        elif crop_fao == 558: # Berries
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 403_333
            crop['seedling_cover'] = 15
            crop['cgc'] = .006; crop['cdc'] = .002
            crop['CCmax'] = .7
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 1500
            crop['gdd_maturity'] = 1900
            crop['gdd_yield_form']  = 500
            crop['gdd_duration_flowering'] = 700
            crop['gdd_duration_yield_form'] = 1000
        elif crop_fao == 560: # Grapes
            crop['tupper'] = 30; crop['tbase'] = 10
            crop['plants_per_ha'] = 403_333
            crop['seedling_cover'] = 15
            crop['cgc'] = .0054; crop['cdc'] = .0014
            crop['CCmax'] = .65
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 1941
            crop['gdd_maturity'] = 3681
            crop['gdd_yield_form']  = 700
            crop['gdd_duration_flowering'] = 300
            crop['gdd_duration_yield_form'] = 1400
        elif crop_fao == 572: # Avocado
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 390_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .002; crop['cdc'] = .0015
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 4000
            crop['gdd_maturity'] = 4200
            crop['gdd_yield_form']  = 400
            crop['gdd_duration_flowering'] = 800
            crop['gdd_duration_yield_form'] = 3200
        elif crop_fao == 603: # Tropical fruits
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 410_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .002; crop['cdc'] = .0015
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 4000
            crop['gdd_maturity'] = 4200
            crop['gdd_yield_form']  = 400
            crop['gdd_duration_flowering'] = 600
            crop['gdd_duration_yield_form'] = 1500
        elif crop_fao == 641: # alfalfa
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 2_000_000
            crop['seedling_cover'] = 40
            crop['cgc'] = .0075; crop['cdc'] = .0048
            crop['CCmax'] = .97
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 2200
            crop['gdd_maturity'] = 2300
            crop['gdd_yield_form'] = 0
            crop['gdd_duration_flowering'] = 0
            crop['gdd_duration_yield_form'] = 2200
        elif crop_fao == 656: # Coffee
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 390000
            crop['seedling_cover'] = 200
            crop['cgc'] = .0018; crop['cdc'] = .0014
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 2800
            crop['gdd_maturity'] = 3100
            crop['gdd_yield_form']  = 900
            crop['gdd_duration_flowering'] = 300
            crop['gdd_duration_yield_form'] = 2000
        elif crop_fao == 661: # Cocoa
            crop['tupper'] = 35; crop['tbase'] = 10
            crop['plants_per_ha'] = 390_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .002; crop['cdc'] = .0015
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 5844
            crop['gdd_maturity'] = 5966
            crop['gdd_yield_form']  = 1461
            crop['gdd_duration_flowering'] = 2192
            crop['gdd_duration_yield_form'] = 4444
        elif crop_fao == 667: # Tea
            crop['tupper'] = 32; crop['tbase'] = 8
            crop['plants_per_ha'] = 220_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .007; crop['cdc'] = .0035
            crop['CCmax'] = .9
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 2800
            crop['gdd_maturity'] = 3100
            crop['gdd_yield_form']  = 900
            crop['gdd_duration_flowering'] = 300
            crop['gdd_duration_yield_form'] = 2000
        elif crop_fao == 821: # Fibres
            crop['tupper'] = 30; crop['tbase'] = 5
            crop['plants_per_ha'] = 250_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .0032; crop['cdc'] = .0012
            crop['CCmax'] = .5
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 1420
            crop['gdd_maturity'] = 3200
            crop['gdd_yield_form']  = 850
            crop['gdd_duration_flowering'] = 280
            crop['gdd_duration_yield_form'] = 1700
        elif crop_fao == 826: # Tobacco
            crop['tupper'] = 30; crop['tbase'] = 7
            crop['plants_per_ha'] = 20_000
            crop['seedling_cover'] = 40
            crop['cgc'] = .0075; crop['cdc'] = .004
            crop['CCmax'] = .75
            crop['gdd_emergence'] = 50
            crop['gdd_max_root'] = 900
            crop['gdd_senescence'] = 1850
            crop['gdd_maturity'] = 2100
            crop['gdd_yield_form']  = 700
            crop['gdd_duration_flowering'] = 600
            crop['gdd_duration_yield_form'] = 1200
        elif crop_fao == 836: # Rubber
            crop['tupper'] = 40; crop['tbase'] = 15
            crop['plants_per_ha'] = 250_000
            crop['seedling_cover'] = 200
            crop['cgc'] = .004463; crop['cdc'] = .002096
            crop['CCmax'] = .85
            crop['gdd_emergence'] = 20
            crop['gdd_max_root'] = 20
            crop['gdd_senescence'] = 1561
            crop['gdd_maturity'] = 2657
            crop['gdd_yield_form']  = 368
            crop['gdd_duration_flowering'] = 290
            crop['gdd_duration_yield_form'] = 2060
        
        else:
            raise ValueError('Non existent crop in the CropDB')
        return crop
    
    def CreateParamRaster30(self, values, path, title):
        "Create 30 arc minute raster with many bands"
        
        parameters = ['gdd_emergence','gdd_max_root','gdd_senescence','gdd_maturity','gdd_yield_form',\
              'gdd_duration_flowering','gdd_duration_yield_form','cdc','cgc','temp_max_avg','temp_min_avg']
            
        cols30 = 720; rows30 = int(cols30/2) # Latitude and Longitude for 5 arcmin grid
        par_num = values.shape[0]
        
        # Create a dataset
        raster = Dataset(path,'w', format='NETCDF4') #'w' stands for write
        raster.title = 'AquaCrop-Earth@lternatives raster: %s' % title
        raster.institution = 'University of Twente, Netherlands'; raster.contact = 'o.mialyk@utwente.nl'
        
        # Dimensions
        param = raster.createDimension('param', par_num)
        lat = raster.createDimension('lat', rows30)
        lon = raster.createDimension('lon', cols30)
        
        # Variables
        lats = raster.createVariable('lat', 'f4', ('lat',))
        lons = raster.createVariable('lon', 'f4', ('lon',))
        
        i = 0
        for par in parameters:
            val = raster.createVariable(par, 'f4', ('lat','lon',), \
                                       zlib = True, complevel = 6, fill_value = -999)
            val[:,:] = values[i,:,:]; i += 1
        # Adding values
        incrm = 0.5
        lats[:] = np.arange(90-incrm/2, -rows30/4, -incrm)
        lats.long_name = 'lat'; lats.units = 'degrees_north'
        lons[:] = np.arange(-180+incrm/2, cols30/4, incrm)
        lons.long_name = 'lon'; lons.units = 'degrees_east'
        raster.close()