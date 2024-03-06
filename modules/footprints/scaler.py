# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 14:22:52 2022

@author: MialykO
"""

#%% Import packages
import numpy as np
import os
import pandas as pd
from  datetime import datetime
import modules.acea.acea_core as ac
import modules.footprints.footprints_core as fc
from osgeo import gdal

#%% Classes
class acea_ha_scaler:
    # Folder paths
    data_folder = ac.GetPath(['data','footprints'])
    
    def __init__(self, conf):
        self.conf = conf; _,_,self.res5 = ac.GetResolution(1)
        self.clock_start_year = int(self.conf.clock_start.split('/')[0])
        self.clock_end_year = int(self.conf.clock_end.split('/')[0])
        self.max_gs = (self.clock_end_year - self.clock_start_year)+1 - self.conf.spinup; self.first_year = self.clock_start_year + self.conf.spinup 
        self.text_years = f'{self.first_year}_{self.clock_end_year}'
        self.crop_code = conf.crop_fao
    
    def GenerateCropAreas(self, ir = True, rf = True):
        "Function to make SPAM2010 historical (1990-2019) using HID and HYDE 3.2"
        import scipy.interpolate
        
        # Define some general variables
        area_threshold = .1 # Mask cells with farm size less than 0.1 ha
        
        # Read physiscal (actual) areas in SPAM
        if rf:
            if self.conf.landuse == "spam2010":
                area_pa_rain5 = gdal.Open(ac.GetPath([self.data_folder, "spam", f'spam2010V2r0_global_A_{self.crop_code}_R.tif']))
                area_pa_rain5 = area_pa_rain5.GetRasterBand(1).ReadAsArray()
                area_pa_rain5 = np.ma.masked_where(area_pa_rain5 <= 0, area_pa_rain5)
            elif self.conf.landuse == "gaez2015":
                area_pa_rain5 = gdal.Open(ac.GetPath([self.data_folder, "gaez", f'GAEZAct2015_HarvArea_Rainfed_{self.crop_code}.tif']))
                area_pa_rain5 = area_pa_rain5.GetRasterBand(1).ReadAsArray()
                area_pa_rain5 = np.ma.masked_where(area_pa_rain5 <= 0, area_pa_rain5)
            else:
                raise Exception("Crop extent is not found") 
                
            area_pa_rain5[(area_pa_rain5 < area_threshold) & (area_pa_rain5>0)] = np.ma.masked
            
            # Get HYDE 3.2 physical cropland ([km2] for 1990, 2000:1:2015), Reference: https://doi.org/10.5194/essd-9-927-2017
            hyde_data_rf = ac.ReadNC(ac.GetPath([self.data_folder, 'cropland','HYDE_3_2_rf_1980_2016.nc']), "hyde")[1:,:,:]*100 # [km2] -> [ha] Get only years since 1990
            hyde_data_rf[(hyde_data_rf<area_threshold) & (hyde_data_rf>0)] = 0
        
        if ir:
            if self.conf.landuse == "spam2010":
                area_pa_ir5 = gdal.Open(ac.GetPath([self.data_folder, "spam", f'spam2010V2r0_global_A_{self.crop_code}_I.tif']))
                area_pa_ir5 = area_pa_ir5.GetRasterBand(1).ReadAsArray()
                area_pa_ir5 = np.ma.masked_where(area_pa_ir5 <= 0, area_pa_ir5)
            elif self.conf.landuse == "gaez2015":
                area_pa_ir5 = gdal.Open(ac.GetPath([self.data_folder, "gaez", f'GAEZAct2015_HarvArea_Irrigated_{self.crop_code}.tif']))
                area_pa_ir5 = area_pa_ir5.GetRasterBand(1).ReadAsArray()
                area_pa_ir5 = np.ma.masked_where(area_pa_ir5 <= 0, area_pa_ir5)
            else:
                raise Exception("Crop extent is not found") 
                
            area_pa_ir5[(area_pa_ir5 < area_threshold) & (area_pa_ir5>0)] = np.ma.masked

            # Get HID data (AEI [ha] for 1990, 1995, 2000, 2005), Reference: https://hess.copernicus.org/articles/19/1521/2015/
            hid_data_ir = ac.ReadNC(ac.GetPath([self.data_folder, 'cropland','HID_aei_ha_1985_2005.nc']), "hid_area")[1:,:,:] # Get only last 4 bands (1990-2005)
            hid_data_ir[np.isnan(hid_data_ir)] = np.ma.masked
            hid_data_ir[(hid_data_ir < area_threshold) & (hid_data_ir > 0)] = 0
        
            # Get HYDE 3.2 physical cropland ([km2] for 1990, 2000:1:2016), Reference: https://doi.org/10.5194/essd-9-927-2017
            hyde_data_ir = ac.ReadNC(ac.GetPath([self.data_folder, 'cropland','HYDE_3_2_aei_1980_2016.nc']), "hyde")[1:,:,:]*100 # [km2] -> [ha]
            hyde_data_ir[(hyde_data_ir<area_threshold) & (hyde_data_ir>0)] = 0
            
        # Make SPAM2010 historical
        if ir: #irrigated
            print('Interpolating irrigated areas')
            cell_count = 0
            hist_area = np.zeros((self.max_gs,self.res5[0],self.res5[1])); hist_area = np.ma.masked_where(hist_area == 0, hist_area)
            area_pa_ir5_list = np.array([*np.where(area_pa_ir5.mask==False)]).transpose()
            
            for cell in area_pa_ir5_list:
                if cell_count % 5000 == 0: print(f'Cells read: {cell_count} out of {len(area_pa_ir5_list)} ({cell_count/len(area_pa_ir5_list)*100:.2f} %)')
                row5 = int(cell[0]); col5 = int(cell[1])
                
                temp_final = np.ones(self.max_gs) * area_pa_ir5[row5, col5] # Set SPAM2010 by default [ha]
                
                if any(hid_data_ir[:, row5, col5] > 0): # Fit SPAM to cropland dynamics
                    
                    temp_hyde = np.zeros(self.max_gs); temp_hid = np.zeros(self.max_gs)
                    
                    # Step 4.a.1 Interpolate HYDE and HID
                    area_interp = scipy.interpolate.interp1d([1990, *range(2000,2017)], hyde_data_ir[:, row5, col5]); y = 0
                    for x in range(1990, 2017): temp_hyde[y] = area_interp(x); y += 1 # Iterpolate years
                    temp_hyde[-3:] = temp_hyde[-4] # Copy 2016 value for 2017-2019
                    area_interp = scipy.interpolate.interp1d([1990, 1995, 2000, 2005], hid_data_ir[:, row5, col5]); y = 0
                    for x in range(1990, 2006): temp_hid[y] = area_interp(x); y += 1 # Iterpolate years
                    
                     # Step 4.a.2 Extrapolate HID by HYDE
                    if any(temp_hyde[-14:]!=0) and temp_hyde[-15]!=0: # Check if HYDE values for 2006-2019 exist
                        temp_hid[-14:] = temp_hyde[-14:] * temp_hid[-15] / temp_hyde[-15]
                    else: # Otherwise copy the 2005 HID value for 2006-2019
                        temp_hid[-14:] = temp_hid[-15]
                    
                    # Step 4.a.3  Scale according to HID
                    temp_vals = temp_hid[18:23] # 2008-2012
                    if any(temp_vals > 0):
                        avg_hid2010 = np.mean(temp_vals[np.nonzero(temp_vals)]) # Check the avg year 2010 for HID (2008-2012)
                    else:
                        avg_hid2010 = np.mean(temp_hid[np.nonzero(temp_hid)][-5:]) # If none, then take the last 5 data points
                    
                    # Scale data to represent the trends in historical irrigated cropland extend
                    temp_final[:] = np.around(np.min(np.vstack((temp_final * (temp_hid/avg_hid2010), temp_hid)),axis=0),4)
                    
                    # if any(temp_final == 0): temp_final[temp_final==0] = area_pa_ir5[row5, col5] # Missing years in HYDE set as SPAM
                    temp_final[(temp_final<area_threshold) & (temp_final>0)] = 0
                    
                hist_area[:, row5, col5] = temp_final*1; cell_count += 1
            
            # Save raster
            raster_path = ac.GetPath([self.data_folder, "phys_area", f'{self.conf.crop_fao}_area_ir_phys_global_annual_{self.text_years}.nc'])
            title = f'Crop: {self.conf.crop_name}, Basis: SPAM2010 irrigated'
            val_name = 'harvested_area'; val_unit = 'ha y-1'
            ac.CreateHistRaster5(hist_area, raster_path, val_name, val_unit, title, self.first_year, self.max_gs)
        
        if rf: #Rainfed
            print('Interpolating rainfed areas')
            cell_count = 0
            hist_area = np.zeros((self.max_gs,self.res5[0],self.res5[1])); hist_area = np.ma.masked_where(hist_area == 0, hist_area)
            area_pa_rain5_list = np.array([*np.where(area_pa_rain5.mask==False)]).transpose()
            
            for cell in area_pa_rain5_list:
                if cell_count % 5000 == 0:print(f'Cells read: {cell_count} out of {len(area_pa_rain5_list)} ({cell_count/len(area_pa_rain5_list)*100:.2f} %)')
                row5 = int(cell[0]); col5 = int(cell[1])
                
                temp_final = np.ones(self.max_gs) * area_pa_rain5[row5, col5] # Set SPAM2010 by default [ha]
                
                if any(hyde_data_rf[:, row5, col5] > 0): # Fit SPAM to cropland dynamics
                    temp_hyde = np.zeros(self.max_gs)
                    
                    # Step 4.b.1 Interpolate HYDE
                    area_interp = scipy.interpolate.interp1d([1990, *range(2000,2017)], hyde_data_rf[:, row5, col5]); y = 0
                    for x in range(1990, 2017, 1): temp_hyde[y] = area_interp(x); y += 1 # Iterpolate years
                    temp_hyde[-3:] = temp_hyde[-4] # Copy 2016 value for 2017-2019
                    
                    # Step 4.b.2 Scale according to HYDE 3.2
                    temp_vals = temp_hyde[18:23] # 2008-2012
                    if any(temp_vals > 0):
                        avg_hyde2010 = np.mean(temp_vals[np.nonzero(temp_vals)]) # Check the non zero avg year 2010 for HYDE (2008-2012)
                    else:
                        avg_hyde2010 = np.mean(temp_hyde[np.nonzero(temp_hyde)][-5:]) # If none, then take the last 5 data points
                    
                    temp_final[:] = np.around(np.min(np.vstack((temp_final * (temp_hyde/avg_hyde2010), temp_hyde)),axis=0),4) # Scale data to represent the trends in historical rainfed cropland extend
                    
                    # if any(temp_final == 0): temp_final[temp_final==0] = area_pa_rain5[row5, col5] # Missing years in HYDE set as SPAM
                    temp_final[(temp_final<area_threshold) & (temp_final>0)] = 0
                    
                hist_area[:, row5, col5] = temp_final*1; cell_count += 1
            
            # Save raster
            raster_path = ac.GetPath([self.data_folder, "phys_area", f'{self.conf.crop_fao}_area_rf_phys_global_annual_{self.text_years}.nc'])
            title = f'Crop: {self.conf.crop_name}, Basis: SPAM2010 rainfed'
            val_name = 'harvested_area'; val_unit = 'ha y-1'
            ac.CreateHistRaster5(hist_area, raster_path, val_name, val_unit, title, self.first_year, self.max_gs)
    
    def ScaleAreasToFAOSTAT(self):
        "Scale physical areas to reported harvested areas in FAOSTAT (revised 18/01/2023)"
        print(f'Scaling areas of {self.conf.crop_name}')
        
        # Create report
        out_report = ac.GetPath(['outputs', 'main','scaling_logs','harvested_area_scaling', f'{self.conf.crop_fao}_NationalScalingReport_{self.text_years}_spam2010.txt'])
        if os.path.exists(out_report): os.remove(out_report)
        ac.printLog(out_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

        # Get inputs
        faostat_ha = pd.read_csv(ac.GetPath([self.data_folder, 'faostat', 'harvested_area_1990_2019.csv']), usecols = ['Area Code','Crop Code','Year','Harvested area(ha)'])
        faostat_ha = faostat_ha.loc[faostat_ha['Crop Code']==self.conf.crop_fao]
        hist_phys_area_rf = ac.ReadNC(ac.GetPath([self.data_folder, "phys_area", f'{self.conf.crop_fao}_area_rf_phys_global_annual_{self.text_years}.nc']),'harvested_area')
        hist_phys_area_ir = ac.ReadNC(ac.GetPath([self.data_folder, "phys_area", f'{self.conf.crop_fao}_area_ir_phys_global_annual_{self.text_years}.nc']),'harvested_area')
        scaled_cells = np.zeros((hist_phys_area_rf.shape)) # To indicate which cells have data
        country_data = fc.GetListOfHistCountries(); country_codes = np.unique(faostat_ha['Area Code'].tolist())
        
        # Start scaling
        country_num = 1
        for country in country_codes: # Iterate countries
            if country == 351: country_num += 1; continue # Skip China, because there are separate codes for mainland and Taiwan
            c_name = country_data.loc[country_data['Country Code'] == country]['Country'].values[0].replace("\n", "")
            ac.printLog(out_report, f'Country: {c_name} ({country}), {country_num} out of {len(country_codes)}')
            ac.printLog(out_report, 'Year \t FAO area [ha] \t Extrapolated Area [ha] \t Scaling factor')
            fao_ha = faostat_ha[faostat_ha['Area Code'] == country]['Harvested area(ha)'].tolist() # Reported national production for cur year
            fao_year = faostat_ha[faostat_ha['Area Code'] == country]['Year'].tolist() # Reported national production for cur year
            c_cells = fc.GetHistCountryCells5(country).values # 5 arcmin grid cells for country X
            
            c_raster =  np.ones((self.res5[0],self.res5[1])); c_raster = np.ma.array(c_raster, mask=(c_raster==1))
            for cell in c_cells: c_raster[cell[1], cell[2]] = 1
            
            # Scale the harvest areas to fit FAOSTAT
            for gs in range(0, self.max_gs, 1):
                cur_year = self.first_year + gs # Current year
                rf_area = hist_phys_area_rf[gs,:,:]; ir_area = hist_phys_area_ir[gs,:,:]
                _scaled_cells = scaled_cells[gs,:,:]
                c_rf_area = np.ma.filled(rf_area[c_raster==1],0); c_ir_area = np.ma.filled(ir_area[c_raster==1],0) # Get only the areas of the country
                # max_area = np.max(cell_area[c_raster==1]) # Max area size for the current country
                # c_tot_area = np.sum(np.min(np.vstack((c_rf_area + c_ir_area, np.ones(len(c_rf_area))*max_area)),axis=0)) # Limit max harvest area and take sum
                c_tot_area = np.sum(c_rf_area + c_ir_area)
                
                if (c_tot_area > 10) and (cur_year in fao_year): # Check if FAO reported data for this year, cut minor producers
                    factor = fao_ha[fao_year.index(cur_year)]/c_tot_area;_scaled_cells[c_raster==1] = 1
                    ac.printLog(out_report, f'{cur_year}\t{int(fao_ha[fao_year.index(cur_year)])}\t{int(c_tot_area)}\t{factor:.3f}')
                    if np.isnan(factor): factor = np.ma.masked
                    rf_area[c_raster==1] = rf_area[c_raster==1]*factor; ir_area[c_raster==1] = ir_area[c_raster==1]*factor
                        # Python overwrites harv_area_rf and harv_area_irr automatically (rf_area and ir_area are reference variables)
            country_num += 1
            
        hist_phys_area_rf[scaled_cells==0] = np.ma.masked; hist_phys_area_ir[scaled_cells==0] = np.ma.masked # Set cells with no FAO data to nan
        
        # Save raster
        raster_folder_path = ac.GetPath([self.data_folder,'harvested_areas', self.conf.crop_fao])
        if not os.path.exists(raster_folder_path): os.makedirs(raster_folder_path)
        
        raster_path = ac.GetPath([raster_folder_path,  f'{self.conf.crop_fao}_area_rf_harv_global_annual_{self.text_years}.nc'])
        title = f'Crop: {self.conf.crop_name}, Basis: SPAM2010 rainfed'
        val_name = 'harvested_area'; val_unit = 'ha y-1'
        ac.CreateHistRaster5(hist_phys_area_rf, raster_path, val_name, val_unit, title, self.first_year, self.max_gs)
        
        raster_path = ac.GetPath([raster_folder_path, f'{self.conf.crop_fao}_area_ir_harv_global_annual_{self.text_years}.nc'])
        title = f'Crop: {self.conf.crop_name}, Basis: SPAM2010 irrigated'
        ac.CreateHistRaster5(hist_phys_area_ir, raster_path, val_name, val_unit, title, self.first_year, self.max_gs)
        
class acea_prod_scaler:
    "To downscale the results and scale the production to FAOSTAT"
    
    def __init__(self, conf):
        self.conf = conf
        self.data_folder = ac.GetPath(['data','footprints'])
        self.raster_folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'rasters'])
        self.gridcells30 = ac.GetAllCells30() # Get gridcells in 30 arc min (id30,rowy30,colx30)
        
        self.clock_start_year = int(self.conf.clock_start.split('/')[0])
        self.clock_end_year = int(self.conf.clock_end.split('/')[0])
        self.max_gs = (self.clock_end_year - self.clock_start_year)+1 - self.conf.spinup; self.first_year = self.clock_start_year + self.conf.spinup
        _,_,self.res5 = ac.GetResolution(1)
        self.text_years = f'{self.first_year}_{self.clock_end_year}'
        
    def ScaleProductionToFAOSTAT(self):
        "Scale production to national FAOSTAT data (revised 18/01/2023)"
        
        print(f'Scaling production of {self.conf.crop_name}')
        out_report = ac.GetPath(['outputs', 'main','scaling_logs','production_scaling', f'{self.conf.crop_fao}_NationalProdScalingReport_{self.text_years}.txt'])
        if os.path.exists(out_report): os.remove(out_report)
        ac.printLog(out_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        
        # 1. Get simulated production
        print('1. Getting simulated production')
        harv_area_rf = fc.GetHistHarvestedAreas5(self.conf.crop_fao, False); harv_area_ir = fc.GetHistHarvestedAreas5(self.conf.crop_fao, True) 
        
        for _,_, files in os.walk(self.raster_folder_path):
            for f in files:
                if 'yield' in f and '5arc'in f:
                    if 'sc1' in f: sim_yields_sc1, _,_,_ = ac.GetOutputsOfScenario(self.conf, 1, False, False) # Rainfed no GW
                    if 'sc2' in f: sim_yields_sc2, _,_,_ = ac.GetOutputsOfScenario(self.conf, 2, False, False) # Rainfed with GW
                    if 'sc3' in f: sim_yields_sc3, _,_,_ = ac.GetOutputsOfScenario(self.conf, 3, False, False) # Irrigated furrow
                    if 'sc4' in f: sim_yields_sc4, _,_,_ = ac.GetOutputsOfScenario(self.conf, 4, False, False) # Irrigated sprinkler
                    if 'sc5' in f: sim_yields_sc5, _,_,_ = ac.GetOutputsOfScenario(self.conf, 5, False, False) # Irrigated drip
                    if 'sc6' in f: sim_yields_sc6, _,_,_ = ac.GetOutputsOfScenario(self.conf, 6, False, False) # Irrigated flood with gw
                        
        # 2. Get FAOSTAT production
        print('2. Getting FAOSTAT production')
        faostat_prod = pd.read_csv(ac.GetPath([self.data_folder,'faostat','production_1990_2019.csv']), usecols = ['Area Code','Crop Code','Year','Production(t)'])
        faostat_prod = faostat_prod.loc[faostat_prod['Crop Code']==self.conf.crop_fao]
        country_codes = np.unique(faostat_prod['Area Code'].tolist())
        country_data = fc.GetListOfHistCountries()
        
        #3. Scale production
        print('3. Scaling production')
        scaling_final = np.ones((self.max_gs,self.res5[0],self.res5[1]))*np.ma.masked
        
        if sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
            ir_syst_data_sur, ir_syst_data_spr,ir_syst_data_drp = fc.GetIrrSystems5() # Get irrigation fractions
            
        country_num = 1
        for country in country_codes: # Iterate countries
            if country == 351: country_num += 1; continue # Skip China, because there are separate codes for mainland and Taiwan
            c_name = country_data.loc[country_data['Country Code'] == country]['Country'].values[0].replace("\n", "")
            
            fao_year = faostat_prod[faostat_prod['Area Code'] == country]['Year'].tolist() # Reported years
            fao_prod = faostat_prod[faostat_prod['Area Code'] == country]['Production(t)'].tolist() # Reported national production for cur year
            
            c_cells = fc.GetHistCountryCells5(country).values # 5 arcmin grid cells for country X
            c_raster =  np.zeros((self.res5[0],self.res5[1])); c_raster = np.ma.array(c_raster, mask=(c_raster==0))
            for cell in c_cells: c_raster[cell[1], cell[2]] = 1
            
            ac.printLog(out_report, f'Country: {c_name} ({country}), {country_num} out of {len(country_codes)}')
            ac.printLog(out_report, 'Year \t FAO production [t] \t Simulated production [t] \t Scaling factor')
            for gs in range(self.max_gs): # Iterate year
                cur_year = self.first_year + gs # Current year
                if (cur_year in fao_year) and fao_prod[fao_year.index(cur_year)] > 0:
                    scaling_final_ref = scaling_final[gs,:,:]
                    prod_raster_rf = np.ones((self.res5[0],self.res5[1]))* np.ma.masked
                    prod_raster_ir =  prod_raster_rf*1; _sim_yields_ir = prod_raster_rf*1
                    
                    if min(self.conf.scenarios) < 3:
                        _harv_area_rf = harv_area_rf[gs,:,:]
                        _sim_yields_sc1 = sim_yields_sc1[gs,:,:]; _sim_yields_sc2 = sim_yields_sc2[gs,:,:]
                        
                        prod_raster_rf[c_raster==1] = np.around(_harv_area_rf[c_raster==1] * \
                                        (np.ma.filled(_sim_yields_sc1[c_raster==1],0) + np.ma.filled(_sim_yields_sc2[c_raster==1],0)),2) # [ton]
                        
                        sum_acea_rf = np.ma.sum(prod_raster_rf[c_raster==1])
                        if np.ma.getmask(sum_acea_rf): sum_acea_rf = 0
                        
                    if max(self.conf.scenarios) > 2:
                        _harv_area_ir = harv_area_ir[gs,:,:]; sum_acea_ir = 0
                        
                        if not np.ma.getmask(np.ma.sum(_harv_area_ir[c_raster==1])):
                            if sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
                                _sim_yields_sc3 = sim_yields_sc3[gs,:,:]
                                _sim_yields_sc4 = sim_yields_sc4[gs,:,:]
                                if 5 not in self.conf.scenarios: # Only strategies 3 and 4
                                    ir_syst_data_sur[(ir_syst_data_sur==0) & (ir_syst_data_spr==0)] = 1 # Assume surface for only drip cells
                                    sys_sum = ir_syst_data_sur + ir_syst_data_spr; ir_syst_data_sur = ir_syst_data_sur/sys_sum; ir_syst_data_spr = ir_syst_data_spr/sys_sum
                                    _sim_yields_ir[c_raster==1] = np.ma.filled(_sim_yields_sc3[c_raster==1] * ir_syst_data_sur[c_raster==1],0) + \
                                                      np.ma.filled(_sim_yields_sc4[c_raster==1] * ir_syst_data_spr[c_raster==1],0)
                                else: # All irrigation strategies
                                    _sim_yields_sc5 = sim_yields_sc5[gs,:,:]
                                    _sim_yields_ir[c_raster==1] = np.ma.filled(_sim_yields_sc3[c_raster==1] * ir_syst_data_sur[c_raster==1],0) + \
                                                      np.ma.filled(_sim_yields_sc4[c_raster==1] * ir_syst_data_spr[c_raster==1],0) + \
                                                      np.ma.filled(_sim_yields_sc5[c_raster==1] * ir_syst_data_drp[c_raster==1],0)
                            else: # only one strategy is considered
                                if 3 in self.conf.scenarios:
                                    _sim_yields_sc3 = sim_yields_sc3[gs,:,:]
                                    _sim_yields_ir[c_raster==1] = np.ma.filled(_sim_yields_sc3[c_raster==1],0)
                                if 6 in self.conf.scenarios:
                                    _sim_yields_sc6 = sim_yields_sc6[gs,:,:]
                                    _sim_yields_ir[c_raster==1] = np.ma.filled(_sim_yields_sc6[c_raster==1],0)
                                
                            prod_raster_ir[c_raster==1] = np.around(_harv_area_ir[c_raster==1] * _sim_yields_ir[c_raster==1],2) # [ton]
                            sum_acea_ir = np.ma.sum(prod_raster_ir[c_raster==1])
                            if np.ma.getmask(sum_acea_ir): sum_acea_ir = 0
                        
                    acea_prod = np.around(sum_acea_rf + sum_acea_ir, 2)
                    if acea_prod < 10: continue
                    else: scaling_factor = np.around(fao_prod[fao_year.index(cur_year)]/acea_prod,5)
                        
                    scaling_final_ref[c_raster==1] = scaling_factor
                    ac.printLog(out_report,f'{cur_year}\t{int(fao_prod[fao_year.index(cur_year)])}\t{int(acea_prod)}\t{scaling_factor:.3f}')
        
            country_num += 1
            
        # Save raster
        print('4. Saving annual scaling factors')
        raster_name = f'scaling_factors_5arc_{self.conf.crop_fao}_global_annual_{self.text_years}.nc'
        raster_path = ac.GetPath([self.data_folder,'scaling_factors', raster_name])
        val_name = f'scaling_factor-{self.conf.crop_name_short}'
        title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: yield scaling factor'
        
        ac.CreateHistRaster5(scaling_final, raster_path, val_name, '-', title, self.first_year, self.max_gs)
        
        print('5. Calculate annual scaling factors with moving average')
        self.ScalingFactorsMovingAvg()
        
        print('6. Create production raster')
        self.CreateProductionRaster()
        
    def ScalingFactorsMovingAvg(self, mid_weight=0.34):
        "Calculate moving average with 3 year window (+1y before and +1y after)"
         # mid_weight = Weight of the middle year [0 to 1]
         
        raster_name = f'scaling_factors_5arc_{self.conf.crop_fao}_global_annual_{self.text_years}.nc'
        raster_path = ac.GetPath([self.data_folder,'scaling_factors', raster_name])
        scaling_factors = ac.ReadNC(raster_path, f'scaling_factor-{self.conf.crop_name_short}') # Get scaling factors for yields
         
        harv_area_rf = fc.GetHistHarvestedAreas5(self.conf.crop_fao, False); harv_area_ir = fc.GetHistHarvestedAreas5(self.conf.crop_fao, True)
        
        scaling_factors[(harv_area_rf.mask==True) & (harv_area_ir.mask==True) & (harv_area_rf!=0) & (harv_area_ir!=0)] = np.ma.masked
        factor_weights = [(1-mid_weight)/2, mid_weight, (1-mid_weight)/2]
        percent_complete = 0
        for row in range(0, self.res5[0]):
            if percent_complete != int(row/self.res5[0]*100):
                percent_complete = int(row/self.res5[0]*100)
                print(f'{percent_complete} % is done')
                
            for col in range(0, self.res5[1]):
                if any(np.ma.getmask(scaling_factors[:, row, col]) == False):
                    temp_factors = scaling_factors[:, row, col]*1; final_factors = np.zeros(temp_factors.shape[0])*np.ma.masked
                    final_factors[0] = temp_factors[0]; final_factors[-1] = temp_factors[-1]
                    for i in range(1,(temp_factors.shape[0]-1)):
                        if all(np.isnan(temp_factors[i-1:i+2])==False):
                            final_factors[i] = np.average(temp_factors[i-1:i+2], weights = factor_weights)
                        elif all(np.isnan(temp_factors[i-1:i+1])==False):
                            final_factors[i] = np.average(temp_factors[i-1:i+1], weights = factor_weights[:2])
                        elif all(np.isnan(temp_factors[i:i+2])==False):
                            final_factors[i] = np.average(temp_factors[i:i+2], weights = factor_weights[-2:])
                    
                    scaling_factors[:, row, col] = np.around(final_factors, 4)
                    
        # Save raster
        scaling_factors[scaling_factors==0] = np.ma.masked
        
        raster_name = f'scaling_factors_5arc_{self.conf.crop_fao}_moving_avg_global_annual_{self.first_year}_{self.clock_end_year}.nc'
        raster_path = ac.GetPath([self.data_folder,'scaling_factors', raster_name])
        val_name = f'scaling_factor-{self.conf.crop_name_short}'
        title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: yield scaling factor'
        
        ac.CreateHistRaster5(scaling_factors, raster_path, val_name, '-', title, self.first_year, self.max_gs)
        
        return scaling_factors
    
    def CreateProductionRaster(self):
        "Create a raster with ranfed/irrigated annual production"
       
        # 1. Get rainfed yields[fresh tonne ha-1]
        _sim_yields_sc1, _,_,_ = ac.GetOutputsOfScenario(self.conf, 1, True, False) # Rainfed no GW
        _sim_yields_sc2, _,_,_ = ac.GetOutputsOfScenario(self.conf, 2, True, False) # Rainfed with GW
        sim_yields_rf = np.ma.filled(_sim_yields_sc1, 0) + np.ma.filled(_sim_yields_sc2, 0)
        
        # 2. Get irrigated yields[fresh t ha-1]
        if sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
            ir_syst_data_sur, ir_syst_data_spr, ir_syst_data_drp = fc.GetIrrSystems5() # Get irrigation fractions
            _sim_yields_sc3, _,_,_ = ac.GetOutputsOfScenario(self.conf, 3, True, False) # Irrigated surface
            _sim_yields_sc4, _,_,_ = ac.GetOutputsOfScenario(self.conf, 4, True, False) # Irrigated sprinkler
            if 5 not in self.conf.scenarios: # Only strategies 3 and 4
                ir_syst_data_sur[(ir_syst_data_sur==0) & (ir_syst_data_spr==0)] = 1 # Assume surface for only drip cells
                sys_sum = ir_syst_data_sur + ir_syst_data_spr; ir_syst_data_sur = ir_syst_data_sur/sys_sum; ir_syst_data_spr = ir_syst_data_spr/sys_sum
                sim_yields_ir = np.ma.filled(_sim_yields_sc3, 0) * ir_syst_data_sur + \
                                np.ma.filled(_sim_yields_sc4, 0) * ir_syst_data_spr
            else:
                _sim_yields_sc5, _,_,_ = ac.GetOutputsOfScenario(self.conf, 5, True, False) # Irrigated drip
                sim_yields_ir = np.ma.filled(_sim_yields_sc3, 0) * ir_syst_data_sur + \
                                np.ma.filled(_sim_yields_sc4, 0) * ir_syst_data_spr + \
                                np.ma.filled(_sim_yields_sc5, 0) * ir_syst_data_drp
        else:
            if 3 in self.conf.scenarios:
                sim_yields_ir, _,_,_ = np.ma.filled(ac.GetOutputsOfScenario(self.conf, 3, True, False),0) # Irrigated surface
            if 6 in self.conf.scenarios: # For rice only
                sim_yields_ir, _,_,_ = np.ma.filled(ac.GetOutputsOfScenario(self.conf, 6, True, False),0) # Irrigated surface

        # 3. Get scaled harvest areas [ha]
        harv_area_rf = fc.GetHistHarvestedAreas5(self.conf.crop_fao, False); harv_area_ir = fc.GetHistHarvestedAreas5(self.conf.crop_fao, True)

        # 4. Calculate production
        prod_rf = harv_area_rf * sim_yields_rf; prod_rf[prod_rf==0] = np.ma.masked
        prod_ir = harv_area_ir * sim_yields_ir; prod_ir[prod_ir==0] = np.ma.masked
        
        # 5. Save rasters
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.text_years}.nc'])
        val_name = f'production-{self.conf.crop_name_short}'
        title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}'
        ac.CreateHistRaster5(prod_rf, _raster_path, val_name, 'fresh t y-1', title, self.first_year, self.max_gs)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.text_years}.nc'])
        ac.CreateHistRaster5(prod_ir, _raster_path, val_name, 'fresh t y-1', title, self.first_year, self.max_gs)