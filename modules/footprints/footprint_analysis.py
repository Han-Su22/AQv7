# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 09:15:56 2022

@author: MialykO
"""
#%% Import packages
import numpy as np, os, importlib, pandas as pd
from  datetime import datetime
import modules.acea.acea_core as ac, modules.footprints.footprints_core as fc

class wf_analyser:
    "To calculate and analyse crop WFs"
    
    def __init__(self, conf):
        self.conf = conf
        self.data_folder = ac.GetPath(['data','footprints'])
        self.raster_folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'rasters'])
        self.clock_start_year = int(self.conf.clock_start.split('/')[0])
        self.clock_end_year = int(self.conf.clock_end.split('/')[0])
        self.max_gs = (self.clock_end_year - self.clock_start_year)+1 - self.conf.spinup; self.first_year = self.clock_start_year + self.conf.spinup
        _,_,self.res5 = ac.GetResolution(1)
        self.text_years = f'{self.first_year}_{self.clock_end_year}'
        
        print(f"============Analysis of {self.conf.crop_name} WFs============")
        
    def CalculateIrReq(self):
        "Calculate potential net irrigation"
        print("Calculating potential net irrigation")
        if max(self.conf.scenarios) > 2: # Check for Irrigated scenarios
            ir_req = np.zeros((self.max_gs,self.res5[0],self.res5[1]))
            if sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
                ir_syst_data_sur, ir_syst_data_spr, ir_syst_data_drp = fc.GetIrrSystems5() # Get irrigation fractions
                print('Scenario 3-5...')
                sim_irreq_sc3 = np.ma.filled(ac.ReadNC(ac.GetPath([self.raster_folder_path, f'acea_5arc_sc3_{self.conf.crop_name_short}_pirnreq_global_annual_{self.text_years}.nc']), f'pirnreq-{self.conf.crop_name_short}'),0)
                sim_irreq_sc4 = np.ma.filled(ac.ReadNC(ac.GetPath([self.raster_folder_path, f'acea_5arc_sc4_{self.conf.crop_name_short}_pirnreq_global_annual_{self.text_years}.nc']), f'pirnreq-{self.conf.crop_name_short}'),0)
                
                if 5 not in self.conf.scenarios: # Only strategies 3 and 4
                    ir_syst_data_sur[(ir_syst_data_sur==0) & (ir_syst_data_spr==0)] = 1 # Assume surface for only drip cells
                    sys_sum = ir_syst_data_sur + ir_syst_data_spr; ir_syst_data_sur = ir_syst_data_sur/sys_sum; ir_syst_data_spr = ir_syst_data_spr/sys_sum
                    ir_req = sim_irreq_sc3*ir_syst_data_sur + sim_irreq_sc4*ir_syst_data_spr
                else:
                    sim_irreq_sc5 = np.ma.filled(ac.ReadNC(ac.GetPath([self.raster_folder_path, f'acea_5arc_sc5_{self.conf.crop_name_short}_pirnreq_global_annual_{self.text_years}.nc']), f'pirnreq-{self.conf.crop_name_short}'),0)
                    ir_req = sim_irreq_sc3*ir_syst_data_sur + sim_irreq_sc4*ir_syst_data_spr + sim_irreq_sc5*ir_syst_data_drp
                    
            else: # only one strategy is considered
                if 3 in self.conf.scenarios:
                    print('Scenario 3...')
                    ir_req = np.ma.filled(ac.ReadNC(ac.GetPath([self.raster_folder_path, f'acea_5arc_sc3_{self.conf.crop_name_short}_pirnreq_global_annual_{self.text_years}.nc']), f'pirnreq-{self.conf.crop_name_short}'),0)
                if 6 in self.conf.scenarios:
                    print('Scenario 6...')
                    ir_req = np.ma.filled(ac.ReadNC(ac.GetPath([self.raster_folder_path, f'acea_5arc_sc6_{self.conf.crop_name_short}_pirnreq_global_annual_{self.text_years}.nc']), f'pirnreq-{self.conf.crop_name_short}'),0)
                
            ir_req = np.ma.masked_where(ir_req == 0, ir_req)
            
            raster_title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: potential net irrigation'
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_pirnreq_global_annual_{self.text_years}.nc'])
            val_name = f'pirnreq-{self.conf.crop_name_short}'
            ac.CreateHistRaster5(ir_req, _raster_path, val_name, 'mm y-1', raster_title, self.first_year, self.max_gs)
        
    def CalculateWF(self):
        "Calculate green and blue WFs"
        print("Calculating green and blue WFs")
        
        raster_title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit water footprint'
        
        # Calculating footprints[m3 t-1]
        print('Calculating footprints...')
        wf_green_rf, wf_green_ir, wf_blue_rf, wf_blue_ir= (np.zeros((self.max_gs,self.res5[0],self.res5[1])) for _ in range(4))
        
        # Rainfed
        if min(self.conf.scenarios) < 3: # Check for rainfed scenarios
            if 1 in self.conf.scenarios:
                print('Scenario 1...')
                sim_yields_sc1, sim_cwu_green_sc1, _, _ = ac.GetOutputsOfScenario(self.conf, 1) # Rainfed no GW, True = scaled
                sim_yields_sc1 = np.ma.masked_where(sim_yields_sc1 == 0, sim_yields_sc1) # Mask 0 values  to avoid division by 0
                if len(sim_yields_sc1[sim_yields_sc1.mask==False])>0:
                    _limit = np.percentile(sim_yields_sc1[sim_yields_sc1.mask==False],1) # Filter out extremely low yield data points
                    sim_yields_sc1[sim_yields_sc1 < _limit] = _limit
                
                wf_green_rf = wf_green_rf + np.ma.filled((sim_cwu_green_sc1 / sim_yields_sc1 * 10), 0) # Only green WF in scenario 1
                
            if 2 in self.conf.scenarios:
                print('Scenario 2...')
                sim_yields_sc2, sim_cwu_green_sc2, _, sim_cwu_blue_cr_sc2 = ac.GetOutputsOfScenario(self.conf, 2) # Rainfed with GW
                sim_yields_sc2 = np.ma.masked_where(sim_yields_sc2 == 0, sim_yields_sc2)
                if len(sim_yields_sc2[sim_yields_sc2.mask==False])>0:
                    _limit = np.percentile(sim_yields_sc2[sim_yields_sc2.mask==False],1) # Filter out extremely low yield data points
                    sim_yields_sc2[sim_yields_sc2 < _limit] = _limit
                
                wf_green_rf = wf_green_rf + np.ma.filled((sim_cwu_green_sc2 / sim_yields_sc2 * 10), 0) # Green WF
                wf_blue_rf = wf_blue_rf + np.ma.filled((sim_cwu_blue_cr_sc2 / sim_yields_sc2 * 10), 0) # Blue from capillary rise
            
            wf_green_rf = np.ma.masked_where(wf_green_rf==0, wf_green_rf); wf_blue_rf = np.ma.masked_where(wf_blue_rf==0, wf_blue_rf)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_rf_global_annual_{self.text_years}.nc'])
            val_name = f'wf_green_rf-{self.conf.crop_name_short}'
            ac.CreateHistRaster5(wf_green_rf, _raster_path, val_name, 'm3 t-1 y-1', raster_title, self.first_year, self.max_gs)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_rf_global_annual_{self.text_years}.nc'])
            val_name = f'wf_blue_rf-{self.conf.crop_name_short}'
            ac.CreateHistRaster5(wf_blue_rf, _raster_path, val_name, 'm3 t-1 y-1', raster_title, self.first_year, self.max_gs)
        
        # Irrigated
        if max(self.conf.scenarios) > 2: # Check for Irrigated scenarios
            if sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
                ir_syst_data_sur, ir_syst_data_spr, ir_syst_data_drp = fc.GetIrrSystems5() # Get irrigation fractions
                print('Scenario 3-5...')
                sim_yields_sc3, sim_cwu_green_sc3, sim_cwu_blue_ir_sc3, _ = ac.GetOutputsOfScenario(self.conf, 3) # Irrigated surface
                sim_yields_sc3 = np.ma.masked_where(sim_yields_sc3 == 0, sim_yields_sc3)
                sim_yields_sc4, sim_cwu_green_sc4, sim_cwu_blue_ir_sc4, _ = ac.GetOutputsOfScenario(self.conf, 4) # Irrigated sprinkler
                sim_yields_sc4 = np.ma.masked_where(sim_yields_sc4 == 0, sim_yields_sc4)
                
                if len(sim_yields_sc3[sim_yields_sc3.mask==False])>0:
                    _limit = np.percentile(sim_yields_sc3[sim_yields_sc3.mask==False],1) # Filter out extremely low yield data points
                    sim_yields_sc3[sim_yields_sc3 < _limit] = _limit
                    sim_yields_sc4[sim_yields_sc4 < _limit] = _limit
                
                if 5 not in self.conf.scenarios: # Only strategies 3 and 4
                    ir_syst_data_sur[(ir_syst_data_sur==0) & (ir_syst_data_spr==0)] = 1 # Assume surface for only drip cells
                    sys_sum = ir_syst_data_sur + ir_syst_data_spr; ir_syst_data_sur = ir_syst_data_sur/sys_sum; ir_syst_data_spr = ir_syst_data_spr/sys_sum
                    wf_green_ir = wf_green_ir + np.ma.filled((sim_cwu_green_sc3 / sim_yields_sc3 * ir_syst_data_sur * 10), 0) + \
                                                np.ma.filled((sim_cwu_green_sc4 / sim_yields_sc4 * ir_syst_data_spr * 10), 0) # Green WF
                                
                    wf_blue_ir = wf_blue_ir   + np.ma.filled((sim_cwu_blue_ir_sc3 / sim_yields_sc3 * ir_syst_data_sur * 10), 0) + \
                                                np.ma.filled((sim_cwu_blue_ir_sc4 / sim_yields_sc4 * ir_syst_data_spr * 10), 0) # Blue from irrigation
                else:
                    sim_yields_sc5, sim_cwu_green_sc5, sim_cwu_blue_ir_sc5, _ = ac.GetOutputsOfScenario(self.conf, 5) # Irrigated drip
                    sim_yields_sc5 = np.ma.masked_where(sim_yields_sc5 == 0, sim_yields_sc5)
                    if len(sim_yields_sc5[sim_yields_sc5.mask==False])>0:
                        _limit = np.percentile(sim_yields_sc5[sim_yields_sc5.mask==False],1) # Filter out extremely low yield data points
                        sim_yields_sc5[sim_yields_sc5 < _limit] = _limit
                
                    wf_green_ir = wf_green_ir + np.ma.filled((sim_cwu_green_sc3 / sim_yields_sc3 * ir_syst_data_sur * 10), 0) + \
                                                np.ma.filled((sim_cwu_green_sc4 / sim_yields_sc4 * ir_syst_data_spr * 10), 0) + \
                                                np.ma.filled((sim_cwu_green_sc5 / sim_yields_sc5 * ir_syst_data_drp * 10), 0)  # Green WF
                                
                    wf_blue_ir = wf_blue_ir   + np.ma.filled((sim_cwu_blue_ir_sc3 / sim_yields_sc3 * ir_syst_data_sur * 10), 0) + \
                                                np.ma.filled((sim_cwu_blue_ir_sc4 / sim_yields_sc4 * ir_syst_data_spr * 10), 0) + \
                                                np.ma.filled((sim_cwu_blue_ir_sc5 / sim_yields_sc5 * ir_syst_data_drp * 10), 0) # Blue from irrigation
            else: # only one strategy is considered
                if 3 in self.conf.scenarios:
                    print('Scenario 3...')
                    sim_yields_sc3or6, sim_cwu_green_sc3or6, sim_cwu_blue_sc3or6, _ = ac.GetOutputsOfScenario(self.conf, 3) # Irrigated furrow
                if 6 in self.conf.scenarios:
                    print('Scenario 6...')
                    sim_yields_sc3or6, sim_cwu_green_sc3or6, sim_cwu_blue_sc3or6, sim_cwu_blue_cr_sc6 = ac.GetOutputsOfScenario(self.conf, 6) # Irrigated flooding with gw
                    
                sim_yields_sc3or6 = np.ma.masked_where(sim_yields_sc3or6 == 0, sim_yields_sc3or6)
                if len(sim_yields_sc3or6[sim_yields_sc3or6.mask==False])>0:
                    _limit = np.percentile(sim_yields_sc3or6[sim_yields_sc3or6.mask==False],1) # Filter out extremely low yield data points
                    sim_yields_sc3or6[sim_yields_sc3or6 < _limit] = _limit
                
                wf_green_ir = wf_green_ir + np.ma.filled((sim_cwu_green_sc3or6 / sim_yields_sc3or6 * 10), 0) # Green WF
                wf_blue_ir = wf_blue_ir + np.ma.filled((sim_cwu_blue_sc3or6 / sim_yields_sc3or6 * 10), 0) # Blue from irrigation
                
                #Save ir CR raster for Rice 
                if 6 in self.conf.scenarios:   
                    wf_blue_ir_cr = np.ma.filled((sim_cwu_blue_cr_sc6 / sim_yields_sc3or6 * 10), 0) # Blue from cr
                    wf_blue_ir_cr = np.ma.masked_where(wf_blue_ir_cr==0, wf_blue_ir_cr)
                    _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_cr_ir_global_annual_{self.text_years}.nc'])
                    val_name = f'wf_blue_cr_ir-{self.conf.crop_name_short}'
                    ac.CreateHistRaster5(wf_blue_ir_cr, _raster_path, val_name, 'm3 t-1 y-1', raster_title, self.first_year, self.max_gs)
        
            wf_green_ir = np.ma.masked_where(wf_green_ir==0, wf_green_ir); wf_blue_ir = np.ma.masked_where(wf_blue_ir==0, wf_blue_ir)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_ir_global_annual_{self.text_years}.nc'])
            val_name = f'wf_green_ir-{self.conf.crop_name_short}'
            ac.CreateHistRaster5(wf_green_ir, _raster_path, val_name, 'm3 t-1 y-1', raster_title, self.first_year, self.max_gs)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_ir_global_annual_{self.text_years}.nc'])
            val_name = f'wf_blue_ir-{self.conf.crop_name_short}'
            ac.CreateHistRaster5(wf_blue_ir, _raster_path, val_name, 'm3 t-1 y-1', raster_title, self.first_year, self.max_gs)
     
    def GetGlobalWFStats(self, save_report=True, save_rasters=True):
        "Get main global statistics"
        
        print('Calculating main global statistics')
        out_report = ac.GetPath(["outputs", "main",'global_wf_stats', f'{self.conf.crop_fao}_WF_GlobalStats_{self.text_years}.txt'])
        if os.path.exists(out_report): os.remove(out_report)
        ac.printLog(out_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        
        # 1. Get ACEA data
        print('1. Getting input data...')
        # Get scaling factors for yields
        scaling_factors = ac.GetScalingFactors(self.conf.crop_fao, self.conf.crop_name_short, self.text_years)
        
        # Get harvest areas
        harv_area_rf = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, False),0); harv_area_ir = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, True),0)
        harv_area = harv_area_ir + harv_area_rf; harv_area = np.ma.masked_where(harv_area == 0, harv_area)
        harv_area_rf_frac = harv_area_rf/harv_area; harv_area_ir_frac = harv_area_ir/harv_area
        
        # Get production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_rf = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Rainfed production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_ir = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Irrigated production
        prod_sum = prod_rf + prod_ir;  prod_sum = np.ma.masked_where(prod_sum == 0, prod_sum) # Total production
        prod_rf_frac = prod_rf/prod_sum; prod_ir_frac = prod_ir/prod_sum
        
        #  Get simulation outputs
        sim_yields_sc1, sim_cwu_green_sc1, _, _ = ac.GetOutputsOfScenario(self.conf, 1, False) # Rainfed no GW, True = scaled
        sim_yields_sc2, sim_cwu_green_sc2, _, sim_cwu_blue_cr_sc2 = ac.GetOutputsOfScenario(self.conf, 2, False) # Rainfed with GW
        
        if sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
            ir_syst_data_sur, ir_syst_data_spr, ir_syst_data_drp = fc.GetIrrSystems5() # Get irrigation fractions
            sim_yields_sc3, sim_cwu_green_sc3, sim_cwu_blue_ir_sc3, _ = ac.GetOutputsOfScenario(self.conf, 3, False)  # Irrigated furrow
            sim_yields_sc4, sim_cwu_green_sc4, sim_cwu_blue_ir_sc4, _ = ac.GetOutputsOfScenario(self.conf, 4, False) # Irrigated sprinkler
            
            if 5 in self.conf.scenarios:
                sim_yields_sc5, sim_cwu_green_sc5, sim_cwu_blue_ir_sc5, _ = ac.GetOutputsOfScenario(self.conf, 5, False) # Irrigated drip
                sim_yields_sc5 = sim_yields_sc5 * ir_syst_data_drp; sim_cwu_green_sc5 = sim_cwu_green_sc5 * ir_syst_data_drp; sim_cwu_blue_ir_sc5 = sim_cwu_blue_ir_sc5 * ir_syst_data_drp
            else:
                ir_syst_data_sur[(ir_syst_data_sur==0) & (ir_syst_data_spr==0)] = 1 # Assume surface for only drip cells
                sys_sum = ir_syst_data_sur + ir_syst_data_spr; ir_syst_data_sur = ir_syst_data_sur/sys_sum; ir_syst_data_spr = ir_syst_data_spr/sys_sum
                
            sim_yields_sc3 = sim_yields_sc3 * ir_syst_data_sur; sim_cwu_green_sc3 = sim_cwu_green_sc3 * ir_syst_data_sur; sim_cwu_blue_ir_sc3 = sim_cwu_blue_ir_sc3 * ir_syst_data_sur
            sim_yields_sc4 = sim_yields_sc4 * ir_syst_data_spr; sim_cwu_green_sc4 = sim_cwu_green_sc4 * ir_syst_data_spr; sim_cwu_blue_ir_sc4 = sim_cwu_blue_ir_sc4 * ir_syst_data_spr
        else:
            if 3 in self.conf.scenarios:
                sim_yields_sc3, sim_cwu_green_sc3, sim_cwu_blue_ir_sc3, _ = ac.GetOutputsOfScenario(self.conf, 3, False)  # Irrigated furrow
            if 6 in self.conf.scenarios:
                sim_yields_sc6, sim_cwu_green_sc6, sim_cwu_blue_ir_sc6, sim_cwu_blue_cr_sc6 = ac.GetOutputsOfScenario(self.conf, 6, False)  # Irrigated furrow
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_pirnreq_global_annual_{self.text_years}.nc'])
        ir_req = np.ma.filled(ac.ReadNC(_raster_path, f'pirnreq-{self.conf.crop_name_short}'),0)
        ir_req = np.ma.masked_where(ir_req < 20, ir_req)
        
        # Get WFs
        wf_green_rf, wf_blue_rf, wf_green_ir, wf_blue_ir_ir, wf_blue_ir_cr = self.GetCropWF()
        if 6 in self.conf.scenarios: wf_blue_ir = wf_blue_ir_ir + wf_blue_ir_cr
        else: wf_blue_ir = wf_blue_ir_ir
        
        # 2. Calculate statistics
        print('2. Calculating statistics...')
        
        #%% Rainfed vs Irrigated
        # Get total unscaled Y
        y_rf = np.ma.filled(sim_yields_sc1, 0) + np.ma.filled(sim_yields_sc2, 0)
        if not sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
            if 3 in self.conf.scenarios: y_ir = np.ma.filled(sim_yields_sc3, 0)
            if 6 in self.conf.scenarios: y_ir = np.ma.filled(sim_yields_sc6, 0)
        else:
            y_ir = np.ma.filled(sim_yields_sc3, 0) + np.ma.filled(sim_yields_sc4, 0)
            if 5 in self.conf.scenarios: y_ir = y_ir + np.ma.filled(sim_yields_sc5, 0)
        
        # Get total CWU
        cwu_rf = np.ma.filled(sim_cwu_green_sc1, 0) + np.ma.filled(sim_cwu_green_sc2, 0) + np.ma.filled(sim_cwu_blue_cr_sc2, 0)
        if not sum(map(lambda x : x>2, self.conf.scenarios)) > 1:
            if 3 in self.conf.scenarios: cwu_ir = np.ma.filled(sim_cwu_green_sc3, 0) + np.ma.filled(sim_cwu_blue_ir_sc3, 0)
            if 6 in self.conf.scenarios:
                cwu_ir = np.ma.filled(sim_cwu_green_sc6, 0) + np.ma.filled(sim_cwu_blue_ir_sc6, 0) + np.ma.filled(sim_cwu_blue_cr_sc6, 0)
                cwu_ir_cr = np.ma.filled(sim_cwu_blue_cr_sc6, 0)
        else: 
            cwu_ir = np.ma.filled(sim_cwu_green_sc3, 0) + np.ma.filled(sim_cwu_blue_ir_sc3, 0) + \
                    np.ma.filled(sim_cwu_green_sc4, 0) + np.ma.filled(sim_cwu_blue_ir_sc4, 0)
            if 5 in self.conf.scenarios: cwu_ir = cwu_ir + np.ma.filled(sim_cwu_green_sc5, 0) + np.ma.filled(sim_cwu_blue_ir_sc5, 0)
        
        y_rf_scaled = y_rf * scaling_factors; y_ir_scaled = y_ir * scaling_factors
        wf_rf = wf_green_rf + wf_blue_rf;  wf_ir = wf_green_ir + wf_blue_ir # Get total rainfed WF
        
        if save_report:
            # Rainfed report
            ac.printLog(out_report, 'Year\tGreen\tBlue\tUnit WF rf[m3 t-1]\tYield[t ha-1]\tScaling factor\tCWU[mm]\tProduction[t]\tHarvested area[ha]')
            for y in range(self.max_gs):
                ac.printLog(out_report, '{}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.3f}\t{:0.2f}\t{:0.2f}\t{:0.2f}'.format(self.first_year+y, \
                        np.ma.average(wf_green_rf[y,:,:], weights = prod_rf[y,:,:]), \
                        np.ma.average(wf_blue_rf[y,:,:], weights = prod_rf[y,:,:]), \
                        np.ma.average(wf_rf[y,:,:], weights = prod_rf[y,:,:]), \
                        np.ma.average(y_rf_scaled[y,:,:], weights = harv_area_rf[y,:,:]), \
                        np.ma.average(scaling_factors[y,:,:], weights = harv_area_rf[y,:,:]), \
                        np.ma.average(cwu_rf[y,:,:], weights = harv_area_rf[y,:,:]), \
                        np.ma.sum(prod_rf[y,:,:]), np.ma.sum(harv_area_rf[y,:,:])
                        ))
            ac.printLog(out_report, '================================')
        
            # Irrigated report
            if 6 in self.conf.scenarios: # Report also blue CR CWU for this scenario
                ac.printLog(out_report, 'Year\tGreen\tBlue\tUnit WF ir[m3 t-1]\tYield[t ha-1]\tScaling factor\tCWU[mm]\tCWU_cr[mm]\tIrrigation[mm]\tProduction[t]\tHarvested area[ha]')
                for y in range(self.max_gs):
                    ac.printLog(out_report, '{}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.3f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}'.format(self.first_year+y, \
                            np.ma.average(wf_green_ir[y,:,:], weights = prod_ir[y,:,:]), \
                            np.ma.average(wf_blue_ir[y,:,:], weights = prod_ir[y,:,:]), \
                            np.ma.average(wf_ir[y,:,:], weights = prod_ir[y,:,:]), \
                            np.ma.average(y_ir_scaled[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.average(scaling_factors[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.average(cwu_ir[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.average(cwu_ir_cr[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.average(ir_req[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.sum(prod_ir[y,:,:]), np.ma.sum(harv_area_ir[y,:,:])
                            ))
            else: 
                ac.printLog(out_report, 'Year\tGreen\tBlue\tUnit WF ir[m3 t-1]\tYield[t ha-1]\tScaling factor\tCWU[mm]\tIrrigation[mm]\tProduction[t]\tHarvested area[ha]')
                for y in range(self.max_gs):
                    ac.printLog(out_report, '{}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.3f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}'.format(self.first_year+y, \
                            np.ma.average(wf_green_ir[y,:,:], weights = prod_ir[y,:,:]), \
                            np.ma.average(wf_blue_ir[y,:,:], weights = prod_ir[y,:,:]), \
                            np.ma.average(wf_ir[y,:,:], weights = prod_ir[y,:,:]), \
                            np.ma.average(y_ir_scaled[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.average(scaling_factors[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.average(cwu_ir[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.average(ir_req[y,:,:], weights = harv_area_ir[y,:,:]), \
                            np.ma.sum(prod_ir[y,:,:]), np.ma.sum(harv_area_ir[y,:,:])
                            ))
            ac.printLog(out_report, '================================')
        
        #%% Green vs Blue
        # Get total green WF
        wf_green_rf = wf_green_rf * prod_rf_frac; wf_green_ir = wf_green_ir * prod_ir_frac; wf_blue_ir_ir = wf_blue_ir_ir * prod_ir_frac
        wf_blue_cr = wf_blue_rf * prod_rf_frac
        if 6 in self.conf.scenarios: wf_blue_cr = wf_blue_cr + wf_blue_ir_cr * prod_ir_frac
        
        wf_green = wf_green_rf + wf_green_ir # Add up rainfed and irrigated
        wf_blue = wf_blue_cr + wf_blue_ir_ir # Add up rainfed and irrigated
        wf_tot = wf_green + wf_blue
        
        # WF report
        if save_report: 
            ac.printLog(out_report, 'Year\tGreen\tBlue cr\tBlue ir\tUnit WF[m3 t-1]\tYield[t ha-1]\tScaling factor\tProduction[t]\tHarvested area[ha]')
            for y in range(self.max_gs):
                ac.printLog(out_report, '{}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.2f}\t{:0.3f}\t{:0.2f}\t{:0.2f}'.format(self.first_year+y, \
                            np.ma.average(wf_green[y,:,:], weights = prod_sum[y,:,:]), \
                            np.ma.average(wf_blue_cr[y,:,:], weights = prod_sum[y,:,:]), \
                            np.ma.average(wf_blue_ir_ir[y,:,:], weights = prod_sum[y,:,:]), \
                            np.ma.average(wf_tot[y,:,:], weights = prod_sum[y,:,:]), \
                            np.ma.sum(prod_sum[y,:,:])/np.ma.sum(harv_area[y,:,:]), \
                            np.ma.average(scaling_factors[y,:,:], weights = harv_area[y,:,:]), \
                            np.ma.sum(prod_sum[y,:,:]), np.ma.sum(harv_area[y,:,:])
                        ))
            ac.printLog(out_report, '================================')
        
        # 3. Save rasters
        print('3. Saving rasters...')
        if save_rasters:
            # Yield unscaled
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_yield_unscaled_rf_global_annual_{self.text_years}.nc'])
            val_name = f'yield_unscaled_rf-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: fresh crop yield unscaled'
            ac.CreateHistRaster5(y_rf, _raster_path, val_name, 't ha-1 y-1', title, self.first_year, self.max_gs)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_yield_unscaled_ir_global_annual_{self.text_years}.nc'])
            val_name = f'yield_unscaled_ir-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: fresh crop yield unscaled'
            ac.CreateHistRaster5(y_ir, _raster_path, val_name, 't ha-1 y-1', title, self.first_year, self.max_gs)
            
            y_final = y_rf * harv_area_rf_frac + y_ir * harv_area_ir_frac
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_yield_unscaled_global_annual_{self.text_years}.nc'])
            val_name = f'yield_unscaled-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: fresh crop yield unscaled'
            ac.CreateHistRaster5(y_final, _raster_path, val_name, 't ha-1 y-1', title, self.first_year, self.max_gs)
            
            # Yield scaled
            # _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_yield_scaled_rf_global_annual_{self.text_years}.nc'])
            # val_name = f'yield_scaled_rf-{self.conf.crop_name_short}'
            # title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: fresh crop yield scaled'
            # ac.CreateHistRaster5(y_rf_scaled, _raster_path, val_name, 't ha-1 y-1', title, self.first_year, self.max_gs)
            
            # _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_yield_scaled_ir_global_annual_{self.text_years}.nc'])
            # val_name = f'yield_scaled_ir-{self.conf.crop_name_short}'
            # title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: fresh crop yield scaled'
            # ac.CreateHistRaster5(y_ir_scaled, _raster_path, val_name, 't ha-1 y-1', title, self.first_year, self.max_gs)
            
            # y_final = y_rf_scaled * harv_area_rf_frac + y_ir_scaled * harv_area_ir_frac
            # _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_yield_scaled_global_annual_{self.text_years}.nc'])
            # val_name = f'yield_scaled-{self.conf.crop_name_short}'
            # title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: fresh crop yield scaled'
            # ac.CreateHistRaster5(y_final, _raster_path, val_name, 't ha-1 y-1', title, self.first_year, self.max_gs)
            
            # WFtot
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_tot_global_annual_{self.text_years}.nc'])
            val_name = f'wf_total-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit water footprint'
            ac.CreateHistRaster5(np.ma.array(wf_tot, mask=(wf_tot==0)), _raster_path, val_name, 'm3 t-1 y-1', title, self.first_year, self.max_gs)
            
            # wf_blue_cr
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_tot_blue_cr_global_annual_{self.text_years}.nc'])
            val_name = f'wf_tot_blue_cr-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit water footprint'
            ac.CreateHistRaster5(np.ma.array(wf_blue_cr, mask=(wf_tot==0)), _raster_path, val_name, 'm3 t-1 y-1', title, self.first_year, self.max_gs)
            
            # wf_blue_ir_ir
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_tot_blue_ir_global_annual_{self.text_years}.nc'])
            val_name = f'wf_tot_blue_ir-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit water footprint'
            ac.CreateHistRaster5(np.ma.array(wf_blue_ir_ir, mask=(wf_tot==0)), _raster_path, val_name, 'm3 t-1 y-1', title, self.first_year, self.max_gs)
            
            # wf_green
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_tot_green_global_annual_{self.text_years}.nc'])
            val_name = f'wf_tot_green-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit water footprint'
            ac.CreateHistRaster5(np.ma.array(wf_green, mask=(wf_tot==0)), _raster_path, val_name, 'm3 t-1 y-1', title, self.first_year, self.max_gs)
    
            
    def AnalyseAllCountriesWF_annual(self):
        "Report the average annual WFs of each country"

        print('Report the average annual WFs of each country')
        print('1. Getting data...')
        # Get production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_rf = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Rainfed production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_ir = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Irrigated production
        prod_sum = prod_rf + prod_ir;  prod_sum = np.ma.masked_where(prod_sum == 0, prod_sum) # Total production
        prod_rf_frac = prod_rf/prod_sum; prod_ir_frac = prod_ir/prod_sum
        faostat_prod = pd.read_csv(ac.GetPath(['data','footprints','faostat','production_1990_2019.csv']), usecols = ['Area Code','Crop Code','Year','Production(t)'])
        faostat_prod = faostat_prod.loc[faostat_prod['Crop Code']==self.conf.crop_fao]
        faostat_ha = pd.read_csv(ac.GetPath(['data','footprints', 'faostat', 'harvested_area_1990_2019.csv']), usecols = ['Area Code','Crop Code','Year','Harvested area(ha)'])
        faostat_ha = faostat_ha.loc[faostat_ha['Crop Code']==self.conf.crop_fao]
        
        # Get harvest areas
        harv_area_rf = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, False),0); harv_area_ir = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, True),0)
        harv_area = harv_area_ir + harv_area_rf; harv_area = np.ma.masked_where(harv_area == 0, harv_area)
        # harv_area_rf_frac = harv_area_rf/harv_area; harv_area_ir_frac = harv_area_ir/harv_area
        harv_area_ir = np.ma.masked_where(harv_area_ir == 0, harv_area_ir)
        
        # Get WFs
        wf_green_rf, wf_blue_rf, wf_green_ir, wf_blue_ir_ir, wf_blue_ir_cr = self.GetCropWF()
        # if 6 in self.conf.scenarios: wf_blue_ir = wf_blue_ir_ir + wf_blue_ir_cr
        # else: wf_blue_ir = wf_blue_ir_ir
        
        print('2. Calculating statistics...')
        wf_green_rf = wf_green_rf * prod_rf_frac; wf_green_ir = wf_green_ir * prod_ir_frac; wf_blue_ir_ir = wf_blue_ir_ir * prod_ir_frac
        wf_blue_cr = wf_blue_rf * prod_rf_frac
        if 6 in self.conf.scenarios: wf_blue_cr = wf_blue_cr + wf_blue_ir_cr * prod_ir_frac
        
        # Get total WF
        wf_green = wf_green_rf + wf_green_ir # Add up rainfed and irrigated
        wf_blue = wf_blue_cr + wf_blue_ir_ir # Add up rainfed and irrigated
        wf_tot = wf_green + wf_blue
        
        print('3. Writing report...')
        
        cols = ['Name', 'FAO_code', 'Year', 'WFg[m3 t-1]', 'WFbcr[m3 t-1]', 'WFbi[m3 t-1]', 'WFt[m3 t-1]', 'Production[t]', 'Harvested area[ha]', 'Irrigated area[fraction]', 'Simulated yield[t ha-1]', 'FAOSTAT yield[t ha-1]']
        final_data = pd.DataFrame(columns=cols)
        # Iterate countries
        country_data = fc.GetListOfHistCountries(); c_count = 0
        for c_code in country_data['Country Code']:
            c_count += 1
            if c_count % 10 == 0: print(f'Analysing countries {int(c_count/len(country_data)*100)}%')
            
            # Get country mask
            gridcells = fc.GetHistCountryCells5(c_code).values
            c_mask = np.ones(np.shape(harv_area)[1:])*np.ma.masked
            for r in gridcells: c_mask[r[1],r[2]] = 1
            if (c_mask.mask==True).all(): continue
         
            # Cut off outspace
            xlim = [min(gridcells[:,1]), max(gridcells[:,1])+1]; ylim = [min(gridcells[:,2]), max(gridcells[:,2])+1]
            c_mask = c_mask[xlim[0]:xlim[1],ylim[0]:ylim[1]]
            prod_sum_masked = prod_sum[:,xlim[0]:xlim[1],ylim[0]:ylim[1]]*c_mask
            if (prod_sum_masked.mask==True).all(): continue
        
            c_name = country_data.loc[country_data['Country Code'] == c_code]['Country'].values[0].replace("\n", "")
            wf_g = np.ma.average(wf_green[:,xlim[0]:xlim[1],ylim[0]:ylim[1]]*c_mask, axis = (1,2), weights = prod_sum_masked)
            wf_bcr = np.ma.average(wf_blue_cr[:,xlim[0]:xlim[1],ylim[0]:ylim[1]]*c_mask, axis = (1,2), weights = prod_sum_masked)
            wf_bi = np.ma.average(wf_blue_ir_ir[:,xlim[0]:xlim[1],ylim[0]:ylim[1]]*c_mask, axis = (1,2), weights = prod_sum_masked)
            wf_t = np.ma.average(wf_tot[:,xlim[0]:xlim[1],ylim[0]:ylim[1]]*c_mask, axis = (1,2), weights = prod_sum_masked)
            prod = np.ma.sum(prod_sum_masked, axis = (1,2)); harv = np.ma.sum(harv_area[:,xlim[0]:xlim[1],ylim[0]:ylim[1]]*c_mask, axis = (1,2))
            harv_ir = np.ma.sum(harv_area_ir[:,xlim[0]:xlim[1],ylim[0]:ylim[1]]*c_mask, axis = (1,2)); harv_ir_frac = harv_ir/harv;
            
            c_end_year = country_data.loc[country_data['Country Code'] == c_code]['End Year'].values[0]
            c_start_year = country_data.loc[country_data['Country Code'] == c_code]['Start Year'].values[0]
            if np.isnan(c_end_year): c_end_year = 2100
            if np.isnan(c_start_year): c_start_year = 1900
            
            for y in range(len(wf_tot)):
                cur_year = self.first_year+y
                if (cur_year <= c_end_year) and (cur_year >= c_start_year) and wf_t[y] > 0:
                    _faostat_prod = faostat_prod[(faostat_prod['Area Code']==c_code) & (faostat_prod['Year']==cur_year)]['Production(t)'].values
                    _faostat_ha = faostat_ha.loc[(faostat_ha['Area Code']==c_code) & (faostat_ha['Year']==cur_year)]['Harvested area(ha)'].values
                    if len(_faostat_prod)==0 or len(_faostat_ha)==0: fao_yield = 0
                    else: fao_yield = (_faostat_prod/_faostat_ha)[0]
                    
                    final_data.loc[len(final_data)] = [c_name, c_code, cur_year, wf_g[y], wf_bcr[y], wf_bi[y], wf_t[y], prod[y], harv[y], harv_ir_frac[y], prod[y]/harv[y], fao_yield]    
            
        # Save the final table
        out_report = ac.GetPath(["outputs", 'main','national_wf_stats_hist_borders', f'{self.conf.crop_fao}_WF_HistNations_annual_{self.text_years}.csv'])
        if os.path.exists(out_report): os.remove(out_report)
        final_data.to_csv(out_report, index=False)
        
    def AnalyseAllCountriesWF_avg(self, avg_period = 3):
        "Report the avg WFs of each country"
        
        print('1. Getting data...')
        # Get production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_rf = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Rainfed production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_ir = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Irrigated production
        prod_sum = prod_rf + prod_ir;  prod_sum = np.ma.masked_where(prod_sum == 0, prod_sum) # Total production
        prod_rf_frac = prod_rf/prod_sum; prod_ir_frac = prod_ir/prod_sum
        
        # Get harvest areas
        harv_area_rf = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, False),0); harv_area_ir = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, True),0)
        harv_area = harv_area_ir + harv_area_rf; harv_area = np.ma.masked_where(harv_area == 0, harv_area)
        # harv_area_rf_frac = harv_area_rf/harv_area; harv_area_ir_frac = harv_area_ir/harv_area
        
        # Get WFs
        wf_green_rf, wf_blue_rf, wf_green_ir, wf_blue_ir_ir, wf_blue_ir_cr = self.GetCropWF()
        # if 6 in self.conf.scenarios: wf_blue_ir = wf_blue_ir_ir + wf_blue_ir_cr
        # else: wf_blue_ir = wf_blue_ir_ir
        
        print('2. Calculating statistics...')
        wf_green_rf = wf_green_rf * prod_rf_frac; wf_green_ir = wf_green_ir * prod_ir_frac; wf_blue_ir_ir = wf_blue_ir_ir * prod_ir_frac
        wf_blue_cr = wf_blue_rf * prod_rf_frac
        if 6 in self.conf.scenarios: wf_blue_cr = wf_blue_cr + wf_blue_ir_cr * prod_ir_frac
        
        # Get total WF
        wf_green = wf_green_rf + wf_green_ir # Add up rainfed and irrigated
        wf_blue = wf_blue_cr + wf_blue_ir_ir # Add up rainfed and irrigated
        wf_tot = wf_green + wf_blue
        
        print('3. Writing report...')
        # Write report
        out_report = ac.GetPath(["outputs", f"{self.conf.project_name}",'reports', f'{self.conf.crop_fao}_WF_NationalReport_{self.text_years}_avg.txt'])
        if os.path.exists(out_report): os.remove(out_report)
        ac.printLog(out_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        ac.printLog(out_report, f'{self.conf.crop_name} Average Water Footprints[m3 t-1]: {self.first_year}-{self.first_year+avg_period-1} vs around 2000 vs {self.clock_end_year-avg_period+1}-{self.clock_end_year}')
        ac.printLog(out_report, 'Country\tFAO_code\tPast unit WF[m3 t-1]\tPast green\tPast blue cr\tPast blue ir\tPast production[t]\tPast harvested area[ha]\
                    \t2000 unit WF[m3 t-1]\t2000 green\t2000 blue cr\t2000 blue ir\t2000 production[t]\t2000 harvested area[ha]\
                    \tRecent unit WF[m3 t-1]\tRecent green\tRecent blue cr\tRecent blue ir\tRecent production[t]\tRecent harvested area[ha]')
        
        # Iterate countries
        country_data = fc.GetListOfCurCountries()
        for c_code in country_data['Country Code']:
            gridcells = fc.GetHistCountryCells5(c_code).values
            
            _wf_green = np.ma.ones((len(gridcells), avg_period)); _wf_green[:,:] = np.ma.masked
            _wf_blue_cr, _wf_blue_ir, _wf_tot, _prod, _harv_area = (_wf_green *1 for _ in range(5))
            wf_green_avg, wf_blue_cr_avg, wf_blue_ir_avg, wf_tot_avg, prod_avg, harv_avg = (np.zeros((3)) for _ in range(6))
            
            # Getting past numbers
            count = 0
            for c5 in gridcells:
                row5 = c5[1]; col5 = c5[2]
                _wf_green[count,:] = wf_green[:avg_period,row5,col5];_wf_blue_ir[count,:] = wf_blue_ir_ir[:avg_period,row5,col5]; _wf_blue_cr[count,:] = wf_blue_cr[:avg_period,row5,col5]
                _wf_tot[count,:] = wf_tot[:avg_period,row5,col5]; _prod[count,:] = prod_sum[:avg_period,row5,col5]; _harv_area[count,:] = harv_area[:avg_period,row5,col5]
                count += 1
            
            if (_wf_green.count() > 0) or (_wf_blue_ir.count() > 0):
                wf_green_avg[0] = np.ma.average(np.ma.average(_wf_green, axis = 0, weights = _prod))
                wf_blue_cr_avg[0] = np.ma.average(np.ma.average(_wf_blue_cr, axis = 0, weights = _prod))
                wf_blue_ir_avg[0] = np.ma.average(np.ma.average(_wf_blue_ir, axis = 0, weights = _prod))
                wf_tot_avg[0] = np.ma.average(np.ma.average(_wf_tot, axis = 0, weights = _prod))
                prod_avg[0] = np.ma.average(np.ma.sum(_prod, axis = 0)); harv_avg[0] = np.ma.average(np.ma.sum(_harv_area, axis = 0))
                
            _wf_green = np.ma.ones((len(gridcells), avg_period)); _wf_green[:,:] = np.ma.masked
            _wf_blue_cr, _wf_blue_ir, _wf_tot, _prod, _harv_area = (_wf_green *1 for _ in range(5))
             
            # Getting 2000 numbers
            count = 0
            for c5 in gridcells:
                row5 = c5[1]; col5 = c5[2]
                _wf_green[count,:] = wf_green[9:12,row5,col5];_wf_blue_ir[count,:] = wf_blue_ir_ir[9:12,row5,col5]; _wf_blue_cr[count,:] = wf_blue_cr[9:12,row5,col5]
                _wf_tot[count,:] = wf_tot[9:12,row5,col5]; _prod[count,:] = prod_sum[9:12,row5,col5]; _harv_area[count,:] = harv_area[9:12,row5,col5]
                count += 1
            
            if (_wf_green.count() > 0) or (_wf_blue_ir.count() > 0):
                wf_green_avg[1] = np.ma.average(np.ma.average(_wf_green, axis = 0, weights = _prod))
                wf_blue_cr_avg[1] = np.ma.average(np.ma.average(_wf_blue_cr, axis = 0, weights = _prod))
                wf_blue_ir_avg[1] = np.ma.average(np.ma.average(_wf_blue_ir, axis = 0, weights = _prod))
                wf_tot_avg[1] = np.ma.average(np.ma.average(_wf_tot, axis = 0, weights = _prod))
                prod_avg[1] = np.ma.average(np.ma.sum(_prod, axis = 0)); harv_avg[1] = np.ma.average(np.ma.sum(_harv_area, axis = 0))
                
            _wf_green = np.ma.ones((len(gridcells), avg_period)); _wf_green[:,:] = np.ma.masked
            _wf_blue_cr, _wf_blue_ir, _wf_tot, _prod, _harv_area = (_wf_green *1 for _ in range(5)) 
            
            # Getting recent numbers
            count = 0
            for c5 in gridcells:
                row5 = c5[1]; col5 = c5[2]
                _wf_green[count,:] = wf_green[-avg_period:,row5,col5];_wf_blue_ir[count,:] = wf_blue_ir_ir[-avg_period:,row5,col5]; _wf_blue_cr[count,:] = wf_blue_cr[-avg_period:,row5,col5]
                _wf_tot[count,:] = wf_tot[-avg_period:,row5,col5]; _prod[count,:] = prod_sum[-avg_period:,row5,col5]; _harv_area[count,:] = harv_area[-avg_period:,row5,col5]
                count += 1
                
            if (_wf_green.count() > 0) or (_wf_blue_ir.count() > 0):
                wf_green_avg[2] = np.ma.average(np.ma.average(_wf_green, axis = 0, weights = _prod))
                wf_blue_cr_avg[2] = np.ma.average(np.ma.average(_wf_blue_cr, axis = 0, weights = _prod))
                wf_blue_ir_avg[2] = np.ma.average(np.ma.average(_wf_blue_ir, axis = 0, weights = _prod))
                wf_tot_avg[2] = np.ma.average(np.ma.average(_wf_tot, axis = 0, weights = _prod))
                prod_avg[2] = np.ma.average(np.ma.sum(_prod, axis = 0)); harv_avg[2] = np.ma.average(np.ma.sum(_harv_area, axis = 0))
            
            if any(wf_green_avg > 0) or any(wf_blue_ir_avg > 0):
                ac.printLog(out_report, "{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                         country_data.loc[country_data['Country Code'] == c_code]['Country'].values[0].replace("\n", ""), c_code,
                         wf_tot_avg[0], wf_green_avg[0], wf_blue_cr_avg[0], wf_blue_ir_avg[0], prod_avg[0], harv_avg[0],
                         wf_tot_avg[1], wf_green_avg[1], wf_blue_cr_avg[1], wf_blue_ir_avg[1], prod_avg[1], harv_avg[1],
                         wf_tot_avg[2], wf_green_avg[2], wf_blue_cr_avg[2], wf_blue_ir_avg[2], prod_avg[2], harv_avg[2],
                    ))
    
    def AnalyseRegionalWF(self):
        "Analyse the WFs of sub regions"
        print("Analysing the WFs of regions")
        
        out_subreg_report = ac.GetPath(["outputs", f"{self.conf.project_name}",'reports', f'{self.conf.crop_fao}_WF_SubRegionalStats_{self.text_years}.txt']) # Log the process
        ac.printLog(out_subreg_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        if os.path.exists(out_subreg_report): os.remove(out_subreg_report)
        
        out_reg_report = ac.GetPath(["outputs", f"{self.conf.project_name}",'reports', f'{self.conf.crop_fao}_WF_RegionalStats_{self.text_years}.txt'])
        ac.printLog(out_reg_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        if os.path.exists(out_reg_report): os.remove(out_reg_report)
        
        out_reg_hist_report= ac.GetPath(["outputs", f"{self.conf.project_name}",'reports', f'{self.conf.crop_fao}_WF_RegionalStats_Hist_{self.text_years}.txt'])
        ac.printLog(out_reg_hist_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        if os.path.exists(out_reg_hist_report): os.remove(out_reg_hist_report)
        
        regions = ac.GetRegionList()
        
        # Get WFs
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_rf_global_annual_{self.text_years}.nc'])
        wf_green_rf = np.ma.filled(ac.ReadNC(_raster_path, f'wf_green_rf-{self.conf.crop_name_short}'), 0) # Rainfed green WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_rf_global_annual_{self.text_years}.nc']) 
        wf_blue_rf = np.ma.filled(ac.ReadNC(_raster_path, f'wf_blue_rf-{self.conf.crop_name_short}'), 0) # Rainfed blue capillary rise WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_ir_global_annual_{self.text_years}.nc'])
        wf_green_ir = np.ma.filled(ac.ReadNC(_raster_path, f'wf_green_ir-{self.conf.crop_name_short}'), 0)# Irrigated green WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_ir_global_annual_{self.text_years}.nc'])
        wf_blue_ir = np.ma.filled(ac.ReadNC(_raster_path, f'wf_blue_ir-{self.conf.crop_name_short}'), 0) # Irrigated blue irrigated WF
        
        # Get production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_rf = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Rainfed production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_ir = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Irrigated production
        prod_sum = prod_rf + prod_ir;  prod_sum = np.ma.masked_where(prod_sum == 0, prod_sum) # Total production
        prod_rf_frac = prod_rf/prod_sum; prod_ir_frac = prod_ir/prod_sum
        
        
        wf_green_rf = wf_green_rf * prod_rf_frac; wf_blue_rf = wf_blue_rf * prod_rf_frac
        wf_green_ir = wf_green_ir * prod_ir_frac; wf_blue_ir = wf_blue_ir * prod_ir_frac
        
        wf_green = wf_green_rf + wf_green_ir # Add up rainfed and irr
        wf_blue = wf_blue_rf + wf_blue_ir # Add up rainfed and irr
        wf_tot = wf_green + wf_blue
        
        # Get harvest areas
        harv_area_rf = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, False),0); harv_area_ir = np.ma.filled(fc.GetHistHarvestedAreas5(self.conf.crop_fao, True),0)
        harv_area = harv_area_ir + harv_area_rf; harv_area = np.ma.masked_where(harv_area == 0, harv_area)
        harv_area_rf_frac = harv_area_rf/harv_area; harv_area_ir_frac = harv_area_ir/harv_area

        #  Get simulation outputs
        sim_yields_sc1, sim_cwu_green_sc1, _, _ = ac.GetOutputsOfScenario(self.conf, 1, False) # Rainfed no GW, True = scaled
        sim_yields_sc2, sim_cwu_green_sc2, _, sim_cwu_blue_cr_sc2 = ac.GetOutputsOfScenario(self.conf, 2, False) # Rainfed with GW
        
        sim_yields_sc3, sim_cwu_green_sc3, sim_cwu_blue_ir_sc3, _ = ac.GetOutputsOfScenario(self.conf, 3, False) # Irrigated surface
        sim_yields_sc4, sim_cwu_green_sc4, sim_cwu_blue_ir_sc4, _ = ac.GetOutputsOfScenario(self.conf, 4, False) # Irrigated sprinkler
        ir_syst_data_sur, ir_syst_data_spr,_ = fc.GetIrrSystems5() # Get irrigation fractions
        sim_yields_sc3 = sim_yields_sc3 * ir_syst_data_sur; sim_cwu_green_sc3 = sim_cwu_green_sc3 * ir_syst_data_sur; sim_cwu_blue_ir_sc3 = sim_cwu_blue_ir_sc3 * ir_syst_data_sur;
        sim_yields_sc4 = sim_yields_sc4 * ir_syst_data_spr; sim_cwu_green_sc4 = sim_cwu_green_sc4 * ir_syst_data_spr; sim_cwu_blue_ir_sc4 = sim_cwu_blue_ir_sc4 * ir_syst_data_spr;
        
        # Get scaling factors for yields
        _raster_path = ac.GetPath([self.data_folder,'scaling_factors', f'scaling_factors_5arc_{self.conf.crop_fao}_moving_avg_global_annual_{self.text_years}.nc'])
        scaling_factors = ac.ReadNC(_raster_path, f'scaling_factor-{self.conf.crop_name_short}')
        
        y_rf = np.ma.filled(sim_yields_sc1, 0) + np.ma.filled(sim_yields_sc2, 0)
        y_ir = np.ma.filled(sim_yields_sc3, 0) + np.ma.filled(sim_yields_sc4, 0)
        y_tot_unscaled = (y_rf * harv_area_rf_frac + y_ir * harv_area_ir_frac)
        y_tot = y_tot_unscaled * scaling_factors
        
        cwu_rf = np.ma.filled(sim_cwu_green_sc1, 0) + np.ma.filled(sim_cwu_green_sc2, 0) + np.ma.filled(sim_cwu_blue_cr_sc2, 0)
        cwu_ir = np.ma.filled(sim_cwu_green_sc3, 0) + np.ma.filled(sim_cwu_blue_ir_sc3, 0) + np.ma.filled(sim_cwu_green_sc4, 0) + np.ma.filled(sim_cwu_blue_ir_sc4, 0)
        cwu_tot =  cwu_rf * harv_area_rf_frac + cwu_ir * harv_area_ir_frac
        
        # Coeffient of variation
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_tot_trended_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_tot_trended = ac.ReadNC(_raster_path, f'wf_total-{self.conf.crop_name_short}')
        else: wf_tot_trended = self.TrendRasterData(wf_tot, 1, _raster_path)
        wf_tot_detrended = wf_tot[0,:,:] + (wf_tot - wf_tot_trended); wf_tot_detrended[wf_tot_detrended<0] = 0
        std_wf_tot = np.ma.std(wf_tot_detrended, axis=0); mean_wf_tot = np.ma.mean(wf_tot_detrended, axis = 0); cv_wf_tot = std_wf_tot/mean_wf_tot*100
        prod_avg = np.ma.average(prod_sum[-5:,:,:], axis=0)
                                  
        # Get region cells
        ac.printLog(out_subreg_report, 'Sub-region\tHarv. area[ha]\tProduction[t]\tIrrigated[%]\tCrop yield[t ha-1]\tCrop yield unscaled[t ha-1]\tCWU[mm]\tWFg\tWFbc\tWFbi\tWFtot_per1[m3 t-1]\tDetrended CV of WFt[%]')
        ac.printLog(out_reg_report, 'Region\tHarv. area[ha]\tProduction[t]\tIrrigated[%]\tCrop yield[t ha-1]\tCrop yield unscaled[t ha-1]\tCWU[mm]\tWFg\tWFbc\tWFbi\tWFtot_per1[m3 t-1]\tDetrended CV of WFt[%]')
        ac.printLog(out_reg_hist_report, 'Region\tHarv area rf\tHarv area ir[ha]\tProduction[t]\tWFtot[m3 t-1]\tCWU[mm]')
        
        for region in regions:
            sub_regions = ac.GetSubRegionListByID(int(region[0]))
            prod1, wf_tot1, wf_green3, yield3, yield_unsc3, cwu3, wf_blue_rf3, wf_blue_ir3, prod3, prod_ir3, harv3 = \
                (np.zeros((5, len(sub_regions))) for _ in range(11))
                
            hist_area_rf, hist_area_ir, hist_wf, hist_prod, hist_cwu = (np.zeros((31, len(sub_regions))) for _ in range(5))
            cv_tot3 = np.zeros(len(sub_regions))
            subr = 0
            for sub_region in sub_regions:
                _, region_cells = ac.GetSubRegionCells5(sub_region[0], sub_region[1])
                cv_tot3[subr] = np.ma.average(cv_wf_tot[region_cells==1], weights = prod_avg[region_cells==1])
                
                idx = 0
                for y in range(0,5):# 1986-1990
                    wf_tot_temp = wf_tot[y,:,:]; prod_temp = prod_sum[y,:,:]
                    wf_tot1[idx, subr] = np.ma.average(wf_tot_temp[region_cells==1], weights = prod_temp[region_cells==1])
                    prod1[idx, subr] = np.ma.sum(prod_temp[region_cells==1])
                    idx += 1
                
                idx = 0
                for y in range(26,31):# 2012-2016
                    prod_temp = prod_sum[y,:,:]; prod_ir_temp = prod_ir[y,:,:]
                    yield_temp = y_tot[y,:,:]; yield_unsc_temp = y_tot_unscaled[y,:,:]; cwu_temp = cwu_tot[y,:,:]; harv_temp = harv_area[y,:,:];
                    wf_green_temp = wf_green[y,:,:]; wf_blue_rf_temp = wf_blue_rf[y,:,:]; wf_blue_ir_temp = wf_blue_ir[y,:,:]
                    
                    prod3[idx, subr] = np.ma.sum(prod_temp[region_cells==1])
                    harv3[idx, subr] = np.ma.sum(harv_temp[region_cells==1])
                    prod_ir3[idx, subr] = np.ma.sum(prod_ir_temp[region_cells==1])/prod3[idx, subr]*100
                    yield3[idx, subr] = np.ma.average(yield_temp[region_cells==1], weights = harv_temp[region_cells==1])
                    yield_unsc3[idx, subr] = np.ma.average(yield_unsc_temp[region_cells==1], weights = harv_temp[region_cells==1])
                    cwu3[idx, subr] = np.ma.average(cwu_temp[region_cells==1], weights = harv_temp[region_cells==1])
                    wf_green3[idx, subr] = np.ma.average(wf_green_temp[region_cells==1], weights = prod_temp[region_cells==1])
                    wf_blue_rf3[idx, subr] = np.ma.average(wf_blue_rf_temp[region_cells==1], weights = prod_temp[region_cells==1])
                    wf_blue_ir3[idx, subr] = np.ma.average(wf_blue_ir_temp[region_cells==1], weights = prod_temp[region_cells==1])
                    idx += 1
                
                sub_region_name = sub_region[3] if isinstance(sub_region[3], str) else sub_region[2]
                ac.printLog(out_subreg_report, '{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(sub_region_name, \
                    np.nanmean(harv3[:, subr]), np.nanmean(prod3[:, subr]), np.nanmean(prod_ir3[:, subr]),\
                    np.nanmean(yield3[:, subr]), np.nanmean(yield_unsc3[:, subr]), np.nanmean(cwu3[:, subr]),\
                    np.nanmean(wf_green3[:, subr]), np.nanmean(wf_blue_rf3[:, subr]), np.nanmean(wf_blue_ir3[:, subr]), np.nanmean(wf_tot1[:, subr]),\
                    cv_tot3[subr]
                ))
                
                for y in range(31):# 1986-2016
                    harv_rf_temp = harv_area_rf[y,:,:]; harv_ir_temp = harv_area_ir[y,:,:]; harv_temp = harv_area[y,:,:]
                    wf_tot_temp = wf_tot[y,:,:]; prod_temp = prod_sum[y,:,:]; cwu_temp = cwu_tot[y,:,:]
                    
                    hist_area_rf[y, subr] = np.ma.sum(harv_rf_temp[region_cells==1])
                    hist_area_ir[y, subr] = np.ma.sum(harv_ir_temp[region_cells==1])
                    hist_wf[y, subr] = np.ma.average(wf_tot_temp[region_cells==1], weights = prod_temp[region_cells==1])
                    hist_prod[y, subr] = np.ma.sum(prod_temp[region_cells==1])
                    hist_cwu[y, subr] = np.ma.average(cwu_temp[region_cells==1], weights = harv_temp[region_cells==1])
                
                subr += 1
            
            prod1 = np.ma.masked_where(np.isnan(prod1), prod1); wf_tot1 = np.ma.masked_where(np.isnan(wf_tot1), wf_tot1)
            prod3 = np.ma.masked_where(np.isnan(prod3), prod3)
            harv3 = np.ma.masked_where(np.isnan(harv3), harv3)
            prod_ir3 = np.ma.masked_where(np.isnan(prod_ir3), prod_ir3)
            yield3 = np.ma.masked_where(np.isnan(yield3), yield3); yield_unsc3 = np.ma.masked_where(np.isnan(yield_unsc3), yield_unsc3)
            cwu3 = np.ma.masked_where(np.isnan(cwu3), cwu3)
            wf_green3 = np.ma.masked_where(np.isnan(wf_green3), wf_green3)
            wf_blue_rf3 = np.ma.masked_where(np.isnan(wf_blue_rf3), wf_blue_rf3)
            wf_blue_ir3 = np.ma.masked_where(np.isnan(wf_blue_ir3), wf_blue_ir3)
            hist_area_rf = np.ma.masked_where(np.isnan(hist_area_rf), hist_area_rf)
            hist_area_ir = np.ma.masked_where(np.isnan(hist_area_ir), hist_area_ir)
            hist_prod = np.ma.masked_where(np.isnan(hist_prod), hist_prod)
            hist_wf = np.ma.masked_where(np.isnan(hist_wf), hist_wf)
            hist_cwu = np.ma.masked_where(np.isnan(hist_cwu), hist_cwu)
            
            ac.printLog(out_subreg_report, '=============')
            ac.printLog(out_reg_report, '{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(region[1], \
                   np.average(np.ma.sum(harv3, axis = 1)), np.average(np.ma.sum(prod3, axis = 1)), \
                    np.nanmean(np.ma.average(prod_ir3, axis = 1, weights = prod3)),\
                    np.nanmean(np.ma.average(yield3, axis = 1, weights = harv3)),\
                    np.nanmean(np.ma.average(yield_unsc3, axis = 1, weights = harv3)),\
                    np.nanmean(np.ma.average(cwu3, axis = 1, weights = harv3)),\
                    np.nanmean(np.ma.average(wf_green3, axis = 1, weights = prod3)), \
                    np.nanmean(np.ma.average(wf_blue_rf3, axis = 1, weights = prod3)), \
                    np.nanmean(np.ma.average(wf_blue_ir3, axis = 1, weights = prod3)), \
                    np.nanmean(np.ma.average(wf_tot1, axis = 1, weights = prod1)), \
                    np.average(cv_tot3, weights = np.ma.mean(prod3, axis = 0))
            ))
            for y in range(31):# 1986-2016
                ac.printLog(out_reg_hist_report, '{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(region[1], \
                        np.ma.sum(hist_area_rf[y,:]), np.ma.sum(hist_area_ir[y,:]), np.ma.sum(hist_prod[y,:]), \
                        np.ma.average(hist_wf[y,:], weights = hist_prod[y,:]), 
                        np.ma.average(hist_cwu[y,:], weights = (hist_area_rf[y,:]+hist_area_ir[y,:])), 
                ))
    
    def GetCropWF(self):
        
        wf_green_rf, wf_blue_rf, wf_green_ir, wf_blue_ir, wf_blue_ir_cr = [0,0,0,0,0]
        
        # Rainfed green WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_rf_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_green_rf = np.ma.filled(ac.ReadNC(_raster_path, f'wf_green_rf-{self.conf.crop_name_short}'), 0)
        
        # Rainfed blue capillary rise WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_rf_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_blue_rf = np.ma.filled(ac.ReadNC(_raster_path, f'wf_blue_rf-{self.conf.crop_name_short}'), 0)
        
        # Irrigated green WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_ir_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_green_ir = np.ma.filled(ac.ReadNC(_raster_path, f'wf_green_ir-{self.conf.crop_name_short}'), 0)
        
        # Irrigated blue irrigated WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_ir_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_blue_ir = np.ma.filled(ac.ReadNC(_raster_path, f'wf_blue_ir-{self.conf.crop_name_short}'), 0)
        
        # Irrigated blue irrigated WF
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_cr_ir_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_blue_ir_cr = np.ma.filled(ac.ReadNC(_raster_path, f'wf_blue_cr_ir-{self.conf.crop_name_short}'), 0)
        
        wf_blue_ir_cr
        
        return wf_green_rf, wf_blue_rf, wf_green_ir, wf_blue_ir, wf_blue_ir_cr
    
    def TrendRasterData(self, data_ori, save=False, fpath=False, limit=True):
        "Detrend historical raster"
        
        cols = data_ori.shape[2]; years = np.arange(self.first_year, self.clock_end_year+1)
        empty_arr = np.zeros((data_ori.shape))
        data_trended = np.ma.masked_where(empty_arr == 0, empty_arr)
            
        for row in range(0, int(cols/2)):
            print('Detrending row %i' % row)
            for col in range(0, cols):
                data_temp = data_ori[:,row,col]
                if all(np.ma.getmask(data_temp)): continue
                if limit:
                    limit_up = np.percentile(data_temp,90); limit_lo = np.percentile(data_temp,10)
                    med_val = np.percentile(data_temp,50)
                    if limit_up > med_val*2: data_temp[data_temp>limit_up] = limit_up
                    if limit_lo < med_val/2: data_temp[data_temp<limit_lo] = limit_lo
                    
                coef, intercept = np.ma.polyfit(years, data_temp, 1)
                data_trended[:,row,col] = years*coef + intercept
        
        data_trended[data_trended<0] = 0
        
        if save:
            val_name = f'wf-{self.conf.crop_name_short}'; title = 'Trended WF total'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}'
            ac.CreateHistRaster5(data_trended, fpath, val_name, 'm3 t-1 y-1', title)
            
        return data_trended

class lf_analyser:
    "To calculate and analyse crop LFs"
    
    # General settings
    # _,_,res5 = ac.GetResolution(1)
    lf_class = [33, 75]
    data_folder = ac.GetPath(['data','footprints'])
    
    def __init__(self, conf):
        self.conf = conf
        self.raster_folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'rasters'])
        self.clock_start_year = int(self.conf.clock_start.split('/')[0]); self.clock_end_year = int(self.conf.clock_end.split('/')[0])
        self.max_gs = (self.clock_end_year - self.clock_start_year)+1 - self.conf.spinup; self.first_year = self.clock_start_year + self.conf.spinup
        self.text_years = f'{self.first_year}_{self.clock_end_year}'
        
        print(f"============Analysis of {self.conf.crop_name} LFs============")
        print("Getting general data...")
        self.GetGeneralData()
        
    def GetGeneralData(self):
        _limit = .1 # Filter out extremely low data points
        
        # Get harvest areas
        harv_area_rf = fc.GetHistHarvestedAreas5(self.conf.crop_fao, False); harv_area_ir = fc.GetHistHarvestedAreas5(self.conf.crop_fao, True)
        harv_area = np.ma.filled(harv_area_ir,0) + np.ma.filled(harv_area_rf,0)
        self.harv_area = np.ma.masked_where(harv_area <= _limit, harv_area)
        self.harv_area_rf = np.ma.masked_where(harv_area_rf <= _limit, harv_area_rf)
        self.harv_area_ir = np.ma.masked_where(harv_area_ir <= _limit, harv_area_ir)
        
        # Get production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_rf = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Rainfed production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        prod_ir = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Irrigated production
        prod = prod_rf + prod_ir
        self.prod = np.ma.masked_where(prod <= _limit, prod)
        self.prod_rf = np.ma.masked_where(prod_rf <= _limit, prod_rf)
        self.prod_ir = np.ma.masked_where(prod_ir <= _limit, prod_ir)
        
    def CalculateLF(self, save_rasters=True):
        "Calculate LFs"
        
        # 1. Get data
        print("1. Getting more data...")
        
        # 2. Calculate LFs
        print('2. Calculating footprints...')
        unitLF_rf = (self.harv_area_rf*10**4)/self.prod_rf
        unitLF_ir = (self.harv_area_ir*10**4)/self.prod_ir
        unitLF = (self.harv_area*10**4)/self.prod
        
        # 3. Save rasters
        if save_rasters:
            print('3. Saving footprints as rasters...')
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_lf_tot_global_annual_{self.text_years}.nc'])
            val_name = f'lf_total-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit land footprint'
            ac.CreateHistRaster5(unitLF, _raster_path, val_name, 'm2 t-1 y-1', title, self.first_year, self.max_gs)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_lf_rf_global_annual_{self.text_years}.nc'])
            val_name = f'lf_rf-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit land footprint'
            ac.CreateHistRaster5(unitLF_rf, _raster_path, val_name, 'm2 t-1 y-1', title, self.first_year, self.max_gs)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_lf_ir_global_annual_{self.text_years}.nc'])
            val_name = f'lf_ir-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: unit land footprint'
            ac.CreateHistRaster5(unitLF_ir, _raster_path, val_name, 'm2 t-1 y-1', title, self.first_year, self.max_gs)
            
    def GetGlobalLFStats(self, save_report=True):
        "Get global statistics"
        
        print("Getting global statistics...")
        # Get data
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_lf_tot_global_annual_{self.text_years}.nc'])
        unitLF = ac.ReadNC(_raster_path, f'lf_total-{self.conf.crop_name_short}')
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_lf_rf_global_annual_{self.text_years}.nc'])
        unitLF_rf = ac.ReadNC(_raster_path, f'lf_rf-{self.conf.crop_name_short}')
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_lf_ir_global_annual_{self.text_years}.nc'])
        unitLF_ir = ac.ReadNC(_raster_path, f'lf_ir-{self.conf.crop_name_short}')
        # scaling_factors = ac.GetScalingFactors(self.conf.crop_fao, self.conf.crop_name_short, self.text_years)
        land_suit_all = ac.ReadNC(ac.GetPath([self.data_folder,'suitability','overall_crop_suit_5arc_1981_2010.nc']))
        
        #%% Rainfed vs Irrigated
        # Rainfed
        # Header of the report
        _text = 'Year\tLF_HS\tLF_MS\tLF_PS\tUnit LF[m2 t-1]\tProduction[t]\tHarvested area[ha]\tHA_HS\tHA_MS\tHA_PS\tAvgS[%]'
        if save_report:
            out_report = ac.GetPath(["outputs", f"{self.conf.project_name}",'reports', f'{self.conf.crop_fao}_LF_GlobalStats_{self.text_years}.txt'])
            if os.path.exists(out_report): os.remove(out_report)
            ac.printLog(out_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            ac.printLog(out_report, _text)
        else: print(_text)
        
        # Body of the report
        for y in range(self.max_gs):
            _harv_area = self.harv_area_rf[y,:,:]; _harv_area_sum = np.ma.sum(_harv_area)
            _unitLF = np.ma.average(unitLF_rf[y,:,:], weights = self.prod_rf[y,:,:]) # Global avg unit LF[m2 t-1]
            # _y_gap = np.ma.average((1-scaling_factors[y,:,:])*100, weights = _harv_area)
            _unitLF_hs = np.ma.average(unitLF_rf[y,:,:][land_suit_all>self.lf_class[1]], weights = self.prod_rf[y,:,:][land_suit_all>self.lf_class[1]])
            _unitLF_ms = np.ma.average(unitLF_rf[y,:,:][(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])], weights = self.prod_rf[y,:,:][(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])]) 
            _unitLF_ps = np.ma.average(unitLF_rf[y,:,:][(land_suit_all<=self.lf_class[0])], weights = self.prod_rf[y,:,:][(land_suit_all<=self.lf_class[0])])
            
            # Suitability
            _harv_area_hs = np.ma.sum(_harv_area[land_suit_all>self.lf_class[1]])/_harv_area_sum*100
            _harv_area_ms = np.ma.sum(_harv_area[(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])])/_harv_area_sum*100
            _harv_area_ps = np.ma.sum(_harv_area[(land_suit_all<=self.lf_class[0])])/_harv_area_sum*100
            _harv_area_wavg = np.ma.average(land_suit_all, weights = _harv_area)
            
            _text = f"{self.first_year+y}\t{_unitLF_hs:.2f}\t{_unitLF_ms:.2f}\t{_unitLF_ps:.2f}\t{_unitLF:.2f}\t{np.ma.sum(self.prod_rf[y,:,:]):.2f}\t{_harv_area_sum:.2f}\t{_harv_area_hs:.2f}\t{_harv_area_ms:.2f}\t{_harv_area_ps:.2f}\t{_harv_area_wavg:.2f}"
            if save_report: ac.printLog(out_report, _text)
            else: print(_text)

        # Irrigated
        # Header of the report
        _text = 'Year\tLF_HS\tLF_MS\tLF_PS\tUnit LF[m2 t-1]\tProduction[t]\tHarvested area[ha]\tHA_HS\tHA_MS\tHA_PS\tAvgS[%]'
        if save_report:
            ac.printLog(out_report, '================================')
            ac.printLog(out_report, _text)
        else: print(_text)
        
        # Body of the report
        for y in range(self.max_gs):
            _harv_area = self.harv_area_ir[y,:,:]; _harv_area_sum = np.ma.sum(_harv_area)
            _unitLF = np.ma.average(unitLF_ir[y,:,:], weights = self.prod_ir[y,:,:]) # Global avg unit LF[m2 t-1]
            # _y_gap = np.ma.average((1-scaling_factors[y,:,:])*100, weights = _harv_area)
            _unitLF_hs = np.ma.average(unitLF_ir[y,:,:][land_suit_all>self.lf_class[1]], weights = self.prod_ir[y,:,:][land_suit_all>self.lf_class[1]])
            _unitLF_ms = np.ma.average(unitLF_ir[y,:,:][(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])], weights = self.prod_ir[y,:,:][(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])]) 
            _unitLF_ps = np.ma.average(unitLF_ir[y,:,:][(land_suit_all<=self.lf_class[0])], weights = self.prod_ir[y,:,:][(land_suit_all<=self.lf_class[0])])
            
            # Suitability
            _harv_area_hs = np.ma.sum(_harv_area[land_suit_all>self.lf_class[1]])/_harv_area_sum*100
            _harv_area_ms = np.ma.sum(_harv_area[(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])])/_harv_area_sum*100
            _harv_area_ps = np.ma.sum(_harv_area[(land_suit_all<=self.lf_class[0])])/_harv_area_sum*100
            _harv_area_wavg = np.ma.average(land_suit_all, weights = _harv_area)
            
            _text = f"{self.first_year+y}\t{_unitLF_hs:.2f}\t{_unitLF_ms:.2f}\t{_unitLF_ps:.2f}\t{_unitLF:.2f}\t{np.ma.sum(self.prod_ir[y,:,:]):.2f}\t{_harv_area_sum:.2f}\t{_harv_area_hs:.2f}\t{_harv_area_ms:.2f}\t{_harv_area_ps:.2f}\t{_harv_area_wavg:.2f}"
            if save_report: ac.printLog(out_report, _text)
            else: print(_text)
            
        #%% Total
        # Header of the report
        _text = 'Year\tLF_HS\tLF_MS\tLF_PS\tUnit LF[m2 t-1]\tProduction[t]\tHarvested area[ha]\tHA_HS\tHA_MS\tHA_PS\tAvgS[%]'
        if save_report: 
            ac.printLog(out_report, '================================')
            ac.printLog(out_report, _text)
        else: print(_text)
            
        # Body of the report
        for y in range(self.max_gs):
            _harv_area = self.harv_area[y,:,:]; _harv_area_sum = np.ma.sum(_harv_area)
            _unitLF = np.ma.average(unitLF[y,:,:], weights = self.prod[y,:,:]) # Global avg unit LF[m2 t-1]
            # _y_gap = np.ma.average((1-scaling_factors[y,:,:])*100, weights = _harv_area)
            _unitLF_hs = np.ma.average(unitLF[y,:,:][land_suit_all>self.lf_class[1]], weights = self.prod[y,:,:][land_suit_all>self.lf_class[1]])
            _unitLF_ms = np.ma.average(unitLF[y,:,:][(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])], weights = self.prod[y,:,:][(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])]) 
            _unitLF_ps = np.ma.average(unitLF[y,:,:][(land_suit_all<=self.lf_class[0])], weights = self.prod[y,:,:][(land_suit_all<=self.lf_class[0])])
            
            # Suitability
            _harv_area_hs = np.ma.sum(_harv_area[land_suit_all>self.lf_class[1]])/_harv_area_sum*100
            _harv_area_ms = np.ma.sum(_harv_area[(land_suit_all>self.lf_class[0]) & (land_suit_all<=self.lf_class[1])])/_harv_area_sum*100
            _harv_area_ps = np.ma.sum(_harv_area[(land_suit_all>0) & (land_suit_all<=self.lf_class[0])])/_harv_area_sum*100
            _harv_area_wavg = np.ma.average(land_suit_all, weights = _harv_area)
            
            _text = f"{self.first_year+y}\t{_unitLF_hs:.2f}\t{_unitLF_ms:.2f}\t{_unitLF_ps:.2f}\t{_unitLF:.2f}\t{np.ma.sum(self.prod[y,:,:]):.2f}\t{_harv_area_sum:.2f}\t{_harv_area_hs:.2f}\t{_harv_area_ms:.2f}\t{_harv_area_ps:.2f}\t{_harv_area_wavg:.2f}"
            if save_report: ac.printLog(out_report, _text)
            else: print(_text)
        