# -*- coding: utf-8 -*-
"""
Created on Fri May 28 10:00:26 2021

@author: MialykO
"""
#%% Import packages
import numpy as np
import os
import modules.acea.acea_core as ac

#%% Classes
class acea_downscaler:
    "To downscale the results and scale the production to FAOSTAT"
    
    def __init__(self, conf):
        self.conf = conf
        self.data_folder = ac.GetPath(["data","acea"])
        self.raster_folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'rasters'])
        self.gridcells30 = ac.GetAllCells30() # Get gridcells in 30 arc min (id30,rowy30,colx30)
        
        self.clock_start_year = int(self.conf.clock_start.split('/')[0])
        self.clock_end_year = int(self.conf.clock_end.split('/')[0])
        self.max_gs = (self.clock_end_year - self.clock_start_year)+1 - self.conf.spinup; self.first_year = self.clock_start_year + self.conf.spinup
        _,_,self.res5 = ac.GetResolution(1)
        self.text_years = f'{self.first_year}_{self.clock_end_year}'
        
    def Downscale30To5arc(self):
        "Downscale results from 30 to 5 arcminutes"

        variables = ac.GetListOfVaribales() # Get variables
        scenarios = ac.GetListOfScenarios() # Get scenarios
        
        harv_area_rf = ac.GetHarvestedAreas5(self.conf.crop_fao, 0, self.conf.landuse)
        harv_area_ir = ac.GetHarvestedAreas5(self.conf.crop_fao, 1, self.conf.landuse)
        
        consider_shallow_gw = True if 2 in self.conf.scenarios else False
        if consider_shallow_gw: gw_data = ac.ReadNC(ac.GetPath([self.data_folder,'soil','GW_monthly_5arcmin_50m.nc']), 'gw_levels')
        
        for _,_, files in os.walk(self.raster_folder_path):
            for f in files:
                filename = f.split('_')
                
                if '30arc'in f and '.nc' in f and '.xml' not in f:
                    rainfed = True if int(''.join(filter(str.isdigit, filename[2]))) < 3 else False
                    
                    for var, var_info in variables.items():
                        if var in f:
                            print(f'Downscaling {f}')
                                        
                            data5 = np.zeros((self.max_gs,self.res5[0],self.res5[1])); data5 = np.ma.array(data5, mask=(data5==0))
                            new_name = 'acea_5arc_' + '_'.join(filename[2:])
                            data30 = ac.ReadNC(ac.GetPath([self.raster_folder_path, f]), f'{var}-{self.conf.crop_name_short}')
                            
                            cell_count = 0
                            for idx, cell in self.gridcells30.iterrows():
                                _temp_data = data30[:,int(cell['rowy30']), int(cell['colx30'])]
                                
                                if np.any(_temp_data.mask==False):
                                    
                                    [x5, y5] = [int(cell['colx30'])*6, int(cell['rowy30'])*6]
                                    x_range5 = [x5, x5+6]; y_range5 = [y5, y5+6]
                                    
                                    if rainfed:
                                        if consider_shallow_gw:
                                            _mask = np.zeros((6,6)); _mask = np.ma.array(_mask, mask=(_mask==0))
                                            shallow_gw_monthly = np.ma.filled(gw_data[:, y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]],-50) # Groundwater levels at 5 arc mins
                                            if 'sc2' in filename:
                                                _mask[np.max(shallow_gw_monthly, 0) >= self.conf.gw_min_level] = 1
                                                _area5 = harv_area_rf[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]] * _mask
                                            else:
                                                _mask[np.max(shallow_gw_monthly, 0) < self.conf.gw_min_level] = 1
                                                _area5 = harv_area_rf[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]] * _mask
                                        else:
                                            _area5 = harv_area_rf[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
                                    else:
                                        _area5 = harv_area_ir[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
                                        
                                    if np.any(_area5.mask == False):
                                        _mask = np.stack([_area5.mask for _ in range(self.max_gs)], axis=0)
                                        data5[:,y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]] = \
                                            np.ma.array(np.array([_temp_data for _ in range(36)]).reshape(6,6,-1).transpose(), mask=(_mask))
                                    
                                        cell_count += 1
                                
                            # Save raster
                            val_name = f'{var}-{self.conf.crop_name_short}'; val_unit = var_info[1]
                            raster_path = ac.GetPath([self.raster_folder_path, new_name])
                            sc_desc = scenarios[int(''.join(filter(str.isdigit, filename[2])))]
                            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Scenario: {sc_desc}, Variable: {var_info[0]}'
                            
                            ac.CreateHistRaster5(data5, raster_path, val_name, val_unit, title, self.first_year, self.max_gs)
                            print(f'\tCells downscaled: {cell_count}')
    

# class acea_validator:
#     "To validate main ACEA outputs"
    
#     def __init__(self, conf):
#         self.conf = conf
        
#         # Get extra configurations
#         self.clock_start_year = int(self.conf.clock_start.split('/')[0])
#         self.clock_end_year = int(self.conf.clock_end.split('/')[0])
#         self.raster_folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'rasters'])
#         self.max_gs = (self.clock_end_year - self.clock_start_year)-1
#         self.first_year = self.clock_end_year-self.max_gs+1
#         self.text_years = f'{self.first_year}_{self.clock_end_year}'
#         self.data_folder = ac.GetPath(["data","acea"])
        
#     def CompareWithFAOSTAT(self):
#         "Compare national annual production to FAOSTAT"
        
#         out_folder = ac.GetPath(["outputs", f"{self.conf.project_name}",'reports'])
#         if not os.path.exists(out_folder): os.makedirs(out_folder)
#         out_report = ac.GetPath([out_folder, f'{self.conf.crop_fao}_FAOcomparisonReport_{self.text_years}.txt'])
#         ac.printLog(out_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        
#         # Get FAOSTAT production
#         faostat_prod = np.array(pd.read_csv(ac.GetPath([self.data_folder,'faostat',f'FAOSTAT_data_{self.conf.crop_fao}_prod_world.csv']), usecols = ['Value'])['Value'].tolist())
#         faostat_ha = np.array(pd.read_csv(ac.GetPath([self.data_folder,'faostat', f'FAOSTAT_data_{self.conf.crop_fao}_ha_world.csv']), usecols = ['Value'])['Value'].tolist())
        
#         # Get ACEA production
#         harv_area_rf = ac.GetHistHarvestedAreas5(self.conf.crop_fao, False); harv_area_ir = ac.GetHistHarvestedAreas5(self.conf.crop_fao, True)
        
#         _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
#         prod_rf = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Rainfed production
#         _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
#         prod_ir = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Irrigated production
        
#         ac.printLog(out_report, 'Year\t FAO yield\t ACEA yield [t ha-1]\t FAO production\t ACEA production [t]\t FAO area\t ACEA area rf\t ACEA area ir [ha]')
        
#         for y in range(self.max_gs):
#             acea_prod = np.sum(prod_rf[y,:,:]) + np.sum(prod_ir[y,:,:])
#             acea_ha_rf = np.ma.sum(harv_area_rf[y,:,:]); acea_ha_ir = np.ma.sum(harv_area_ir[y,:,:])
#             acea_ha = acea_ha_rf + acea_ha_ir
            
#             ac.printLog(out_report, '{}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(self.first_year+y, \
#                                                                  faostat_prod[y]/faostat_ha[y], acea_prod/acea_ha, faostat_prod[y], acea_prod, faostat_ha[y], acea_ha_rf, acea_ha_ir))
#         ac.printLog(out_report, '================================')
        
    # def TempAllCells(self):
        
    #     scenarios = ac.GetListOfScenarios() # Get scenarios
    #     variables = ac.GetListOfVaribales() # Get variables
        
    #     for sc, sc_desc in scenarios.items():
    #         res = {}
    #         for var, var_info in variables.items(): # Do it dynamically
    #             res_file_path = ac.GetPath([self.raster_folder_path,f'acea_30arc_sc{sc}_{self.conf.crop_name_short}_{var}_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
                
    #             # Check if file exist
    #             if os.path.exists(res_file_path):
    #                 res[var] = ac.ReadNC(res_file_path, f'{var}-{self.conf.crop_name_short}')
                
    #         # Read data
    #         if len(res)>0:
    #             print(f'Results for scenario {sc}')
    #             vars_real = [v for v in res]
                
    #             print('\t\tYear\t' + '\t'.join(vars_real))
    #             for y in range(0,self.years_to_consider):
    #                 _values = []
    #                 for var in vars_real:
    #                     _values.append(f'{np.around(np.ma.average(res[var][y,:,:]),2)}')
    #                 print(f'\t\t{self.clock_start_year +2 + y}\t' + '\t'.join(_values))
