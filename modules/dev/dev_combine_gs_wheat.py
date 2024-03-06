# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 15:34:36 2022

@author: MialykO
"""
import modules.acea.acea_core as ac
import os
import numpy as np

project1 = 'swh_phd_mialyk_2022'
project2 = 'wwh_phd_mialyk_2022'
gs_mask = r"O:\Oleks\ACEA_v2\data\acea\harvested_areas\15\winter_and_spring_wheat_areas_phase3.nc4"
climate_name = 'gswp3-w5e5_obsclim/gswp3-w5e5_obsclim_01-01-1981_31-12-2019'
crop_name = 'Wheat'
gs_to_consider = 30; first_year = 1990
new_crop_name = 'wh'
new_folder = r'O:\Oleks\ACEA_v2\outputs\wh_phd_mialyk_2022\rasters'
folder_path_gs1 = ac.GetPath(["outputs", project1,'rasters'])
folder_path_gs2 = ac.GetPath(["outputs", project2,'rasters'])


gs_mask =  ac.ReadNC(gs_mask, 'swh_mask')
gs_mask = np.ones((30, 360, 720)) * gs_mask

scenarios = ac.GetListOfScenarios() # Get scenarios
variables = ac.GetListOfVaribales() # Get variables
for _,_, files in os.walk(folder_path_gs1):
    for f in files:
        filename = f.split('_')
        sc_desc = filename[2]
        crop_name_short = filename[3]
        for var, var_info in variables.items():
            if var in f:
                new_name = '_'.join(filename[:3]) + '_' + new_crop_name +'_' + '_'.join(filename[4:])
                data30 = ac.ReadNC(ac.GetPath([folder_path_gs1, f]), f'{var}-{crop_name_short}')
                crop_name_short2 = 'wwh'; f2 = '_'.join(filename[:3]) + '_' + crop_name_short2 + '_' + '_'.join(filename[4:])
                data30_gs2 = ac.ReadNC(ac.GetPath([folder_path_gs2, f2]), f'{var}-{crop_name_short2}')
                
                data30[gs_mask==0] = data30_gs2[gs_mask==0] 
                    
                # Save
                sc_desc = scenarios[int(''.join(filter(str.isdigit, filename[2])))]
                val_name = f'{var}-{new_crop_name}'; val_unit = var_info[1]
                raster_path = ac.GetPath([new_folder, new_name])
                title = f'Crop: {crop_name}, Climate: {climate_name}, Scenario: {sc_desc}, Variable: {var_info[0]}'
                
                ac.CreateHistRaster30(data30, raster_path, val_name, val_unit, title, first_year, gs_to_consider)
                