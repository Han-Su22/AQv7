# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 15:34:36 2022

@author: MialykO
"""
import modules.acea.acea_core as ac
import os
import numpy as np

project1 = 'ri1_phd_mialyk_2022'
project2 = 'ri2_phd_mialyk_2022'
gs_mask = r"O:\Oleks\ACEA_v2\data\acea\crop_calendar\ri1_ggcmi_ir_crop_calendar.nc"
climate_name = 'gswp3-w5e5_obsclim/gswp3-w5e5_obsclim_01-01-1981_31-12-2019'
crop_name = 'Rice'
gs_to_consider = 30; first_year = 1990
new_crop_name = 'ri'
new_folder = r"O:\Oleks\ACEA_v2\outputs\ri_phd_mialyk_2022\rasters"
folder_path_gs1 = ac.GetPath(["outputs", project1,'rasters'])
folder_path_gs2 = ac.GetPath(["outputs", project2,'rasters'])


gs_mask =  ac.ReadNC(gs_mask, 'fraction_of_harvested_area')
# gs_mask2 =  ac.ReadNC(r"O:\Oleks\ACEA_v2\data\acea\crop_calendar\ri1_ggcmi_rf_crop_calendar.nc", 'fraction_of_harvested_area')
# gs_mask3 =  ac.ReadNC(r"O:\Oleks\ACEA_v2\data\acea\crop_calendar\ri2_ggcmi_ir_crop_calendar.nc", 'fraction_of_harvested_area')
# gs_mask4 =  ac.ReadNC(r"O:\Oleks\ACEA_v2\data\acea\crop_calendar\ri2_ggcmi_rf_crop_calendar.nc", 'fraction_of_harvested_area')

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
                data30_gs1 = ac.ReadNC(ac.GetPath([folder_path_gs1, f]), f'{var}-{crop_name_short}')
                crop_name_short2 = 'ri2'; f2 = '_'.join(filename[:3]) + '_' + crop_name_short2 + '_' + '_'.join(filename[4:])
                data30_gs2 = ac.ReadNC(ac.GetPath([folder_path_gs2, f2]), f'{var}-{crop_name_short2}')
                
                # Procedure
                data30_1 = np.ones((30, 360, 720)) * np.ma.masked
                data30_2 = np.ones((30, 360, 720)) * np.ma.masked
                
                data30_1[gs_mask>0] = data30_gs1[gs_mask>0] * gs_mask[gs_mask>0]
                gs_mask_rev = 1- gs_mask
                data30_2[gs_mask>0] = data30_gs2[gs_mask>0] * gs_mask_rev[gs_mask>0]
                
                data30 = np.ma.filled(data30_1,0) + np.ma.filled(data30_2,0)
                
                data30 = np.ma.masked_where(data30==0, data30)
                
                
                # Save
                sc_desc = scenarios[int(''.join(filter(str.isdigit, filename[2])))]
                val_name = f'{var}-{new_crop_name}'; val_unit = var_info[1]
                raster_path = ac.GetPath([new_folder, new_name])
                title = f'Crop: {crop_name}, Climate: {climate_name}, Scenario: {sc_desc}, Variable: {var_info[0]}'
                
                ac.CreateHistRaster30(data30, raster_path, val_name, val_unit, title, first_year, gs_to_consider)
                