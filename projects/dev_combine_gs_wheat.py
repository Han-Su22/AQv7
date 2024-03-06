# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 15:34:36 2022

@author: HanS, MialykO
"""
import modules.acea.acea_core as ac
import os
import numpy as np

project1 = 'SpringWheat'
project2 = 'WinterWheat'
gs_mask_ = r"D:\ACEA_v2_2\data\acea\harvested_areas\15\winter_and_spring_wheat_areas_phase3.nc4"
climate_name = 'gswp3-w5e5_obsclim/gswp3-w5e5_obsclim_01-01-1981_31-12-2019'
crop_name = 'Wheat'
gs_to_consider = 5; first_year = 2008
new_crop_name = 'wht'

gs_mask =  ac.ReadNC(gs_mask_, 'swh_mask')
gs_mask = np.ones((gs_to_consider, 360, 720)) * gs_mask


new_folder = r'D:\ACEA_v2_2\outputs\Wheat'


irrscenarios=['Highinput','HighinputNoSF','Highvirt','Lowinput','Lowvirt']

for sci in irrscenarios:

    new_path = ac.GetPath([new_folder, sci,'rasters'])
    if not os.path.exists(new_path): os.makedirs(new_path)
    folder_path_gs1 = ac.GetPath(["outputs", project1,sci,'rasters'])
    folder_path_gs2 = ac.GetPath(["outputs", project2,sci,'rasters'])
    
    scenarios = ac.GetListOfScenarios() # Get scenarios
    variables = ac.GetListOfVaribales() # Get variables
    for _,_, files in os.walk(folder_path_gs1):
        for f in files:
            filename = f.split('_')
            
            if filename[-1].split('.')[-1]!='nc':
                continue
            
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
                    raster_path = ac.GetPath([new_path, new_name])
                    title = f'Crop: {crop_name}, Climate: {climate_name}, Scenario: {sc_desc}, Variable: {var_info[0]}'
                    
                    ac.CreateHistRaster30(data30, raster_path, val_name, val_unit, title, first_year, gs_to_consider)


new_path = ac.GetPath([new_folder, 'rastersSF'])
if not os.path.exists(new_path): os.makedirs(new_path)
folder_path_gs1 = ac.GetPath(["outputs", project1,'rastersSF'])
folder_path_gs2 = ac.GetPath(["outputs", project2,'rastersSF'])


gs_mask =  ac.ReadNC(gs_mask_, 'swh_mask')
gs_mask = np.ones((360, 720)) * gs_mask

for _,_, files in os.walk(folder_path_gs1):
    for f in files:
        filename = f.split('_')
        
        if filename[-1].split('.')[-1]!='nc':
            continue
        
        if len(filename)==4:
            varname=filename[-1].split('.')[0]
        else:
            varname='cur_'+filename[-1].split('.')[0]
        
        data30 = ac.ReadNC(ac.GetPath([folder_path_gs1, f]), varname)
        f2 = 'WinterWheat_'+'_'.join(filename[1:])
        data30_gs2 = ac.ReadNC(ac.GetPath([folder_path_gs2, f2]), varname)
        
        data30[gs_mask==0] = data30_gs2[gs_mask==0]


        new_name = crop_name+'_'+filename[1] + '_' + filename[2] +'_' + varname+'.nc'

        # Save
        val_name = varname; val_unit = '1'
        raster_path = ac.GetPath([new_path, new_name])
        title = varname
        
        ac.CreateRaster30(data30, raster_path, val_name, val_unit, title)

        



                    