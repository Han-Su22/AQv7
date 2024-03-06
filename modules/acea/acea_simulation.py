# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 13:47:06 2021

@author: MialykO
"""
#%% Import packages
import pickle
import os
import datetime
import numpy as np
import traceback
import time
import pandas as pd
from functools import partial
from osgeo import gdal

import modules.acea.acea_core as ac
from modules.acea.acea_crops import CropParameters

os.environ["DEVELOPMENT"] = '1'

from aquacrop.core import AquaCropModel
from aquacrop.entities.soil import Soil as SoilClass
from aquacrop.entities.inititalWaterContent import InitialWaterContent as InitWCClass
from aquacrop.entities.irrigationManagement import IrrigationManagement as IrrMngtClass
from aquacrop.entities.groundWater import GroundWater as GwClass
from aquacrop.entities.fieldManagement import FieldMngt as FieldMngtClass
from aquacrop.entities.co2 import CO2


#%% Simulation
def RunACEA(model, cell_id):
    'Call AquaCrop for a cell'
    
    #from modules.aquacrop_v6_02.core import AquaCropModel
    #from modules.aquacrop_v6_02.classes import SoilClass,InitWCClass,IrrMngtClass,GwClass, FieldMngtClass
    #from modules.acea.acea_crops import CropParameters
    
    # Current cell
    cell = model.gridcells30[model.gridcells30['id30']==cell_id] 
    
    # 1. Get weather
    with open(ac.GetPath(["data", "acea",'climate',f'{model.conf.climate_name}_{cell_id}.pckl']), 'rb') as f:
        _tmax, _tmin, _prec, _et0 = pickle.load(f); _et0[_et0==0] = .001
        data = {'MinTemp': _tmin,
                'MaxTemp': _tmax,
                'Precipitation': _prec,
                'ReferenceET': _et0,
                'Date': pd.date_range(start=np.datetime64(str(model.conf.climate_start)+'-01-01'),\
                                      end=np.datetime64(str(model.conf.climate_end)+'-12-31'))}
        wdf = pd.DataFrame(data)
        
    # 2. Get soil
    init_wc = InitWCClass(wc_type='Pct',value=[model.conf.init_wc]) # Define initial soil water conditions
    
    use_soilgrids=model.conf.use_soilgrids
    if use_soilgrids==0:
        soil_data = SoilClass('custom', dz = model.conf.soil_dz,calc_cn=1,adj_rew=0) # define soil
        soil_texture = model.GetSoilTexture(cell,use_soilgrids=use_soilgrids)
        soil_data.add_layer_from_texture(thickness=sum(model.conf.soil_dz),
                                              Sand=soil_texture[0],Clay=soil_texture[1],
                                              OrgMat=2.5,penetrability=100)
    else:
        soil_data = SoilClass('custom', dz = model.conf.soil_dz,calc_cn=1,adj_rew=0) # define soil
        soil_texture = model.GetSoilTexture(cell,use_soilgrids=use_soilgrids)
        #'0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm'
        #[.1, .1, .1, .3, .4, .6, .7, .7], be cautious when dz changes
        for i in range(len(model.conf.soil_dz)):
            if i<2:
                soil_data.add_layer_from_texture(thickness=model.conf.soil_dz[i],
                                                      Sand=(soil_texture[0][i]+soil_texture[0][i+1])/2,\
                                                          Clay=(soil_texture[1][i]+soil_texture[1][i+1])/2,\
                                                      OrgMat=(soil_texture[2][i]+soil_texture[2][i+1])/2,\
                                                          penetrability=(soil_texture[3][i]+soil_texture[3][i+1])/2)
            elif i<5:
                soil_data.add_layer_from_texture(thickness=model.conf.soil_dz[i],
                                                      Sand=soil_texture[0][i],\
                                                          Clay=soil_texture[1][i],
                                                      OrgMat=soil_texture[2][i],penetrability=soil_texture[3][i])
            else:
                soil_data.add_layer_from_texture(thickness=model.conf.soil_dz[i],
                                                      Sand=soil_texture[0][5],\
                                                          Clay=soil_texture[1][5],
                                                      OrgMat=soil_texture[2][5],penetrability=soil_texture[3][5])
                
    #print(soil_texture)
    
    #print(soil_data.profile[['th_wp','th_fc','th_s','Ksat']])

    # 3. Get scenarios
    scenarios = model.AdjustScenarios(cell) # Check which scenarios to run
    
    # 4. Run scenarios
    sim_success = False
    if len(scenarios) > 0:
        if not model.conf.multi_core: print(f'> Cell {cell_id}:')
        for sc in scenarios:
            try:
                gw = GwClass(water_table='N') # Define initial groundwater conditions
                
                # Get irrigation management
                if sc < 3: # Rainfed
                    irmngt = IrrMngtClass(irrigation_method=0); irrigated = 0
                    if model.conf.virtual_irrigation in ['Lowvirt','Highvirt','HighvirtNoSF']:
                        max_ir = min(500, int(soil_data.profile['Ksat'][0]))
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 100, MaxIrr=max_ir)
                    
                    if sc == 2: # Rainfed with gw
                        gw_dates, gw_values = model.GetGroundwaterLevels(cell)
                        if len(gw_values)>0: gw = GwClass(water_table='Y', method = "Variable", dates=gw_dates, values=gw_values)
                        else: continue # If no gw levels, then skip this scenario
                else: # Irrigated
                    irrigated = 1
                    max_ir = min(500, int(soil_data.profile['Ksat'][0])) # Limit irrigation to Ksat of the top layer to avoid runoff
                    
                    if sc == 3: # Irrigated with furrow irrigation, no gw
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 80, MaxIrr=max_ir)
                    elif sc == 4: # Irrigated with sprinklers, no gw
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 100, MaxIrr=max_ir)
                    elif sc == 5: # Irrigated with drip, no gw
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 30, MaxIrr=max_ir)
                    elif sc == 6: # Irrigated with flood, with gw (only for rice)
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 100, MaxIrr=max_ir)
                        gw_dates, gw_values = model.GetGroundwaterLevels(cell, False)
                        if len(gw_values)>0: gw = GwClass(water_table='Y', method = "Variable", dates=gw_dates, values=gw_values)
                        
                    else:
                         if not model.conf.multi_core: print('Error simulation scenario does not exist')
                         return
                
                if model.conf.soil_fertility==1:
                    if model.conf.tuned==1:
                        try:
                            ctry=model.gridcells30ctry.loc[cell['id30'],'faostat'].values
                            yield_scaling=0
                            i_=0
                            for c in ctry:
                                if model.scaling.loc[c,'scaling']<np.inf:
                                    yield_scaling+=model.scaling.loc[c,'scaling']
                                    i_+=1
                            if i_>0:
                                yield_scaling=yield_scaling/i_
                            else:
                                yield_scaling=1
                        except:
                            yield_scaling=1
                    else:
                        yield_scaling=1

                    folder_path = ac.GetPath(["outputs", f"{model.conf.project_name}",'Para_SF'])
                    if sc<3:
                        file_name = ac.GetPath([folder_path, f'{model.conf.project_name}_{model.conf.crop_fao}_{cell_id}_rainfed.npz'])
                    else:
                        file_name = ac.GetPath([folder_path, f'{model.conf.project_name}_{model.conf.crop_fao}_{cell_id}_irr.npz'])
                        
                    try:
                        data = np.load(file_name)
                        cur_search=data['cur_search']
                        if cur_search>-1:
                            
                            sf_es=data['sf_es']#in some extream cases, sf could be negative, need to be fixed in next version
                            Ksexpf_es=data['Ksexpf_es']
                            fcdecline_es=data['fcdecline_es']
                            Kswp_es=data['Kswp_es']
                            Ksccx_es=data['Ksccx_es']
                            relbio_es=data['relbio_es']

                            if model.conf.virtual_irrigation=='Lowvirt' or model.conf.virtual_irrigation=='Lowinput':
                                #if low input, use the calibrated value as defult

                                relbio=relbio_es[0]
                                
                            else:
                                #use certain soil fertility stress for high input from GAEZ
                                
                                relbio,sea_=model.GetfromGAEZ(cell, irrigated,'relbio_high')
                                if sea_==-1:
                                    relbio=1
                                    
                                #ensure irrigation have higher biomass
                                if relbio<relbio_es[0]:
                                    relbio=relbio_es[0]

                            relbio=relbio*yield_scaling
                            
                            if relbio>1:
                                relbio=1
                                
                            loc_=np.argmin(np.abs(relbio_es[0:100]-relbio))
                            
                            Ksccx=Ksccx_es[loc_]
                            Ksexpf=Ksexpf_es[loc_]
                            Kswp=Kswp_es[loc_]
                            fcdecline=fcdecline_es[loc_]
                            sfertstress=sf_es[loc_]
                            
                            if sfertstress>0.8:
                                loc_=np.argmin(np.abs(sf_es[0:100]-0.8))
                                
                                Ksccx=Ksccx_es[loc_]
                                Ksexpf=Ksexpf_es[loc_]
                                Kswp=Kswp_es[loc_]
                                fcdecline=fcdecline_es[loc_]
                                sfertstress=sf_es[loc_]

                        else:
                            if not model.conf.multi_core: print('No calibrated parameters!')
                            return
                    except:
                        print('Run soil fertility calibration first!')

                else:
                    Ksccx=1
                    Ksexpf=1
                    Kswp=1
                    fcdecline=0
                    sfertstress=0
                    Ksccx_es=np.zeros(10000)+1
                    Ksexpf_es=np.zeros(10000)+1
                    Kswp_es=np.zeros(10000)+1
                    fcdecline_es=np.zeros(10000)
                    sf_es=np.zeros(10000)
                    relbio_es=np.zeros(10000)+1

                #not used, in case will use for future update
                if model.conf.virtual_irrigation=='Lowvirt' or model.conf.virtual_irrigation=='Lowinput':
                    HI0,search_hi=model.GetfromGAEZ(cell, irrigated,'hi_low')
                else:
                    HI0,search_hi=model.GetfromGAEZ(cell, irrigated,'hi_high')
                CCx,search_ccx=model.GetfromGAEZ(cell, irrigated,'cc_high')
                if model.conf.get_GAEZ_ccx_hi==False:
                    search_hi=-1
                    search_ccx=-1

                #'HI0','CCx'    
                # Get crop
                planting_date, harvest_date = model.GetCropCycleDates(cell, irrigated)
                
                if search_hi>-1 and search_ccx>-1:
                    crop_data = CropParameters(model.conf, cell, irrigated, model.crop_params_rf, model.crop_params_ir,\
                                                planting_date, harvest_date=harvest_date, \
                                                PlantPop=0,\
                                                Ksccx=Ksccx,\
                                                Ksexpf=Ksexpf,\
                                                Kswp=Kswp,\
                                                fcdecline=fcdecline,\
                                                sfertstress=sfertstress,sf_es=sf_es,Ksexpf_es=Ksexpf_es,fcdecline_es=fcdecline_es,
                                                Kswp_es=Kswp_es,Ksccx_es=Ksccx_es,relbio_es=relbio_es,HI0=HI0,CCx=CCx)
                elif search_hi>-1:
                    crop_data = CropParameters(model.conf, cell, irrigated, model.crop_params_rf, model.crop_params_ir,\
                                                planting_date, harvest_date=harvest_date, \
                                                PlantPop=0,\
                                                Ksccx=Ksccx,\
                                                Ksexpf=Ksexpf,\
                                                Kswp=Kswp,\
                                                fcdecline=fcdecline,\
                                                sfertstress=sfertstress,sf_es=sf_es,Ksexpf_es=Ksexpf_es,fcdecline_es=fcdecline_es,
                                                Kswp_es=Kswp_es,Ksccx_es=Ksccx_es,relbio_es=relbio_es,HI0=HI0)
                elif search_ccx>-1:
                    crop_data = CropParameters(model.conf, cell, irrigated, model.crop_params_rf, model.crop_params_ir,\
                                                planting_date, harvest_date=harvest_date, \
                                                PlantPop=0,\
                                                Ksccx=Ksccx,\
                                                Ksexpf=Ksexpf,\
                                                Kswp=Kswp,\
                                                fcdecline=fcdecline,\
                                                sfertstress=sfertstress,sf_es=sf_es,Ksexpf_es=Ksexpf_es,fcdecline_es=fcdecline_es,
                                                Kswp_es=Kswp_es,Ksccx_es=Ksccx_es,relbio_es=relbio_es,CCx=CCx)
                else:
                    crop_data = CropParameters(model.conf, cell, irrigated, model.crop_params_rf, model.crop_params_ir,\
                                                planting_date, harvest_date=harvest_date, \
                                                PlantPop=0,\
                                                Ksccx=Ksccx,\
                                                Ksexpf=Ksexpf,\
                                                Kswp=Kswp,\
                                                fcdecline=fcdecline,\
                                                sfertstress=sfertstress,sf_es=sf_es,Ksexpf_es=Ksexpf_es,fcdecline_es=fcdecline_es,
                                                Kswp_es=Kswp_es,Ksccx_es=Ksccx_es,relbio_es=relbio_es)
                
                
                if model.conf.crop_phenology in ['transient','average']: crop_data.calibrate_phenology(cell,model.conf)   

                _co2_file = ac.GetPath(["data", "acea", 'co2', f'{model.conf.co2_name}.txt'])
                _co2_data = pd.read_csv(_co2_file,header=1,delim_whitespace=True,names=["year", "ppm"])
                _co2 = CO2(co2_data=_co2_data)
                
                # Get field management
                fldmngt = FieldMngtClass(mulches=model.conf.mulching,bunds=model.conf.bunds,
                                         mulch_pct=model.conf.mulching_area*100,f_mulch=model.conf.mulching_factor,
                                         z_bund=model.conf.bunds_dz)
                
                fal_fldmngt = FieldMngtClass(mulches=False, bunds=model.conf.bunds, z_bund=model.conf.bunds_dz)
                # Create and run model

                AOS_Py = AquaCropModel(sim_start_time=f'{model.conf.clock_start}',
                                       sim_end_time=f'{model.conf.clock_end}',
                                       weather_df=wdf, 
                                       soil=soil_data,
                                       crop=crop_data,
                                       irrigation_management=irmngt, 
                                       field_management=fldmngt, 
                                       fallow_field_management=fal_fldmngt, 
                                       initial_water_content=init_wc, 
                                       groundwater=gw, 
                                       off_season=model.conf.off_season,
                                       co2_concentration=_co2,)
                #AOS_Py._initialize()
                AOS_Py.run_model(initialize_model=True, till_termination=True)
            
                if not model.conf.multi_core: model.DisplaySummary(AOS_Py._outputs, sc)
                model.SaveRawResults(AOS_Py._outputs, sc, cell); sim_success = True
            except Exception:
                if not model.conf.multi_core: print (traceback.print_exc())
                continue
                
        if sim_success: return cell_id
    return


def InitACEA(model, cell_id):
    'Call AquaCrop for a cell'
    
    #from modules.aquacrop_v6_02.core import AquaCropModel
    #from modules.aquacrop_v6_02.classes import SoilClass,InitWCClass,IrrMngtClass,GwClass, FieldMngtClass
    #from modules.acea.acea_crops import CropParameters
    
    # Current cell
    cell = model.gridcells30[model.gridcells30['id30']==cell_id] 
    
    # 1. Get weather
    with open(ac.GetPath(["data", "acea",'climate',f'{model.conf.climate_name}_{cell_id}.pckl']), 'rb') as f:
        _tmax, _tmin, _prec, _et0 = pickle.load(f); _et0[_et0==0] = .001
        data = {'MinTemp': _tmin,
                'MaxTemp': _tmax,
                'Precipitation': _prec,
                'ReferenceET': _et0,
                'Date': pd.date_range(start=np.datetime64(str(model.conf.climate_start)+'-01-01'),\
                                      end=np.datetime64(str(model.conf.climate_end)+'-12-31'))}
        wdf = pd.DataFrame(data)

    # 2. Get soil
    init_wc = InitWCClass(wc_type='Pct',value=[model.conf.init_wc]) # Define initial soil water conditions
    #print(cell)
    use_soilgrids=model.conf.use_soilgrids
    if use_soilgrids==0:
        soil_data = SoilClass('custom', dz = model.conf.soil_dz,calc_cn=1,adj_rew=0) # define soil
        soil_texture = model.GetSoilTexture(cell,use_soilgrids=use_soilgrids)
        soil_data.add_layer_from_texture(thickness=sum(model.conf.soil_dz),
                                              Sand=soil_texture[0],Clay=soil_texture[1],
                                              OrgMat=2.5,penetrability=100)
    else:
        soil_data = SoilClass('custom', dz = model.conf.soil_dz,calc_cn=1,adj_rew=0) # define soil
        soil_texture = model.GetSoilTexture(cell,use_soilgrids=use_soilgrids)
        #'0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm'
        #[.1, .1, .1, .3, .4, .6, .7, .7], be cautious when dz changes
        for i in range(len(model.conf.soil_dz)):
            if i<2:
                soil_data.add_layer_from_texture(thickness=model.conf.soil_dz[i],
                                                      Sand=(soil_texture[0][i]+soil_texture[0][i+1])/2,\
                                                          Clay=(soil_texture[1][i]+soil_texture[1][i+1])/2,\
                                                      OrgMat=(soil_texture[2][i]+soil_texture[2][i+1])/2,\
                                                          penetrability=(soil_texture[3][i]+soil_texture[3][i+1])/2)
            elif i<5:
                soil_data.add_layer_from_texture(thickness=model.conf.soil_dz[i],
                                                      Sand=soil_texture[0][i],\
                                                          Clay=soil_texture[1][i],
                                                      OrgMat=soil_texture[2][i],penetrability=soil_texture[3][i])
            else:
                soil_data.add_layer_from_texture(thickness=model.conf.soil_dz[i],
                                                      Sand=soil_texture[0][5],\
                                                          Clay=soil_texture[1][5],
                                                      OrgMat=soil_texture[2][5],penetrability=soil_texture[3][5])

    # 3. Get scenarios
    scenarios = model.AdjustScenarios(cell) # Check which scenarios to run

    # 4. Run scenarios
    sim_success = False
    if len(scenarios) > 0:
        if not model.conf.multi_core: print(f'> Cell {cell_id}:')
        for sc in scenarios:
            try:
                gw = GwClass(water_table='N') # Define initial groundwater conditions
                # Get irrigation management
                if sc < 3: # Rainfed
                    irmngt = IrrMngtClass(irrigation_method=0); irrigated = 0
                    if sc == 2: # Rainfed with gw
                        gw_dates, gw_values = model.GetGroundwaterLevels(cell)
                        if len(gw_values)>0: gw = GwClass(water_table='Y', method = "Variable", dates=gw_dates, values=gw_values)
                        else: continue # If no gw levels, then skip this scenario
                else: # Irrigated
                    irrigated = 1
                    max_ir = min(500, int(soil_data.profile['Ksat'][0])) # Limit irrigation to Ksat of the top layer to avoid runoff
                    
                    if sc == 3: # Irrigated with furrow irrigation, no gw
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 80, MaxIrr=max_ir)
                    elif sc == 4: # Irrigated with sprinklers, no gw
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 100, MaxIrr=max_ir)
                    elif sc == 5: # Irrigated with drip, no gw
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 30, MaxIrr=max_ir)
                    elif sc == 6: # Irrigated with flood, with gw (only for rice)
                        irmngt = IrrMngtClass(irrigation_method=1, SMT = model.conf.irr_thresholds, WetSurf = 100, MaxIrr=max_ir)
                        gw_dates, gw_values = model.GetGroundwaterLevels(cell, False)
                        if len(gw_values)>0: gw = GwClass(water_table='Y', method = "Variable", dates=gw_dates, values=gw_values)
                        
                    else:
                         if not model.conf.multi_core: print('Error simulation scenario does not exist')
                         return
                if sc==1:
                    RelativeBio,Ksccx_in,f3low,f3high,cur_search=model.GetBiomassCanopySF(cell,rainfed=True)
                else:
                    RelativeBio,Ksccx_in,f3low,f3high,cur_search=model.GetBiomassCanopySF(cell,rainfed=False)

                #if no data for calibration, skip calibration
                if RelativeBio==-1:
                    if sc==1:
                        model.SaveParametersSF([],f3low,f3high,cur_search, cell,'rainfed')
                    else:
                        model.SaveParametersSF([],f3low,f3high,cur_search, cell,'irr')
                    if not model.conf.multi_core: print('No GAEZ data available')
                    continue
                
                #only for rubber, GAEZ have no rubber for soil fertility calibration
                if model.conf.project_name == "Rubber":
                    RelativeBio=1

                CCx,search_ccx=model.GetfromGAEZ(cell, irrigated,'cc_high')
                
                if model.conf.get_GAEZ_ccx_hi==False:
                    search_ccx=-1
                
                #'HI0','CCx'    
                # Get crop
                planting_date, harvest_date = model.GetCropCycleDates(cell, irrigated)
                
                if search_ccx>-1:
                    crop_data = CropParameters(model.conf, cell, irrigated, model.crop_params_rf, model.crop_params_ir,\
                                            planting_date, harvest_date=harvest_date, \
                                            PlantPop=0,\
                                            need_calib=1,\
                                            RelativeBio=RelativeBio,\
                                            Ksccx_in=Ksccx_in,\
                                            fcdecline_in=1,CCx=CCx)
                else:
                    crop_data = CropParameters(model.conf, cell, irrigated, model.crop_params_rf, model.crop_params_ir,\
                                            planting_date, harvest_date=harvest_date, \
                                            PlantPop=0,\
                                            need_calib=1,\
                                            RelativeBio=RelativeBio,\
                                            Ksccx_in=Ksccx_in,\
                                            fcdecline_in=1)
                
                if model.conf.crop_phenology in ['transient','average']: crop_data.calibrate_phenology(cell,model.conf)   
                
                _co2_file = ac.GetPath(["data", "acea", 'co2', f'{model.conf.co2_name}.txt'])
                _co2_data = pd.read_csv(_co2_file,header=1,delim_whitespace=True,names=["year", "ppm"])
                _co2 = CO2(co2_data=_co2_data)

                # Get field management
                fldmngt = FieldMngtClass(mulches=model.conf.mulching,bunds=model.conf.bunds,
                                         mulch_pct=model.conf.mulching_area*100,f_mulch=model.conf.mulching_factor,
                                         z_bund=model.conf.bunds_dz)
                
                fal_fldmngt = FieldMngtClass(mulches=False, bunds=model.conf.bunds, z_bund=model.conf.bunds_dz)
                # Create and run model
                AOS_Py = AquaCropModel(sim_start_time=r'2010/'+planting_date,
                                       sim_end_time=f'{model.conf.clock_end}',
                                       weather_df=wdf, 
                                       soil=soil_data,
                                       crop=crop_data,
                                       irrigation_management=irmngt, 
                                       field_management=fldmngt, 
                                       fallow_field_management=fal_fldmngt, 
                                       initial_water_content=init_wc, 
                                       groundwater=gw, 
                                       off_season=model.conf.off_season,
                                       co2_concentration=_co2,)
                AOS_Py._initialize()
                if sc==1:
                    model.SaveParametersSF(AOS_Py,f3low,f3high,cur_search, cell,'rainfed')
                else:
                    model.SaveParametersSF(AOS_Py,f3low,f3high,cur_search, cell,'irr')
                sim_success = True
                
            except Exception:
                if not model.conf.multi_core: print (traceback.print_exc())
                
                try:
                    if sc==1:
                        model.SaveParametersSF([],f3low,f3high,-1, cell,'rainfed')
                    else:
                        model.SaveParametersSF([],f3low,f3high,-1, cell,'irr')
                except:
                    if sc==1:
                        model.SaveParametersSF([],-1,-1,-1, cell,'rainfed')
                    else:
                        model.SaveParametersSF([],-1,-1,-1, cell,'irr')

                continue
                
        if sim_success: return cell_id
    return


class acea_model:
    # Folder paths
    data_folder = ac.GetPath(["data", "acea"])

    def __init__(self, conf):
        print('========= Welcome to ACEA v2.2 =========')
        self.conf = conf
        self.tic = time.time() # Timer for performance assessment
        
        # Get extra configurations
        self.clock_start_year = int(self.conf.clock_start.split('/')[0])
        self.clock_end_year = int(self.conf.clock_end.split('/')[0])
        
        print(f'Crop: {self.conf.crop_name}')
        print(f'Area: {self.conf.gridcells}')
        print(f'Time period: {self.clock_start_year}-{self.clock_end_year}')
        print(f'Climate: {self.conf.climate_name}')
        print(f'CO2: {self.conf.co2_name}')
        print('\t')
        self.GetInputData()
        
#%% Core functions
    def GetInputData(self):
        print('Getting input data...')
        
        # Get grid
        self.gridcells30 = ac.GetAllCells30() # Get gridcells in 30 arc min (id30,rowy30,colx30)
        #if min(self.conf.scenarios) < 3: # Check if rainfed scenarios are present
        #    self.harvested_areas_rain5 = ac.GetHarvestedAreas5(self.conf.crop_fao, 0, self.conf.landuse)
        #if max(self.conf.scenarios) > 2: # Check if irrigation scenarios are present
        #    self.harvested_areas_ir5 = ac.GetHarvestedAreas5(self.conf.crop_fao, 1, self.conf.landuse)
        self.harvested_areas_rain5 = ac.GetHarvestedAreas5(self.conf.crop_fao, 0, self.conf.landuse)
        self.harvested_areas_ir5 = ac.GetHarvestedAreas5(self.conf.crop_fao, 1, self.conf.landuse)
        
        # Get crop calendar
        self.crop_planting_rf = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_rf_crop_calendar.nc']), 'planting_day')
        self.crop_planting_ir = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_ir_crop_calendar.nc']), 'planting_day')
        self.crop_harvest_rf = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_rf_crop_calendar.nc']), 'maturity_day')
        self.crop_harvest_ir = ac.ReadNC(ac.GetPath([self.data_folder,'crop_calendar',f'{self.conf.crop_name_short}_ir_crop_calendar.nc']), 'maturity_day')
        
        # Get crop params
        parameters = ['gdd_emergence','gdd_max_root','gdd_senescence','gdd_maturity','gdd_yield_form', \
                      'gdd_duration_flowering','gdd_duration_yield_form','cdc','cgc','temp_max_avg','temp_min_avg']
        self.crop_params_rf = dict(); self.crop_params_ir = dict()
        for par in parameters:
            self.crop_params_rf[par] = ac.ReadNC(ac.GetPath([self.data_folder,'crop_parameters',f'{self.conf.crop_name_short}_rf_crop_parameters.nc']), par)
            self.crop_params_ir[par] = ac.ReadNC(ac.GetPath([self.data_folder,'crop_parameters',f'{self.conf.crop_name_short}_ir_crop_parameters.nc']), par)
        
        # Get soil texture
        self.clay_data = ac.ReadNC(ac.GetPath([self.data_folder,'soil','HWSD_soil_data_on_cropland_v2.3.nc']), 'clay')
        self.sand_data = ac.ReadNC(ac.GetPath([self.data_folder,'soil','HWSD_soil_data_on_cropland_v2.3.nc']), 'sand')
        if 2 in self.conf.scenarios or 6 in self.conf.scenarios:
            self.gw_data = ac.ReadNC(ac.GetPath([self.data_folder,'soil','GW_monthly_5arcmin_50m.nc']), 'gw_levels')
            
        
        self.gaez_lut_high=[]
        self.gaez_lut_low=[]
        self.gaez_yield_high=[]
        self.gaez_yield_low=[]
        self.gaez_f3_high=[]
        self.gaez_f3_low=[]
        
        for cultivar in self.conf.crop_gaez:
            ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','LUT',f'LUT_{cultivar}_irr_High_1981_2010.tif']), gdal.GA_ReadOnly)
            rb = ds.GetRasterBand(1)
            self.gaez_lut_high.append(rb.ReadAsArray())
            
            ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','LUT',f'LUT_{cultivar}_irr_Low_1981_2010.tif']), gdal.GA_ReadOnly)
            rb = ds.GetRasterBand(1)
            self.gaez_lut_low.append(rb.ReadAsArray())

            ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','PotentialYield',f'PotentialYield_{cultivar}_irr_High_1981_2010.tif']), gdal.GA_ReadOnly)
            rb = ds.GetRasterBand(1)
            self.gaez_yield_high.append(rb.ReadAsArray())
            
            ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','PotentialYield',f'PotentialYield_{cultivar}_irr_Low_1981_2010.tif']), gdal.GA_ReadOnly)
            rb = ds.GetRasterBand(1)
            self.gaez_yield_low.append(rb.ReadAsArray())
            
            ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','F3',f'F3_{cultivar}_irr_High_1981_2010.tif']), gdal.GA_ReadOnly)
            rb = ds.GetRasterBand(1)
            self.gaez_f3_high.append(rb.ReadAsArray())
            
            ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','F3',f'F3_{cultivar}_irr_Low_1981_2010.tif']), gdal.GA_ReadOnly)
            rb = ds.GetRasterBand(1)
            self.gaez_f3_low.append(rb.ReadAsArray())

        ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','SoilNutrient','NutrientRetention_High_1981_2010.tif']), gdal.GA_ReadOnly)
        rb = ds.GetRasterBand(1)
        self.gaez_soilconst_high=rb.ReadAsArray()
        
        ds = gdal.Open(ac.GetPath([self.data_folder,'gaez','SoilNutrient','NutrientAvailability_Low_1981_2010.tif']), gdal.GA_ReadOnly)
        rb = ds.GetRasterBand(1)
        self.gaez_soilconst_low=rb.ReadAsArray()

        self.gaez_LAIHI=pd.read_excel(ac.GetPath([self.data_folder,'gaez','GAEZ_LAI_HI.xlsx']),index_col='Crop/LUT')
        
        try:
            self.soil_clay=[]
            for item in ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']:
                ds = gdal.Open(ac.GetPath([self.data_folder,'soil','soilgrids',f'clay_{item}_mean_5000.tif']), gdal.GA_ReadOnly)
                rb = ds.GetRasterBand(1)
                self.soil_clay.append(rb.ReadAsArray())
                
            self.soil_sand=[]
            for item in ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']:
                ds = gdal.Open(ac.GetPath([self.data_folder,'soil','soilgrids',f'sand_{item}_mean_5000.tif']), gdal.GA_ReadOnly)
                rb = ds.GetRasterBand(1)
                self.soil_sand.append(rb.ReadAsArray())
                
            self.soil_soc=[]
            for item in ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']:
                ds = gdal.Open(ac.GetPath([self.data_folder,'soil','soilgrids',f'soc_{item}_mean_5000.tif']), gdal.GA_ReadOnly)
                rb = ds.GetRasterBand(1)
                self.soil_soc.append(rb.ReadAsArray())
                
            self.soil_bd=[]#bulk density
            for item in ['0-5cm','5-15cm','15-30cm','30-60cm','60-100cm','100-200cm']:
                ds = gdal.Open(ac.GetPath([self.data_folder,'soil','soilgrids',f'bdod_{item}_mean_5000.tif']), gdal.GA_ReadOnly)
                rb = ds.GetRasterBand(1)
                self.soil_bd.append(rb.ReadAsArray())
            
            ds = gdal.Open(ac.GetPath([self.data_folder,'soil','soilgrids','id30forsoilgrids.tif']), gdal.GA_ReadOnly)
            rb = ds.GetRasterBand(1)
            self.soilgrids_id30=rb.ReadAsArray()
        except:
            pass
        
        if self.conf.tuned==1:
            self.scaling=pd.read_csv(r'outputs\{0}_scaling_tuned{1}.csv'.format(self.conf.crop_name_4code,0),index_col='faostat')
            self.gridcells30ctry = np.array([ac.ReadNC(ac.GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'id30').flatten(),
                   ac.ReadNC(ac.GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'faostat').flatten()])
    
            self.gridcells30ctry = pd.DataFrame(self.gridcells30ctry.transpose(), columns=['id30','faostat'])
            self.gridcells30ctry = self.gridcells30ctry.where(self.gridcells30ctry['id30']>0).dropna()[['id30', 'faostat']].drop_duplicates().astype(int)
            self.gridcells30ctry=self.gridcells30ctry.set_index('id30')

    def calib_soilfertility(self):
        
        # Get already simulated cells (if some cells didn't run by any reason)
        simulated_ids1 = []; simulated_ids2 = []; simulated_ids3 = []; simulated_ids4 = []; simulated_ids5 = []; simulated_ids6 = []
        if not self.conf.rerun_simulated_cells:
            try:
                simulated_files = os.listdir(ac.GetPath(["outputs", f"{self.conf.project_name}",'Para_SF']))
                for f in simulated_files:
                    _cell_id = int(''.join(filter(str.isdigit, f.split('_')[-2])))
                    _scenario = f.split('_')[-1].split('.')[0]
                    if _scenario == 'rainfed': simulated_ids1.append(_cell_id);simulated_ids2.append(_cell_id)
                    elif _scenario == 'irr': 
                        simulated_ids3.append(_cell_id);simulated_ids4.append(_cell_id)
                        simulated_ids5.append(_cell_id);simulated_ids6.append(_cell_id)
                simulated_ids1 = set(simulated_ids1); simulated_ids2 = set(simulated_ids2)
                simulated_ids3 = set(simulated_ids3); simulated_ids4 = set(simulated_ids4)
                simulated_ids5 = set(simulated_ids5); simulated_ids6 = set(simulated_ids6)
            except:
                self.conf.rerun_simulated_cells = True
        
        # Get gridcells
        gridcells = ac.GetSimulationCells30(self.conf.gridcells) # List of gridcells to consider
        for cell_id in gridcells*1: # Check which cells are valid
            cell = self.gridcells30[self.gridcells30['id30']==cell_id]
            if not self.CheckCell(cell, simulated_ids1, simulated_ids2, simulated_ids3, simulated_ids4, simulated_ids5, simulated_ids6): gridcells.remove(cell_id)

        # Run specific crop model
        print(f'Running your calibrations... ({len(gridcells)} cells)')
        if self.conf.crop_model == 'AquaCrop': # AOS6
            RunACEA_partial = partial(InitACEA, self)
            
            # Iterate gridcells
            if self.conf.multi_core:
                print('< Running in parallel >')
                import random
                random.shuffle(gridcells)
                from multiprocessing import get_context
                cells_count = []
                with get_context("spawn").Pool(self.conf.CPUs) as p:
                    r = p.map_async(RunACEA_partial, gridcells, callback=cells_count.extend)
                    r.wait()
                
            else:
                cells_count = list(filter(None, map(RunACEA_partial, gridcells)))

        else:
            print(f'Crop model {self.conf.crop_model} cannot be called')
        
        print('---------------------------')
        print(f'Done in {(time.time() - self.tic):.1f}s! Number of run cells: {len(set(cells_count))}')

    
    def Run(self):
        
        # Get already simulated cells (if some cells didn't run by any reason)
        simulated_ids1 = []; simulated_ids2 = []; simulated_ids3 = []; simulated_ids4 = []; simulated_ids5 = []; simulated_ids6 = []
        if not self.conf.rerun_simulated_cells:
            try:
                simulated_files = os.listdir(ac.GetPath(["outputs", f"{self.conf.project_name}",'annual']))
                for f in simulated_files:
                    _cell_id = int(''.join(filter(str.isdigit, f.split('_')[-1])))
                    _scenario = int(''.join(filter(str.isdigit, f.split('_')[-2])))
                    if _scenario == 1: simulated_ids1.append(_cell_id)
                    elif _scenario == 2: simulated_ids2.append(_cell_id)
                    elif _scenario == 3: simulated_ids3.append(_cell_id)
                    elif _scenario == 4: simulated_ids4.append(_cell_id)
                    elif _scenario == 5: simulated_ids5.append(_cell_id)
                    elif _scenario == 6: simulated_ids6.append(_cell_id)
                simulated_ids1 = set(simulated_ids1); simulated_ids2 = set(simulated_ids2)
                simulated_ids3 = set(simulated_ids3); simulated_ids4 = set(simulated_ids4)
                simulated_ids5 = set(simulated_ids5); simulated_ids6 = set(simulated_ids6)
            except:
                self.conf.rerun_simulated_cells = True
        
        # Get gridcells
        gridcells = ac.GetSimulationCells30(self.conf.gridcells) # List of gridcells to consider
        for cell_id in gridcells*1: # Check which cells are valid
            cell = self.gridcells30[self.gridcells30['id30']==cell_id]
            if not self.CheckCell(cell, simulated_ids1, simulated_ids2, simulated_ids3, simulated_ids4, simulated_ids5, simulated_ids6): gridcells.remove(cell_id)
        
        # Run specific crop model
        print(f'Running your simulations... ({len(gridcells)} cells)')
        if self.conf.crop_model == 'AquaCrop': # AOS6
            RunACEA_partial = partial(RunACEA, self)
            
            # Iterate gridcells
            if self.conf.multi_core:
                print('< Running in parallel >')
                import random
                random.shuffle(gridcells)
                from multiprocessing import get_context
                cells_count = []
                with get_context("spawn").Pool(self.conf.CPUs) as p:
                    r = p.map_async(RunACEA_partial, gridcells, callback=cells_count.extend)
                    r.wait()
                
                # from multiprocessing import Pool
                # async_results = []
                # pool = Pool(processes=self.conf.CPUs)

                # # Run AquaCrop-OS
                # for cell_id in gridcells: # Iterate raster cells
                #     r = pool.apply_async(RunACEA_v6, (self, int(cell_id)), callback=_pool_callback)
                #     async_results.append(r)
                    
                # for result in async_results:
                #     r.wait()
                    
                # pool.close()
                
            else:
                cells_count = list(filter(None, map(RunACEA_partial, gridcells)))

        else:
            print(f'Crop model {self.conf.crop_model} cannot be called')
        
        print('---------------------------')
        print(f'Done in {(time.time() - self.tic):.1f}s! Number of run cells: {len(set(cells_count))}')
        
    def Rasterise(self):
        "Save main results to gridded raster"
        
        print('Saving results as global rasters...')
        
        # Preparation
        folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",f"{self.conf.virtual_irrigation}",'annual'])
        gs_to_consider = (self.clock_end_year - self.clock_start_year)+1 - self.conf.spinup
        _,_,res = ac.GetResolution() # Get 30 arcmin resolution
        
        variables = ac.GetListOfVaribales() # Get variables
        scenarios = ac.GetListOfScenarios() # Get scenarios
        
        for sc, sc_desc in scenarios.items():
            cells_to_save = []
            cwu_green, cwu_blue_cr, cwu_blue_ir, yields, gdd, aet, biom, pirnreq, soilevap, transp, runoff, plantday, plantyear, anthday, matyday, harvyear, = \
                (np.ones((gs_to_consider, *res)) * np.ma.masked for _ in range(16)) # Create empty arrays with results
            
            # Read each cell
            for _,_, files in os.walk(folder_path):
                
                for f in files:
                    filename = f.split('_')
                    #print(filename)
                    if f'sc{sc}' in f:
                        cell_id = int(filename[-1].split('.')[0])
                        cell = self.gridcells30[self.gridcells30['id30']==cell_id]
                        if len(cell)==0: print(f'Cell {cell_id} is not in a grid'); continue
                        row = cell['rowy30']; col = cell['colx30']
                        data = np.load(ac.GetPath([folder_path,f]), allow_pickle=True) # Read the output file

                        # Crop growth variables to save
                        insert_limit = len(data['general'][-gs_to_consider:,5].reshape(-1,1))
                        yields[:insert_limit, row, col] = data['general'][-gs_to_consider:,5].reshape(-1,1) # Crop yields [dry ton ha-1 gs-1]
                        biom[:insert_limit, row, col] = data['general'][-gs_to_consider:,6].reshape(-1,1) # Crop biomass [dry ton ha-1 gs-1]
                        gdd[:insert_limit, row, col] = data['general'][-gs_to_consider:,7].reshape(-1,1) # Accumulated GDDs [GDDs]
                        
                        _p_days = data['general'][-gs_to_consider:,2].reshape(-1,1)[:,0] # Planting day
                        _a_days = data['general'][-gs_to_consider:,3].reshape(-1,1)[:,0] # Anthesis day
                        _m_days = data['general'][-gs_to_consider:,4].reshape(-1,1)[:,0] # Harvest day
                        
                        _p_d = [d.timetuple().tm_yday for d in _p_days]
                        ant_year = [d.year for d in _a_days]
                        plantday[:insert_limit, row, col] = [[d.timetuple().tm_yday] for d in _p_days]
                        plantyear[:insert_limit, row, col] = [[d.year] for d in _p_days]; harvyear[:insert_limit, row, col] = [[d.year] for d in _m_days]
                        
                        _p_d_ori = _p_d*1
                        for y in range(insert_limit):
                            if ant_year[y] != int(plantyear[y,row, col]):
                                _p_d[y] = _p_d[y] - ((datetime.datetime(_p_days[y].year,12,31) - datetime.datetime(_p_days[y].year,1,1)).days+1)
                            
                        anthday[:insert_limit, row, col] = [[a] for a in np.array([d.timetuple().tm_yday for d in _a_days]) - np.array(_p_d)] # Anthesis
                        
                        for y in range(insert_limit):
                            if int(harvyear[y,row, col]) != int(plantyear[y,row, col]):
                                _p_d_ori[y] = _p_d_ori[y] - ((datetime.datetime(_p_days[y].year,12,31) - datetime.datetime(_p_days[y].year,1,1)).days+1)
                                
                        matyday[:insert_limit, row, col] = [[a] for a in np.array([d.timetuple().tm_yday for d in _m_days]) - np.array(_p_d_ori)] # Harvest
                        
                        if int(np.min(anthday[:, row, col])) < 0: raise print("Error with anthesis")
                        elif np.max(anthday[:, row, col]) > 360: print(f'Weird anthesis date {np.max(anthday[:, row, col])} for {cell_id}')
                            
                        if int(np.min(matyday[:, row, col])) < 0 : raise Exception("Error with maturity")
                            
                        # Water balance variables to save
                        pirnreq[:insert_limit, row, col] = data['general'][-gs_to_consider:,10].reshape(-1,1) # Irrigation demand [mm] or [kg m-2 gs-1]
                        soilevap[:insert_limit, row, col] = data['general'][-gs_to_consider:,14].reshape(-1,1) # Evaporation [mm] or [kg m-2 gs-1]
                        transp[:insert_limit, row, col] = data['general'][-gs_to_consider:,15].reshape(-1,1) # Transpiration [mm] or [kg m-2 gs-1]
                        aet[:insert_limit, row, col] = soilevap[:insert_limit, row, col] + transp[:insert_limit, row, col] # Evapotranspiration [mm] or [kg m-2 gs-1]
                        runoff[:insert_limit, row, col] = data['general'][-gs_to_consider:,16].reshape(-1,1) # Runoff [mm] or [kg m-2 gs-1]
                        
                        # CWU [mm] variables to save
                        cwu_green[:insert_limit, row, col] = data['cwu'][-gs_to_consider:,0].reshape(-1,1)
                        cwu_blue_ir[:insert_limit, row, col] = data['cwu'][-gs_to_consider:,1].reshape(-1,1)
                        cwu_blue_cr[:insert_limit, row, col] = data['cwu'][-gs_to_consider:,2].reshape(-1,1)
                        
                        cells_to_save.append(cell_id)
                        #if len(cells_to_save)%1000==0: print(len(cells_to_save))
                        
            # biom[biom == 0] = .001
            # hi[biom!=1.e20] = np.around(yields[biom!=1.e20]/biom[biom!=1.e20],3) # Harvest index•
            
            # Save rasters
            if len(cells_to_save)>0:
                
                raster_folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",f"{self.conf.virtual_irrigation}",'rasters'])
                if not os.path.exists(raster_folder_path): os.makedirs(raster_folder_path)
                
                first_year = self.clock_start_year + self.conf.spinup
                for var, var_info in variables.items(): # Do it dynamically
                    if var == 'yield': var = 'yields' # Python claims yield as a name already
                    if (sc < 3) and (var in ['cwu_blue_ir', 'pirnreq']): continue
                    if (sc not in [2,6]) and (var == 'cwu_blue_cr'): continue
                
                    if var in vars():
                        data = vars()[var]
                        if var == 'yields': var = 'yield' 
                        
                        val_name = f'{var}-{self.conf.crop_name_short}'; val_unit = var_info[1]
                        raster_path = ac.GetPath([raster_folder_path, f'acea_30arc_sc{sc}_{self.conf.crop_name_short}_{var}_global_annual_{first_year}_{self.clock_end_year}.nc'])
                        title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Scenario: {sc_desc}, Variable: {var_info[0]}'
                        
                        ac.CreateHistRaster30(data, raster_path, val_name, val_unit, title, first_year, gs_to_consider)
        
#%% General supportive functions
    def CheckCell(self, cell,  simulated_ids1, simulated_ids2, simulated_ids3, simulated_ids4, simulated_ids5, simulated_ids6):
        "Check is this cell is valid"
        
        if len(cell)!=1:
            return False
        
        if self.conf.real_cropland: # Check if the crop grows there
            scenarios = self.AdjustScenarios(cell)
            if len(scenarios)==0: return False
        else:scenarios = self.conf.scenarios
            
        # Check if weather file exsits
        if not os.path.exists(ac.GetPath(["data", "acea",'climate',f'{self.conf.climate_name}_{int(cell["id30"])}.pckl'])):
            return False
        
        # Check if any crop calendar exsits
        if not any(list(self.GetCropCycleDates(cell, False)+self.GetCropCycleDates(cell, True))):
            return False
        
        if not self.conf.rerun_simulated_cells: # Check if the cell was already simulated
            _res = 0
            for _sc in scenarios:
                if _sc == 1 and int(cell["id30"]) in simulated_ids1: _res += 1
                if _sc == 2:
                    if int(cell["id30"]) in simulated_ids2: _res += 1
                    elif len(self.GetGroundwaterLevels(cell)[1]) == 0: _res += 1
                if _sc == 3 and int(cell["id30"]) in simulated_ids3: _res += 1
                if _sc == 4 and int(cell["id30"]) in simulated_ids4: _res += 1
                if _sc == 5 and int(cell["id30"]) in simulated_ids5: _res += 1
                if _sc == 6 and int(cell["id30"]) in simulated_ids6: _res += 1
            if _res == len(scenarios): return False
            
        return True
    
    def AdjustScenarios(self, cell):
        "Check is this cell is valid"
        
        scenarios = np.array(self.conf.scenarios)
        
        if self.conf.real_cropland: # if real cropland, select only real scenarios
            [x5, y5] = [int(cell['colx30'])*6, int(cell['rowy30'])*6]
            x_range5 = [x5, x5+6]; y_range5 = [y5, y5+6]
        
            if min(scenarios) < 3:
                if not((self.harvested_areas_rain5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]].mask==False).any()):
                    scenarios = scenarios[scenarios>2] # No rainfed production in this cell
                    if len(scenarios)==0: return list(scenarios)
            if max(scenarios) > 2:
                if not((self.harvested_areas_ir5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]].mask==False).any()):
                    scenarios = scenarios[scenarios<3] # No irrigated production in this cell
                
        return list(scenarios)
    
    
    def class_soil_texturebased(self,sand_,clay_):
        
        c=clay_
        sa=sand_
        si=100-c-sa
        
        if si + (1.5 * c) < 15:
            return 'Sand'
        elif si + 1.5 * c >= 15 and si + 2 * c < 30:
            return 'Loamy Sand'
        elif ((c >= 7 and c < 20) or (sa > 52) and ((si + 2*c) >= 30) or (c < 7 and si < 50 and (si+2*c) >= 30)):
            return 'Sandy Loam'
        elif ((c >= 7 and c < 27) and (si >= 28 and si < 50) and (sa <= 52)):
            return 'Loam'
        elif ((si >= 50 and (c >= 12 and c < 27)) or ((si >= 50 and si < 80) and c < 12)):
            return 'Silt Loam'
        elif (si >= 80 and c < 12):
            return 'Silt'
        elif ((c >= 20 and c < 35) and (si < 28) and (sa > 45)):
            return 'Sandy Clay Loam'
        elif ((c >= 27 and c < 40) and (sa > 20 and sa <= 45)):
            return 'Clay Loam'
        elif ((c >= 27 and c < 40) and (sa  <= 20)):
            return 'Silty Clay Loam'
        elif (c >= 35 and sa > 45):
            return 'Sandy Clay'
        elif (c >= 40 and si >= 40):
            return 'Silty Clay'
        elif (c >= 40 and sa <= 45 and si < 40):
            return 'Clay'

    def penetrability_estimator(self,soil_class_,bd_):
        if soil_class_ in ['Sand','Loamy Sand']:
            l1=1.60
            l2=1.69
            l3=1.80
            
        elif soil_class_ in ['Sandy Loam','Loam']:
            l1=1.40
            l2=1.63
            l3=1.80
            
        elif soil_class_ in ['Sandy Clay Loam','Clay Loam']:
            l1=1.40
            l2=1.60
            l3=1.75
            
        elif soil_class_ in ['Silt','Silt Loam']:
            l1=1.30
            l2=1.60
            l3=1.75
            
        elif soil_class_ in ['Silty Clay Loam']:
            l1=1.40
            l2=1.55
            l3=1.65
            
        elif soil_class_ in ['Sandy Clay','Silty Clay']:
            l1=1.10
            l2=1.49
            l3=1.58
            
        elif soil_class_ in ['Clay']:
            l1=1.10
            l2=1.39
            l3=1.47
            
        else:
            return np.nan
        
        if bd_<=l1:
            return 100
        elif bd_<=l2:
            return 75
        elif bd_<=l3:
            return 50
        elif bd_>l3:
            return 35
        else:
            return np.nan

        
    def GetSoilTexture(self, cell,use_soilgrids=0):
        "Get sand and clay fractions in a specific cell"
        if use_soilgrids==0:
            sand = self.sand_data[int(cell['rowy30']),int(cell['colx30'])]
            clay = self.clay_data[int(cell['rowy30']),int(cell['colx30'])]
            # silt = 100 - sand - clay
            
            # Limit sand and clay to values acceptable for 
            while sand > 93: sand -= 1.; clay += 1.
            while clay > 93: clay -= 1.; sand += 1.
            
            while clay < 4: clay += 1.
                
            if (sand+clay) > 100: raise ValueError('Something wrong with the soil') # Check if together less than 100%
            
            return sand, clay
        else:
            loc_=self.soilgrids_id30==cell['id30'].values[0]
            
            clay=[]
            sand=[]
            om=[]
            bd=[]
            pene=[]
            #print(cell['id30'].values[0])
            
            for i in range(6):
            
                temp1=self.soil_clay[i][loc_]
                
                
                temp2=self.soil_sand[i][loc_]
                
                
                if np.sum(temp1>0) and np.sum(temp2>0):
                    clay.append(temp1.mean(where=temp1>0)/10)
                    sand.append(temp2.mean(where=temp2>0)/10)
                else:

                    clay.append(self.clay_data[int(cell['rowy30']),int(cell['colx30'])])
                    sand.append(self.sand_data[int(cell['rowy30']),int(cell['colx30'])])
                    while sand[-1] > 93: sand[-1] -= 1.; clay[-1] += 1.
                    while clay[-1] > 93: clay[-1] -= 1.; sand[-1] += 1.
                    
                    while clay[-1] < 4: clay[-1] += 1.
                    if (sand[-1]+clay[-1]) > 100: raise ValueError('Something wrong with the soil')


                temp=self.soil_soc[i][loc_]
                if np.sum(temp>0):
                    om.append(temp.mean(where=temp>0)/10000*2*100)#convert to organic matter %w
                else:
                    om.append(2.5)
                
                temp=self.soil_bd[i][loc_]
                if np.sum(temp>0):
                    bd.append(temp.mean(where=temp>0)/100)#bulk density in Mg.m-3, need to be translated into penetrability table2.41b
                else:
                    bd.append(-1)
                    
                soil_class=self.class_soil_texturebased(sand[-1],clay[-1])
                pene.append(self.penetrability_estimator(soil_class,bd[-1]))
                
                if np.isnan(pene[-1]):
                    pene[-1]=100
            #assume 100% penetrability for all soil layers due to large uncertanties in this parameter
            for i in range(6):
                if pene[i]<100:
                    pene[i]=100

            #for i in [4,3,2,1,0]:
            #    if pene[i]<pene[i+1]:
            #        pene[i]=pene[i+1]
            #pene=[100,100,100,100,100,100]
            #print(sand, clay, om, pene,bd)
            return sand, clay, om, pene
            
            
        
    def GetGroundwaterLevels(self, cell, rainfed=True, option=1):
        "Calculate appx GW levels in a specific 30 arcmin cell based on 5 arc data"
        
        dates = list(); values = list()
        [x5, y5] = [int(cell['colx30'])*6, int(cell['rowy30'])*6]
        x_range5 = [x5, x5+6]; y_range5 = [y5, y5+6]
        if rainfed:crop_cells = self.harvested_areas_rain5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
        else:crop_cells = self.harvested_areas_ir5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
        crop_locations = crop_cells*1; crop_locations[crop_locations>0] = 1
        shallow_gw_monthly = self.gw_data[:, y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]  # Groundwater levels at 5 arc mins
        shallow_gw_monthly = shallow_gw_monthly*crop_locations # Mask cells without a crop
        
        # Mask deep cells
        for row in range(6):
            for col in range(6):
                _gw_cell5 = shallow_gw_monthly[:,row,col]*1
                if all(_gw_cell5.mask == True): continue
                if np.ma.max(_gw_cell5) <= self.conf.gw_min_level:
                    shallow_gw_monthly[:,row,col] = np.ma.masked
        
        # Limit values close to surface (farmers drained the field) t 5 arc min
        shallow_gw_monthly[(shallow_gw_monthly>self.conf.gw_max_level) & (shallow_gw_monthly.mask == False)] = self.conf.gw_max_level
        
        if option == 1: # Just an average of cells with a crop
            shallow_gw_monthly = np.around(np.ma.average(shallow_gw_monthly.reshape((12,-1)), axis = 1),3)
        elif option == 2: # Take a weighted average of cells with a crop
            shallow_gw_monthly = np.around(np.ma.average(shallow_gw_monthly.reshape((12,-1)), axis = 1, weights = crop_cells.flatten(),),3) # Weighted avg of all 5 arcmin cells
        
        if not all(shallow_gw_monthly.mask):
            for y in range(self.clock_start_year, self.clock_end_year+1):
                for m in range(1,12+1):
                    dates.append(f'{y}/{m}/15')
                    values.append(-shallow_gw_monthly[m-1])
           
            # Add start and end dates
            dates = [f'{self.clock_start_year}/1/1']+dates;values = [-shallow_gw_monthly[0]]+values
            dates.append(f'{y}/{m}/31');values.append(-shallow_gw_monthly[0])
        return dates, values
    
    def GetBiomassCanopySF_once(self, cell, rainfed=True,search_=0):
        "Calculate biomass and canopy cover under soil fertility stress in a specific 30 arcmin cell based on 5 arc data"
        
        cur_search=search_
        [x5, y5] = [int(cell['colx30'])*6, int(cell['rowy30'])*6]
        x_range5 = [x5-search_, x5+6+search_]; y_range5 = [y5-search_, y5+6+search_]
        if rainfed:crop_cells = self.harvested_areas_rain5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
        else:crop_cells = self.harvested_areas_ir5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
        
        crop_locations = crop_cells*1; crop_locations[crop_locations>0] = 1
        
        #average over space and cultivars
        RelBiomass=[]
        RelCanopy=[]
        f3low=[]
        f3high=[]
        
        for cultivars in range(len(self.conf.crop_gaez)):
        
            lut_high=self.gaez_lut_high[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            lut_low=self.gaez_lut_low[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            yield_high=self.gaez_yield_high[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            yield_low=self.gaez_yield_low[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            f3_high=self.gaez_f3_high[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            f3_low=self.gaez_f3_low[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            
            soilconst_high=self.gaez_soilconst_high[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            soilconst_low=self.gaez_soilconst_low[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations



            for row in range(lut_high.shape[0]):
                for col in range(lut_high.shape[1]):

                    
                    if lut_high[row,col]>0 and lut_low[row,col]>0 and yield_high[row,col]>0 and yield_low[row,col]>0 and\
                        f3_high[row,col]>0 and f3_low[row,col]>0 and soilconst_high[row,col]>0 and\
                            soilconst_high[row,col]<=10 and soilconst_low[row,col]>0 and soilconst_low[row,col]<=10:

                                HI_high=self.gaez_LAIHI.loc[int(lut_high[row,col]),'Harvest index (high inputs)']
                                HI_low=self.gaez_LAIHI.loc[int(lut_low[row,col]),'Harvest index (low inputs)']
                                LAI_high=self.gaez_LAIHI.loc[int(lut_high[row,col]),'Maximum leaf area index (high inputs)']
                                LAI_low=self.gaez_LAIHI.loc[int(lut_low[row,col]),'Maximum leaf area index (low inputs)']
                                
                                #Biomass_high=yield_high[row,col]/f3_high[row,col]*10000*soilconst_high[row,col]/10/HI_high
                                Biomass_high=yield_high[row,col]/f3_high[row,col]*10000/HI_high
                                Biomass_low=yield_low[row,col]/f3_low[row,col]*10000*soilconst_low[row,col]/10/HI_low
                                
                                #Fix extreme cases for rubber
                                if Biomass_low>Biomass_high:
                                    Biomass_low=Biomass_high

                                #a*(1-np.exp(b*LAI))**c
                                CC_high=self.conf.LAI_CC_a*(1-np.exp(self.conf.LAI_CC_b*LAI_high))**self.conf.LAI_CC_c
                                CC_low=self.conf.LAI_CC_a*(1-np.exp(self.conf.LAI_CC_b*LAI_low))**self.conf.LAI_CC_c
                                
                                RelBiomass.append(Biomass_low/Biomass_high)
                                RelCanopy.append(CC_low/CC_high)
                                f3low.append(f3_low[row,col]/10000)
                                f3high.append(f3_high[row,col]/10000)
        #if GAEZ does not have the data, then find the neighbours; if neighbours also don't have the data, return empty

        if RelBiomass==[]:
            return -1,-1,-1,-1,-1
        else:
            return np.mean(RelBiomass), np.mean(RelCanopy),np.mean(f3low), np.mean(f3high),cur_search
        
    def GetBiomassCanopySF(self, cell, rainfed=True):
        cur_search=0
        bio,cnp,f3low,f3high,search=self.GetBiomassCanopySF_once(cell, rainfed=rainfed,search_=cur_search)
        while search==-1 and cur_search+6<=120:
            cur_search+=6
            bio,cnp,f3low,f3high,search=self.GetBiomassCanopySF_once(cell, rainfed=rainfed,search_=cur_search)
        return bio,cnp,f3low,f3high,search
    
    def GetfromGAEZ_once(self,cell, irrigated,search_,varname):
        "get GAEZ data in a specific 30 arcmin cell based on 5 arc data"
        
        cur_search=search_
        [x5, y5] = [int(cell['colx30'])*6, int(cell['rowy30'])*6]
        x_range5 = [x5-search_, x5+6+search_]; y_range5 = [y5-search_, y5+6+search_]
        if irrigated==0:crop_cells = self.harvested_areas_rain5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
        else:crop_cells = self.harvested_areas_ir5[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]
        
        crop_locations = crop_cells*1; crop_locations[crop_locations>0] = 1
        
        #average over space and cultivars
        values=[]
        
        for cultivars in range(len(self.conf.crop_gaez)):
        
            lut_high=self.gaez_lut_high[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            lut_low=self.gaez_lut_low[cultivars][y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations

            soilconst_high=self.gaez_soilconst_high[y_range5[0]:y_range5[1], x_range5[0]:x_range5[1]]*crop_locations
            
            for row in range(6+2*search_):
                for col in range(6+2*search_):

                    if varname=='relbio_high':
                        if soilconst_high[row,col]>0 and soilconst_high[row,col]<=10:
                            values.append(soilconst_high[row,col]/10)
                            
                    elif varname=='hi_high':
                        if lut_high[row,col]>0:
                            HI=self.gaez_LAIHI.loc[int(lut_high[row,col]),'Harvest index (high inputs)']
                            values.append(HI)
                        
                    elif varname=='hi_low':
                        if lut_low[row,col]>0:
                            HI=self.gaez_LAIHI.loc[int(lut_low[row,col]),'Harvest index (low inputs)']
                            values.append(HI)
                        
                    elif varname=='cc_high':
                        if lut_high[row,col]>0:
                            LAI=self.gaez_LAIHI.loc[int(lut_high[row,col]),'Maximum leaf area index (high inputs)']
                            CC=self.conf.LAI_CC_a*(1-np.exp(self.conf.LAI_CC_b*LAI))**self.conf.LAI_CC_c
                            values.append(CC)

                                
        #if GAEZ does not have the data, then find the neighbours; if neighbours also don't have the data, return empty

        if values==[]:
            return -1,-1
        else:
            return np.mean(values), cur_search
    
    
    def GetfromGAEZ(self,cell, irrigated,varname):
        cur_search=0
        value,search=self.GetfromGAEZ_once(cell, irrigated,search_=cur_search,varname=varname)
        while search==-1 and cur_search+6<=120:
            cur_search+=6
            value,search=self.GetfromGAEZ_once(cell, irrigated,search_=cur_search,varname=varname)
        return value,search
        

    def GetCropCycleDates(self, cell, irrigated):
        "Get cell-specific data in a specific cell"
        if irrigated:
            try:
                planting_day = int(self.crop_planting_ir[cell['rowy30'],cell['colx30']][0])
                harvest_day = int(self.crop_harvest_ir[cell['rowy30'],cell['colx30']][0])
            except:
                planting_day = False; harvest_day = False
                
        else:
            try:
                planting_day = int(self.crop_planting_rf[cell['rowy30'],cell['colx30']][0])
                harvest_day = int(self.crop_harvest_rf[cell['rowy30'],cell['colx30']][0])
            except:
                planting_day = False; harvest_day = False
        
        if (planting_day != False) and (harvest_day!= False):
            planting_day = datetime.datetime(2001, 1, 1) + datetime.timedelta(planting_day - 1)
            planting_day = f'{planting_day.month}/{planting_day.day}'
            
            harvest_day = datetime.datetime(2001, 1, 1) + datetime.timedelta(harvest_day - 1)
            harvest_day = f'{harvest_day.month}/{harvest_day.day}'
        
        return planting_day, harvest_day
    
    def DisplaySummary(self, outputs, sc):
        "Show short summary of a cell run"
        avg_yield = np.around(np.mean(outputs.final_stats['Dry yield (tonne/ha)'].values),1)
        avg_cwu_g = np.around(np.mean(outputs.final_watercolor['ET green (mm)'].values),1)
        avg_cwu_b = np.around(np.mean(outputs.final_watercolor['ET blue irr (mm)'].values + outputs.final_watercolor['ET blue cr (mm)'].values),1)
        avg_ir = np.around(np.mean(outputs.final_stats['Irrigation (mm)'].values),1)
        avg_cr = np.around(np.mean(outputs.final_stats['Capillary rise (mm)'].values),1)
        avg_e = np.around(np.mean(outputs.final_stats['Evaporation (mm)'].values),1)
        avg_t = np.around(np.mean(outputs.final_stats['Transpiration (mm)'].values),1)
        avg_aet = np.around(np.mean(outputs.final_stats['Evaporation (mm)'].values + outputs.final_stats['Transpiration (mm)'].values),1)
        print(f"\tScenario {sc}: yield = {avg_yield} [t ha-1], CWUg = {avg_cwu_g}, CWUb = {avg_cwu_b}, irrigation = {avg_ir}, CR = {avg_cr}, E = {avg_e}, T = {avg_t}, AET = {avg_aet}")
        
    def SaveRawResults(self, outputs, sc, cell):
        "Save summary of a cell run"
        
        cell_id = int(cell["id30"])
        # Create results folder
        if self.conf.project_save_annual or self.conf.project_save_daily:
            folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",f"{self.conf.virtual_irrigation}"])
            if not os.path.exists(folder_path): os.makedirs(folder_path)
        
        # Save the summary of each growing season
        if self.conf.project_save_annual:
            folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",f"{self.conf.virtual_irrigation}",'annual'])
            if not os.path.exists(folder_path): os.makedirs(folder_path)
            
            watercolor = np.around(outputs.final_watercolor,3) # Crop water use [mm]
            #sws = np.around(outputs.final_watercolor.values[:,4:],3) # Soil water storage [mm]
            general = np.around(outputs.final_stats,3)  # Crop growth and water balance

            file_name = ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_annual_sc{sc}_{cell_id}.npz'])
            with open(file_name, 'wb') as f:
                #np.savez_compressed(f, cwu=cwu,sws=sws,general=general) 
                np.savez_compressed(f, watercolor=watercolor,general=general)
        
        # Save the daily outputs of the whole simulation period
        if self.conf.project_save_daily:
            folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",f"{self.conf.virtual_irrigation}",'daily'])
            if not os.path.exists(folder_path): os.makedirs(folder_path)
            
            flux = np.around(outputs.water_flux,3) # Soil water fluxes [mm]
            crop = np.around(outputs.crop_growth,3) # Crop growth
            wc = np.around(outputs.water_storage,3)  # Soil water storage [mm]
            
            et_color = np.around(outputs.ET_color,3)  # Soil water storage [mm]
            s_color = np.around(outputs.S_color,3)  # Soil water storage [mm]
            
            file_name = ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_daily_sc{sc}_{cell_id}.npz'])
            
            with open(file_name, 'wb') as f:
                np.savez_compressed(f, flux = flux, crop = crop, wc = wc,
                                    et_color =et_color, s_color = s_color)
                
                
    def SaveParametersSF(self, outputs,f3low,f3high,cur_search, cell,rainfed):
        "Save calibrated parameters of a cell, per crop"
        
        cell_id = int(cell["id30"])
        # Create results folder

        folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'Para_SF'])
        if not os.path.exists(folder_path): os.makedirs(folder_path)
        
        # Save the summary of each growing season
        if outputs==[]:
            cur_search=cur_search
            file_name = ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{cell_id}_{rainfed}.npz'])
            with open(file_name, 'wb') as f:
                np.savez_compressed(f, cur_search=cur_search) 

        else:             
            Ksccx=outputs.crop.Ksccx_in
            Ksexpf=outputs.crop.Ksexpf_es[0]
            fcdecline=outputs.crop.fcdecline_es[0]
            Kswp=outputs.crop.Kswp_es[0]
            sfertstress=1-outputs.crop.RelativeBio
            cur_search=cur_search
            
            sf_es=outputs.crop.sf_es
            Ksexpf_es=outputs.crop.Ksexpf_es
            fcdecline_es=outputs.crop.fcdecline_es
            Kswp_es=outputs.crop.Kswp_es
            Ksccx_es=outputs.crop.Ksccx_es
            relbio_es=outputs.crop.relbio_es
            f3low=f3low
            f3high=f3high
    
            file_name = ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{cell_id}_{rainfed}.npz'])
            with open(file_name, 'wb') as f:
                np.savez_compressed(f, Ksccx=Ksccx,Ksexpf=Ksexpf,fcdecline=fcdecline,Kswp=Kswp,sfertstress=sfertstress,\
                                    cur_search=cur_search,sf_es=sf_es,Ksexpf_es=Ksexpf_es,fcdecline_es=fcdecline_es,\
                                        Kswp_es=Kswp_es,Ksccx_es=Ksccx_es,relbio_es=relbio_es,f3low=f3low,f3high=f3high) 
         
    def RasterParametersSF(self):
        "Rasterize saved calibrated parameters per crop"
        
        for item in ['irr','rainfed']:
            
            folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'Para_SF'])
            if not os.path.exists(folder_path): os.makedirs(folder_path)
        
            Ksccx=np.zeros((360,720))-1
            Ksexpf=np.zeros((360,720))-1
            fcdecline=np.zeros((360,720))-1
            Kswp=np.zeros((360,720))-1
            sfertstress=np.zeros((360,720))-1
            cur_search=np.zeros((360,720))-1
            f3low=np.zeros((360,720))-1
            f3high=np.zeros((360,720))-1
            
            file_name = os.listdir(folder_path)
            for file in file_name:
                
                #print(file)
                #print(file.split('_'))
                if file.split('_')[3][0:-4]==item:
                    cell_id=int(file.split('_')[2])
                    cell = self.gridcells30[self.gridcells30['id30']==cell_id]
                    
                    col=cell['colx30']
                    row=cell['rowy30']    
            
                    file_name = ac.GetPath([folder_path, file])
        
                    data = np.load(file_name)
                    cur_search[row,col]=data['cur_search']
                    if cur_search[row,col]>-1:
                        Ksccx[row,col]=data['Ksccx']
                        Ksexpf[row,col]=data['Ksexpf']
                        Kswp[row,col]=data['Kswp']
                        fcdecline[row,col]=data['fcdecline']
                        sfertstress[row,col]=data['sfertstress']
                        f3low[row,col]=data['f3low']
                        f3high[row,col]=data['f3high']
            
            Ksccx[Ksccx==-1]=np.NaN  
            Ksccx=np.ma.masked_invalid(Ksccx)
            Ksexpf[Ksexpf==-1]=np.NaN  
            Ksexpf=np.ma.masked_invalid(Ksexpf)
            Kswp[Kswp==-1]=np.NaN  
            Kswp=np.ma.masked_invalid(Kswp)
            fcdecline[fcdecline==-1]=np.NaN  
            fcdecline=np.ma.masked_invalid(fcdecline)
            sfertstress[sfertstress==-1]=np.NaN  
            sfertstress=np.ma.masked_invalid(sfertstress)
            f3low[f3low==-1]=np.NaN  
            f3low=np.ma.masked_invalid(f3low)
            f3high[f3high==-1]=np.NaN  
            f3high=np.ma.masked_invalid(f3high)
            cur_search[cur_search==-1]=np.NaN  
            cur_search=np.ma.masked_invalid(cur_search)
            
            folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'rastersSF'])
            if not os.path.exists(folder_path): os.makedirs(folder_path)
            print('Saving soil fertility parameters')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_Ksccx.nc'])
            ac.CreateRaster30(Ksccx, path, 'Ksccx', '1', 'Ksccx')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_Ksexpf.nc'])
            ac.CreateRaster30(Ksexpf, path, 'Ksexpf', '1', 'Ksexpf')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_Kswp.nc'])
            ac.CreateRaster30(Kswp, path, 'Kswp', '1', 'Kswp')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_fcdecline.nc'])
            ac.CreateRaster30(fcdecline, path, 'fcdecline', '1', 'fcdecline')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_sfertstress.nc'])
            ac.CreateRaster30(sfertstress, path, 'sfertstress', '1', 'sfertstress')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_cur_search.nc'])
            ac.CreateRaster30(cur_search, path, 'cur_search', '1', 'cur_search')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_f3low.nc'])
            ac.CreateRaster30(f3low, path, 'f3low', '1', 'f3low')
            path=ac.GetPath([folder_path, f'{self.conf.project_name}_{self.conf.crop_fao}_{item}_f3high.nc'])
            ac.CreateRaster30(f3high, path, 'f3high', '1', 'f3high')
            
            
            
            
            
            
