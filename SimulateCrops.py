# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 08:02:29 2022

@author: SuH
"""
from modules.acea.acea_simulation import acea_model
import importlib

#'Apple','ArabicaCoffee',"Avocado",'Banana','Barley','Bean','BerriesNes','Cashew','Cassava','Cereals'\
#'Chickpea','Cocoa','Coconut','Cotton','Cowpea','FibreCropsNes','ForageSilageAlfalfa',\
#'FruitTropical','Grapes','Groundnut','Lentils','Maize','Oat','Oilpalm','Oilseeds','Olive',\
#'Onion','Oranges','Pea','PearlMillet','Pepper','PigeonPea','Plantains','Potato','Pulses',\
#'Rapeseed','Rice1','Rice2','RobustaCoffee','Rubber','Sesame','Sheanuts','SmallMillet','Sorghum',\
#'Soybean','SpringWheat','Sugarbeet','Sugarcane','Sunflower','SweetPotato','Taro','Tea',\
#'Tobacco','Tomato','Vegetables','VegetablesLeafy','VegetablesTuber','WinterWheat','Yams'


croplist=['Rubber']
    
    
for crop_temp in croplist:
    module_=importlib.import_module(f'projects.{crop_temp}')
    project_conf_=module_.project_conf()

    #from projects.FibreCropsNes import project_conf
    
    #import warnings
    #warnings.filterwarnings("ignore")
    
    #irrigation scenarios include 5
    irr5list=['Soybean','Banana','ArabicaCoffee','RobustaCoffee','Coconut','Olive','Cocoa','Cotton','Bean',\
              'Oilpalm','Lentils','Chickpea','Plantains','Cowpea','PigeonPea','Tea','Tobacco','Tomato',\
            'Onion','Pea','Apple','Avocado','VegetablesTuber','VegetablesLeafy',\
                'Vegetables','Sheanuts','Rubber','Pulses','Pepper','Oranges','Oilseeds','Grapes',\
                    'FruitTropical','FibreCropsNes','Cashew','BerriesNes']
    
    #calibration soil fertility and export calibrated parameters
    #project_conf_=project_conf()
    irr5=False
    if project_conf_.project_name in irr5list:
        irr5=True
    irr_rice=False
    if project_conf_.project_name in ["Rice1","Rice2"]:
        irr_rice=True
    
    #Yield under soil fertility stress for low-input rainfed, real
    project_conf_.virtual_irrigation='Lowinput'
    project_conf_.scenarios = [1,3]
    
    project_conf_.soil_fertility=1
    
    model = acea_model(project_conf_)
    
    model.calib_soilfertility()
    model.RasterParametersSF()
    
    project_conf_.scenarios = [1,2]
    model = acea_model(project_conf_)
    
    print(project_conf_.virtual_irrigation)
    model.Run()
    model.Rasterise()
    
    
    #Yield under soil fertility stress if no water stress for rainfed locations
    project_conf_.virtual_irrigation='Lowvirt'
    project_conf_.scenarios = [1,2]
    project_conf_.soil_fertility=1
    
    model = acea_model(project_conf_)
    
    print(project_conf_.virtual_irrigation)
    model.Run()
    model.Rasterise()
    
    #Yield under no soil fertility stress if no water stress for rainfed locations
    project_conf_.virtual_irrigation='Highvirt'
    project_conf_.scenarios = [1,2]
    project_conf_.soil_fertility=1
    
    model = acea_model(project_conf_)
    
    print(project_conf_.virtual_irrigation)
    model.Run()
    model.Rasterise()
    
    #Yield under no soil fertility stress for high-input rainfed and irrigated
    project_conf_.virtual_irrigation='Highinput'
    project_conf_.scenarios = [1,2,3,4]
    if irr5: project_conf_.scenarios = [1,2,3,4,5]
    if irr_rice: project_conf_.scenarios = [1,2,6]
    project_conf_.soil_fertility=1 
    
    model = acea_model(project_conf_)
    
    print(project_conf_.virtual_irrigation)
    model.Run()
    model.Rasterise()   
    
    
    #Yield under no soil fertility stress if no water stress for rainfed locations, trun off soil fertility stress, not necessary
    project_conf_.virtual_irrigation='HighvirtNoSF'
    project_conf_.scenarios = [1,2]
    project_conf_.soil_fertility=0
    
    model = acea_model(project_conf_)
    
    print(project_conf_.virtual_irrigation)
    #model.Run()
    #model.Rasterise()
    
    #Yield under no soil fertility stress for high-input rainfed and irrigated, trun off soil fertility stress
    #the difference with Olkes is to use regional CCx and HI for highinput, and a different soil texture
    project_conf_.virtual_irrigation='HighinputNoSF'
    project_conf_.scenarios = [1,2,3,4]
    if irr5: project_conf_.scenarios = [1,2,3,4,5]
    if irr_rice: project_conf_.scenarios = [1,2,6]
    project_conf_.soil_fertility=0 
    
    model = acea_model(project_conf_)
    #project_conf_.gricells
    print(project_conf_.virtual_irrigation)
    model.Run()
    model.Rasterise()
    
    
    import numpy as np
    import pandas as pd
    from sqlalchemy import create_engine
    from modules.acea.acea_yieldscaling import read_spam,read_acea,combine_spam_acea,add_farm_size
    engine = create_engine('postgresql://postgres:0220@localhost/postgres')
    import xarray as xr
    
    import warnings
    warnings.filterwarnings("ignore")
    
    #run from here
    project_name=project_conf_.project_name
    
    crop_name=project_conf_.crop_name_4code
    crop_name3=project_conf_.crop_name_short
    crop_code=project_conf_.crop_fao

    #import sys
    
    if crop_code in [15,27]:
        if crop_code==15:
            project_name='Wheat'
            crop_name3='wht'
            
        elif crop_code==27:
            project_name='Rice'
            crop_name3='ric'
        
        print('Please make sure that you merge the splitted crops and then run the post-processing manually from xr import!')
        continue
        
        
    if crop_code in [10000,10001,641]:
        
        print('Not a FAO or SPAM crop!')
        continue
    
    #in the output results, all the yield will be converted to fresh yield as FAO reported
    irr_sys=xr.open_dataset(r'data\acea\harvested_areas\cell_share_irrigation_systems_5arcmin_ss.nc')
    dry_fresh = np.load(r'data\acea\crop_parameters\crop_dry_to_fresh_ratios.npz', allow_pickle=True)['crop_ratios'] 
    dry_fresh_=dry_fresh[dry_fresh[:,0]==crop_code][0][1]
    
    ti=read_spam(crop_name,'i')
    th=read_spam(crop_name,'h')
    tl=read_spam(crop_name,'l')
    ts=read_spam(crop_name,'s')
    
    for year in [2008,2009,2010,2011,2012]:
    
        for item in ['yield']:
            print(year,item,'For scaling')
        
            high_sc1=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            high_sc2=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr_rice:
                high_sc3=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            else:
                high_sc3=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc3_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc4_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            low_sc1=read_acea(r'outputs\{0}\Lowinput\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            low_sc2=read_acea(r'outputs\{0}\Lowinput\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr5:
                high_sc5=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc5_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr5:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+str(year),'{0}-{1}'.format(item,crop_name3),irrigated=2,dry_fresh=dry_fresh_,acea3=high_sc5)
            else:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+str(year),'{0}-{1}'.format(item,crop_name3),irrigated=1,dry_fresh=dry_fresh_)
            th=combine_spam_acea(th,high_sc1,high_sc2,item+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            tl=combine_spam_acea(tl,low_sc1,low_sc2,item+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            ts=combine_spam_acea(ts,low_sc1,low_sc2,item+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            
            
            high_sc1=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            high_sc2=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr_rice:
                high_sc3=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            else:
                high_sc3=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc3_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc4_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr5:
                high_sc5=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc5_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
    
            if irr5:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),irrigated=2,dry_fresh=dry_fresh_,acea3=high_sc5)
            else:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),irrigated=1,dry_fresh=dry_fresh_)
            
            th=combine_spam_acea(th,high_sc1,high_sc2,item+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            tl=combine_spam_acea(tl,high_sc1,high_sc2,item+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            ts=combine_spam_acea(ts,high_sc1,high_sc2,item+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
                
    #add high low f3
    f3low=xr.open_dataset(r'outputs\{0}\rastersSF\{1}_{2}_rainfed_f3low.nc'.format(project_name,project_name,crop_code))
    f3high=xr.open_dataset(r'outputs\{0}\rastersSF\{1}_{2}_rainfed_f3high.nc'.format(project_name,project_name,crop_code))
    f3high_irr=xr.open_dataset(r'outputs\{0}\rastersSF\{1}_{2}_irr_f3high.nc'.format(project_name,project_name,crop_code))
    
    ti=combine_spam_acea(ti,f3high_irr,f3high_irr,'f3','f3high')
    th=combine_spam_acea(th,f3high,f3high,'f3','f3high')
    tl=combine_spam_acea(tl,f3low,f3low,'f3','f3low')
    ts=combine_spam_acea(ts,f3low,f3low,'f3','f3low')
    
    
    for item in [ti,th,tl,ts]:
        for year in [2008,2009,2010,2011,2012]:
            item['production'+str(year)]=item[crop_name]*item['yield'+str(year)]
            
            item['production'+str(year)+'NOSF']=item[crop_name]*item['yield'+str(year)+'NOSF']
            item['production'+str(year)+'f3']=item[crop_name]*item['yield'+str(year)]*item['f3']
            item['NOSF'+str(year)]=item['yield'+str(year)+'NOSF']==item['yield'+str(year)]
    
    for item in [ti,th,tl,ts]:
        item=item.dropna()
        
    ti_=ti.groupby(['iso3','NOSF'+str(2010)]).sum()
    th_=th.groupby(['iso3','NOSF'+str(2010)]).sum()
    tl_=tl.groupby(['iso3','NOSF'+str(2010)]).sum()
    ts_=ts.groupby(['iso3','NOSF'+str(2010)]).sum()
    
    
    country=set([item[0] for item in ti_.index]+[item[0] for item in th_.index]+\
                [item[0] for item in tl_.index]+[item[0] for item in ts_.index])
    yield_=pd.DataFrame(index=list(country))
    
    for ctry in country:
        area=0
        production=0
        area1=0
        production1=0
        production2=0
        for year in [2008,2009,2010,2011,2012]:
    
            for item in [ti_,th_,ts_,tl_]:
                if ctry in item.index or (ctry,True) in item.index or (ctry,False) in item.index:
    
                    if (ctry,True) in item.index:
                        area+=item.loc[(ctry,True),crop_name]
                        production+=item.loc[(ctry,True),'production'+str(year)]  
                        production2+=item.loc[(ctry,True),'production'+str(year)+'f3']
                    if (ctry,False) in item.index:
                        area1+=item.loc[(ctry,False),crop_name]
                        production1+=item.loc[(ctry,False),'production'+str(year)]
                        production2+=item.loc[(ctry,False),'production'+str(year)+'f3']
        
        if area+area1>0:
    
            yield_.loc[ctry,'yield']=(production+production1)/(area+area1)
            yield_.loc[ctry,'yieldf3']=production2/(area+area1)
            yield_.loc[ctry,'area_SPAM']=area+area1
            yield_.loc[ctry,'pro_NOSF']=production
            yield_.loc[ctry,'pro_SF']=production1
                
    
    fao_=pd.read_csv(r'data\acea\harvested_areas\faostat_production_2008_2012.csv')
    fao_=fao_.loc[fao_['Item Code (FAO)']==project_conf_.crop_fao,:]
    fao_=fao_.set_index('Area Code (ISO3)')
    countries=fao_.index.unique()
    fao_=fao_.reset_index().groupby(['Area Code (ISO3)','Element']).sum(numeric_only=True)
    
    ctry_code=pd.read_table(r'data\acea\grid\regions_countries.txt')
    ctry_code=ctry_code.set_index('ISO3')
    
    for name in countries:
        if name in ctry_code.index:
            if name in yield_.index:
                
                yield_.loc[name,'ISO3']=name
                yield_.loc[name,'Name']=ctry_code.loc[name,'Short name']
                yield_.loc[name,'faostat']=ctry_code.loc[name,'FAOSTAT']
                yield_.loc[name,'yield']=yield_.loc[name,'yield']
                yield_.loc[name,'yieldf3']=yield_.loc[name,'yieldf3']
                yield_.loc[name,'area_SPAM']=yield_.loc[name,'area_SPAM']
                try:
                    yield_.loc[name,'faoyield']=fao_.loc[(name,'Production'),'Value']/fao_.loc[(name,'Area harvested'),'Value']
                    yield_.loc[name,'area_FAO']=fao_.loc[(name,'Area harvested'),'Value']
                    yield_.loc[name,'scaling']=yield_.loc[name,'faoyield']/yield_.loc[name,'yield']
                    yield_.loc[name,'scalingf3']=yield_.loc[name,'faoyield']/yield_.loc[name,'yieldf3']
                    yield_.loc[name,'scaling_area']=yield_.loc[name,'area_FAO']/yield_.loc[name,'area_SPAM']
                    
                except:
                    pass
                
                
    yield_.to_csv(r'outputs\{0}\{1}_scaling.csv'.format(project_name,crop_name),index=True)
    
    ti=read_spam(project_conf_.crop_name_4code,'i')
    th=read_spam(project_conf_.crop_name_4code,'h')
    tl=read_spam(project_conf_.crop_name_4code,'l')
    ts=read_spam(project_conf_.crop_name_4code,'s')
    
    
    irr_sys=xr.open_dataset(r'data\acea\harvested_areas\cell_share_irrigation_systems_5arcmin_ss.nc')
    
    
    for year in [2008,2009,2010,2011,2012]:
    
        for item in ['yield','aet','cwu_green','transp']:
            
            print(year,item,'Actual')
        
            high_sc1=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            high_sc2=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr_rice:
                high_sc3=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            else:
                high_sc3=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc3_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc4_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            low_sc1=read_acea(r'outputs\{0}\Lowinput\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            low_sc2=read_acea(r'outputs\{0}\Lowinput\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr5:
                high_sc5=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc5_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                
            if irr5:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+'_'+str(year),'{0}-{1}'.format(item,crop_name3),irrigated=2,dry_fresh=dry_fresh_,acea3=high_sc5)
            else:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+'_'+str(year),'{0}-{1}'.format(item,crop_name3),irrigated=1,dry_fresh=dry_fresh_)
    
            th=combine_spam_acea(th,high_sc1,high_sc2,item+'_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            tl=combine_spam_acea(tl,low_sc1,low_sc2,item+'_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            ts=combine_spam_acea(ts,low_sc1,low_sc2,item+'_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            
            tl=combine_spam_acea(tl,high_sc1,high_sc2,item+'_hrain_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            ts=combine_spam_acea(ts,high_sc1,high_sc2,item+'_hrain_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            
            high_sc1=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            high_sc2=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr_rice:
                high_sc3=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc6_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            else:
                high_sc3=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc3_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                high_sc4=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc4_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            if irr5:
                high_sc5=read_acea(r'outputs\{0}\HighinputNoSF\rasters\acea_30arc_sc5_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
                
            if irr5:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+'_'+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),irrigated=2,dry_fresh=dry_fresh_,acea3=high_sc5)
            else:
                ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,item+'_'+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),irrigated=1,dry_fresh=dry_fresh_)
    
            th=combine_spam_acea(th,high_sc1,high_sc2,item+'_'+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            tl=combine_spam_acea(tl,high_sc1,high_sc2,item+'_'+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            ts=combine_spam_acea(ts,high_sc1,high_sc2,item+'_'+str(year)+'NOSF','{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            
        high_sc2=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc2_{1}_cwu_blue_cr_global_annual_2008_2012.nc'.format(project_name,crop_name3),time=year)
        
        if irr_rice:
            high_sc3=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc6_{1}_cwu_blue_ir_global_annual_2008_2012.nc'.format(project_name,crop_name3),time=year)
            high_sc4=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc6_{1}_cwu_blue_ir_global_annual_2008_2012.nc'.format(project_name,crop_name3),time=year)
        
        else:
            high_sc3=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc3_{1}_cwu_blue_ir_global_annual_2008_2012.nc'.format(project_name,crop_name3),time=year)
            high_sc4=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc4_{1}_cwu_blue_ir_global_annual_2008_2012.nc'.format(project_name,crop_name3),time=year)
        
        low_sc2=read_acea(r'outputs\{0}\Lowinput\rasters\acea_30arc_sc2_{1}_cwu_blue_cr_global_annual_2008_2012.nc'.format(project_name,crop_name3),time=year)
        
        if irr5:
            high_sc5=read_acea(r'outputs\{0}\Highinput\rasters\acea_30arc_sc5_{1}_cwu_blue_ir_global_annual_2008_2012.nc'.format(project_name,crop_name3),time=year)
            
        if irr5:
            ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,'cwu_blue_ir'+'_'+str(year),'cwu_blue_ir-{0}'.format(crop_name3),irrigated=2,acea3=high_sc5)
        else:
            ti,irr=combine_spam_acea(ti,high_sc3,high_sc4,'cwu_blue_ir'+'_'+str(year),'cwu_blue_ir-{0}'.format(crop_name3),irrigated=1)
    
        th=combine_spam_acea(th,high_sc2,high_sc2,'cwu_blue_cr'+'_'+str(year),'cwu_blue_cr-{0}'.format(crop_name3))
        tl=combine_spam_acea(tl,low_sc2,low_sc2,'cwu_blue_cr'+'_'+str(year),'cwu_blue_cr-{0}'.format(crop_name3))
        ts=combine_spam_acea(ts,low_sc2,low_sc2,'cwu_blue_cr'+'_'+str(year),'cwu_blue_cr-{0}'.format(crop_name3))
                
        for item in ['yield','aet','cwu_green','transp']:
            print(year,item,'Scenario')
            
            high_sc1=read_acea(r'outputs\{0}\Lowvirt\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            high_sc2=read_acea(r'outputs\{0}\Lowvirt\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            tl=combine_spam_acea(tl,high_sc1,high_sc2,item+'_lv_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            ts=combine_spam_acea(ts,high_sc1,high_sc2,item+'_lv_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
    
            high_sc1=read_acea(r'outputs\{0}\Highvirt\rasters\acea_30arc_sc1_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            high_sc2=read_acea(r'outputs\{0}\Highvirt\rasters\acea_30arc_sc2_{1}_{2}_global_annual_2008_2012.nc'.format(project_name,crop_name3,item),time=year)
            
            tl=combine_spam_acea(tl,high_sc1,high_sc2,item+'_hv_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            ts=combine_spam_acea(ts,high_sc1,high_sc2,item+'_hv_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
            th=combine_spam_acea(th,high_sc1,high_sc2,item+'_hv_'+str(year),'{0}-{1}'.format(item,crop_name3),dry_fresh=dry_fresh_)
    
        ti['band1']=irr['band1']
        ti['band2']=irr['band2']
        
        
    #add high low f3
    f3low=xr.open_dataset(r'outputs\{0}\rastersSF\{1}_{2}_rainfed_f3low.nc'.format(project_name,project_name,crop_code))
    f3high=xr.open_dataset(r'outputs\{0}\rastersSF\{1}_{2}_rainfed_f3high.nc'.format(project_name,project_name,crop_code))
    f3high_irr=xr.open_dataset(r'outputs\{0}\rastersSF\{1}_{2}_irr_f3high.nc'.format(project_name,project_name,crop_code))
    
    ti=combine_spam_acea(ti,f3high_irr,f3high_irr,'f3high_irr','f3high')
    th=combine_spam_acea(th,f3high,f3high,'f3high_rainfed','f3high')
    tl=combine_spam_acea(tl,f3low,f3low,'f3low_rainfed','f3low')
    ts=combine_spam_acea(ts,f3low,f3low,'f3low_rainfed','f3low')
    
    
    ti.to_csv(r'outputs\{0}\{1}_ti.csv'.format(project_name,crop_name),index=False)
    th.to_csv(r'outputs\{0}\{1}_th.csv'.format(project_name,crop_name),index=False)
    tl.to_csv(r'outputs\{0}\{1}_tl.csv'.format(project_name,crop_name),index=False)
    ts.to_csv(r'outputs\{0}\{1}_ts.csv'.format(project_name,crop_name),index=False)
    
    
    farm_size_path=r'D:\ACEA_v2_2\data_for_analysis\ESSDmap\01Farm-size- and Crop-specific map\csv\perFarmSize_perCrop_perFarmingSystem_separatedFile\SPAMbased_downscaled_map'
    
    if crop_code in [667]:
        pass
    else:
        fs_ti=pd.read_csv(farm_size_path+r'\spambased_farmsizemap_{0}_i.csv'.format(crop_name))
        ti_fs=add_farm_size(fs_ti,ti)
        ti_fs.to_csv(r'outputs\{0}\{1}_ti_fs.csv'.format(project_name,crop_name),index=False)
        
    fs_th=pd.read_csv(farm_size_path+r'\spambased_farmsizemap_{0}_h.csv'.format(crop_name))
    fs_tl=pd.read_csv(farm_size_path+r'\spambased_farmsizemap_{0}_l.csv'.format(crop_name))
    fs_ts=pd.read_csv(farm_size_path+r'\spambased_farmsizemap_{0}_s.csv'.format(crop_name))
    
    th_fs=add_farm_size(fs_th,th)
    tl_fs=add_farm_size(fs_tl,tl)
    ts_fs=add_farm_size(fs_ts,ts)
    
    th_fs.to_csv(r'outputs\{0}\{1}_th_fs.csv'.format(project_name,crop_name),index=False)
    tl_fs.to_csv(r'outputs\{0}\{1}_tl_fs.csv'.format(project_name,crop_name),index=False)
    ts_fs.to_csv(r'outputs\{0}\{1}_ts_fs.csv'.format(project_name,crop_name),index=False)
    
    ti=[];th=[];tl=[];ts=[]
    ti_fs=[];th_fs=[];tl_fs=[];ts_fs=[]
    
    if crop_code in [79,656]:
        
        print('At least 2 simulated crops belonging to one FAO crop, please update the scaling factor in the seperate *.py!')
    
    



