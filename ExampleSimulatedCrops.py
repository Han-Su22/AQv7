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


croplist=['WinterWheatEXP']
    
    
for crop_temp in croplist:
    module_=importlib.import_module(f'projects.{crop_temp}')
    project_conf_=module_.project_conf()
    project_conf_.use_soilgrids=0
    
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



