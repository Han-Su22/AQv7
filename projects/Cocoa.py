# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:51:00 2021

@author: SuH, MialykO
"""

class project_conf:
    
    # General
    project_name = "Cocoa"
    crop_model = 'AquaCrop_v6' # Only this one available for now
    scenarios = [1,2,3,4,5]  # Scenarios to consider
    multi_core = True
    CPUs = 7
    
    # Scope
    resolution = 1 # Spatial resolution
    gridcells = 'world' # Cells to run
    clock_start = '2006/01/01' # First day of simulation (incl)
    clock_end = '2012/12/31' # Last day of simulation (incl)
    spinup = 2 # Number of spinup growing seasons to generate soil moisture (clipped in post-processing)
    real_cropland = True
    landuse = 'spam2010'
        
    # Crop settings
    crop_name = 'Cocoa' # Full name of a crop
    crop_name_short = 'coc' # Crop's code name from GGCMI
    crop_fao = 661 # Crop's code from FAO
    crop_gs_increase = 0 # % to prolong the growing season length if maturity is not reached
    crop_perennial = True
    crop_phenology = False
    
    # Climate settings
    climate_name = 'gswp3-w5e5_obsclim/gswp3-w5e5_obsclim_01-01-1981_31-12-2019' # Name of climate file
    climate_start = 1981; climate_end = 2019
    co2_name = 'GlobalHistoricalCO2_NOAA_1980_2019' # Name of CO2 file
    
    # Soil settings
    soil_dz = [.1, .1, .1, .3, .4, .6, .7, .7] # Soil compartments [m], 3m in total
    init_wc = 50 # Initital soil moisture [%]
    gw_max_level = -1. # Limit shallow groundwater to this number (only for scenarios with groundwater)
    gw_min_level = -3. # Skip the cells with groundwater levels below this (only for scenarios with groundwater
    
    # Field management settings
    off_season = True # Consider or not fallow period (affects soil moisture)
    irr_thresholds = [50]*4 # % of TAW when to trigger irrigation events
        # 100 = soil is always at FC within the root zone (not efficient), so optimal range is 30-60
        # Stages: 1. Initial (Germination), 2. Crop Development (reaching CCmax), 3. Mid-season (Flowering and yield-formation), 4. Late-season (senescence)
        # Stages 1 and 4 are generally less water-stress sensetive, so max irrigation should be applied for stages 2 and 3
    bunds = False # Surface bunds are present
    bunds_dz = 0.3 # Bund height(0.3m is default)
    mulching = False # Apply mulching
    mulching_area = 0.9 # Fraction of area of soil surface covered by mulches (1 is max, the more the better)
    mulching_factor = 0.8 # Soil evaporation adjustment factor due to effect of mulches (1 is max, the more the better)
    
    # Result management
    project_save_annual = True
    project_save_daily = True
    
    # Additional
    rerun_simulated_cells = True
    
    #needed for soil fertility stress
    crop_name_4code='coco'#SPAM crop
    #C:\Users\SuH\OneDrive - Universiteit Twente\WaterResearch\AquaCrop\GAEZ\Crop_table.xlsx
    crop_gaez=['Cocoa','Cocoa cumoun','Cocoa hybrid']#only useful for soil fertility calibration
    #a*(1-np.exp(b*LAI))**c
    LAI_CC_a=1
    LAI_CC_b=-0.6
    LAI_CC_c=1
    get_GAEZ_ccx_hi=True
    #these will be changed in the run file
    virtual_irrigation='Highvirt' #Virt or Real or noSF with soil_fertility=0, add virtual sprinklers irrigation to rainfed agriculture, only irrigation changes, only yield is meaningful
    soil_fertility=0#0 close; 1 simulation with soil fertility stress. For the calibration, data will be read automatically when calib_soilfertility()
    tuned=0#tune soil fertility stress according to FAOSTAT, will resimulate with new tuned soil fertility