# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 15:12:52 2022

@author: MialykO
"""

class project_conf:
    # General
    crop_name = 'Chestnut' # Full name of a crop
    crop_name_short = 'csn' # Crop's code name (3-letter format)
    scenarios = [1,2,3,4,5]  # Scenarios to consider
    crop_fao = 220 # FAO crop code 
    
    # Deriving from
    cor_gen_crop_fao = 515 # Corresponding generic FAO crop code 
    cor_gen_crop_name_short = 'apl' # Corresponding generic crop's code name (3-letter format)
    
    # Other required info
    project_name = crop_name
    climate_name = 'gswp3-w5e5_obsclim/gswp3-w5e5_obsclim_01-01-1981_31-12-2019' # Name of climate file
    clock_start = '1988/01/01' # First day of simulation (incl)
    clock_end = '2019/12/31' # Last day of simulation (incl)
    spinup = 2 # Number of spinup growing seasons to generate soil moisture (clipped in post-processing)
    
    #needed for soil fertility stress
    crop_name_4code='rest'#SPAM crop
