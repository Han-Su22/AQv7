# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 15:37:03 2021

@author: MialykO
"""
import numpy as np
import pandas as pd

class CropParameters:
    "Get Crop Class with cell-specific data"
    
    def __init__(self, conf, cell, irrigated, crop_params_rf, crop_params_ir, planting_date, harvest_date=None, PlantPop=0,**kwargs):
        
        # General parameters
        crop_code = conf.crop_fao        
        self.GDDmethod = 3 # Growing degree day calculation method
        self.fshape_b = 13.8135 # Shape factor describing the reduction in biomass production for insufficient growing degree days
        self.PctZmin = 70 # Initial percentage of minimum effective rooting depth
        self.fshape_ex = -6 # Shape factor describing the effects of water stress on root expansion
        self.ETadj = 1 # Adjustment to water stress thresholds depending on daily ET0 (0 = No, 1 = Yes)
        self.beta = 12 # Reduction (%) to p_lo3 when early canopy senescence is triggered
        self.a_Tr = 1 # Exponent parameter for adjustment of Kcx once senescence is triggered
        self.GermThr = 0.2 # Proportion of total water storage needed for crop to germinate
        self.CCmin = 0.05 # Minimum canopy size below which yield formation cannot occur
        self.MaxFlowPct = 100/3 # Proportion of total flowering time (%) at which peak flowering occurs
        self.HIini = 0.01 # Initial harvest index
        self.bsted = 0.000138 # WP co2 adjustment parameter given by Steduto et al. 2007
        self.bface = 0.001165 # WP co2 adjustment parameter given by FACE experiments
        self.Aer = 5 # Vol (%) below saturation at which stress begins to occur due to deficient aeration
        self.LagAer = 3 # Number of days lag before aeration stress affects crop growth
        self.phenology_calibration = False
        #soil fertility stress, default values would lead to no soil fertility stress
        self.Ksccx=1#0-1, 1 for no stress
        self.Ksexpf=1#0-1, 1 for no stress
        self.Kswp=1#0-1, 1 for no stress
        self.fcdecline=0#0-0.01, 0 for no stress
        self.sfertstress=0#0-1, 0 for no stress, the total soil fertility stress, as input for calibration, it is 1-realtiveBio
        self.TR_ET0_fertstress = 0
        self.CGC_CD=-1
        
        #for soil fertility stress calibration
        self.need_calib=0 #1 yes,default 1 output; 2 all possibilities
        self.RelativeBio=1#0-1, 1 for no stress
        self.Ksccx_in=1#0-1, 1 for no stress
        self.fcdecline_in=0 #0 small, 1 medium, 2 large
        self.Ksccx_es=np.zeros(10000)+1
        self.Ksexpf_es=np.zeros(10000)+1
        self.Kswp_es=np.zeros(10000)+1
        self.fcdecline_es=np.zeros(10000)
        self.sf_es=np.zeros(10000)
        self.relbio_es=np.zeros(10000)+1
        self.Bio_top=np.zeros(10000)

        
        # Cell-specific parameters
        if irrigated: data = crop_params_ir
        else: data = crop_params_rf
        
        rowy30 = int(cell['rowy30']); colx30 = int(cell['colx30'])
        try:
            tmax = int(data['temp_max_avg'][rowy30, colx30])
            tmin = int(data['temp_min_avg'][rowy30, colx30])
        except:
            raise ValueError('Masked cell, no crop parameters')
        
        # Phenology
        self.planting_date = planting_date # Planting Date (mm/dd)
        self.harvest_date = harvest_date # Latest Harvest Date (mm/dd)
        self.crop_gs_increase = conf.crop_gs_increase # Possible increase of season if maturity is not reached
        self.crop_perennial = conf.crop_perennial # Check if crop is perennial or annual
        
        if conf.crop_phenology in ['transient','average']:
            self._Emergence = int(data['gdd_emergence'][rowy30, colx30]) # Growing degree/Calendar days from sowing to emergence/transplant recovery
            self._MaxRooting = int(data['gdd_max_root'][rowy30, colx30]) # Growing degree/Calendar days from sowing to maximum rooting
            self._Senescence = int(data['gdd_senescence'][rowy30, colx30]) # Growing degree/Calendar days from sowing to senescence
            self._Maturity = int(data['gdd_maturity'][rowy30, colx30]) # Growing degree/Calendar days from sowing to maturity
            self._HIstart = int(data['gdd_yield_form'][rowy30, colx30]) # Growing degree/Calendar days from sowing to start of yield formation
            self._Flowering = int(data['gdd_duration_flowering'][rowy30, colx30]) # Duration of flowering in growing degree/calendar days (-999 for non-fruit/grain crops)
            self._YldForm = int(data['gdd_duration_yield_form'][rowy30, colx30]) # Duration of yield formation in growing degree/calendar days
            self._CDC = float(data['cdc'][rowy30, colx30]) # Canopy decline coefficient (fraction per GDD/calendar day)
            self._CGC = float(data['cgc'][rowy30, colx30]) # Canopy growth coefficient (fraction per GDD)
            
        else:
            self.Emergence = int(data['gdd_emergence'][rowy30, colx30]) # Growing degree/Calendar days from sowing to emergence/transplant recovery
            self.MaxRooting = int(data['gdd_max_root'][rowy30, colx30]) # Growing degree/Calendar days from sowing to maximum rooting
            self.Senescence = int(data['gdd_senescence'][rowy30, colx30]) # Growing degree/Calendar days from sowing to senescence
            self.Maturity = int(data['gdd_maturity'][rowy30, colx30]) # Growing degree/Calendar days from sowing to maturity
            self.HIstart = int(data['gdd_yield_form'][rowy30, colx30]) # Growing degree/Calendar days from sowing to start of yield formation
            self.Flowering = int(data['gdd_duration_flowering'][rowy30, colx30]) # Duration of flowering in growing degree/calendar days (-999 for non-fruit/grain crops)
            self.YldForm = int(data['gdd_duration_yield_form'][rowy30, colx30]) # Duration of yield formation in growing degree/calendar days
            self.CDC = float(data['cdc'][rowy30, colx30]) # Canopy decline coefficient (fraction per GDD/calendar day)
            self.CGC = float(data['cgc'][rowy30, colx30]) # Canopy growth coefficient (fraction per GDD)
            
        if crop_code == 56:
            self.Name = 'Maize'
            self.CropType = 3 # Crop Type (1 = Leafy vegetable, 2 = Root/tuber, 3 = Fruit/grain)
            self.PlantMethod = 1 # Planting method (0 = Transplanted, 1 =  Sown)
            self.CalendarType = 2 # Calendar Type (1 = Calendar days, 2 = Growing degree days)
            self.SwitchGDD = 0 # Convert calendar to GDD mode if inputs are given in calendar days (0 = No; 1 = Yes)
         
            if tmax >= 40: # Adjust for extreme environments
                self.Tmax_up = tmax*1 # Maximum air temperature (degC) above which pollination begins to fail
                self.Tmax_lo = tmax + 5 # Maximum air temperature (degC) at which pollination completely fails
            else:
                self.Tmax_up = 40; self.Tmax_lo = 45 # Default
                
            if tmin <= 10:
                self.Tmin_up = tmin*1 # Minimum air temperature (degC) below which pollination begins to fail (default)
                self.Tmin_lo = tmin - 5 # Minimum air temperature (degC) at which pollination completely fails (default)
            else:
                self.Tmin_up = 10; self.Tmin_lo = 5 # Default
                
            self.Tbase = 8 # Base temperature (degC) below which growth does not progress
            self.Tupp = 30 # Upper temperature (degC) above which crop development no longer increases
            self.PolHeatStress = 1 # Pollination affected by heat stress (0 = No, 1 = Yes)
            self.PolColdStress = 1 # Pollination affected by cold stress (0 = No, 1 = Yes)
            self.TrColdStress = 1 # Transpiration affected by cold temperature stress (0 = No, 1 = Yes)
            self.GDD_up = 12 # Minimum growing degree days (degC/day) required for full crop transpiration potential
            self.GDD_lo = 0 # Growing degree days (degC/day) at which no crop transpiration occurs
            self.Zmin = 0.3 # Minimum effective rooting depth (m)
            self.Zmax = 2.3 # Maximum rooting depth (m)
            self.fshape_r = 1.3 # Shape factor describing root expansion
            self.SxTopQ = 0.0104 # Maximum root water extraction at top of the root zone (m3/m3/day)
            self.SxBotQ = 0.0026 # Maximum root water extraction at the bottom of the root zone (m3/m3/day)
            self.SeedSize = 6.5 # Soil surface area (cm2) covered by an individual seedling at 90% emergence
            self.PlantPop = 75_000 if PlantPop == 0 else PlantPop # Number of plants per hectare
            self.CCx = 0.96 # Maximum canopy cover (fraction of soil cover)
            self.Kcb = 1.05 # Crop coefficient when canopy growth is complete but prior to senescence
            self.fage = 0.3 #  Decline of crop coefficient due to ageing (%/day)
            self.WP = 33.7 # Water productivity normalized for ET0 and C02 (g/m2)
            self.WPy = 100 # Adjustment of water productivity in yield formation stage (% of WP)
            self.fsink = 0.5 # Crop performance under elevated atmospheric CO2 concentration (%/100)
            self.HI0 = 0.48 # Reference harvest index
            self.dHI_pre = 0 # Possible increase of harvest index due to water stress before flowering (%)
            self.a_HI = 7 # Coefficient describing positive impact on harvest index of restricted vegetative growth during yield formation
            self.b_HI = 3 # Coefficient describing negative impact on harvest index of stomatal closure during yield formation
            self.dHI0 = 15 # Maximum allowable increase of harvest index above reference value
            self.Determinant = 1 # Crop Determinancy (0 = Indeterminant, 1 = Determinant)
            self.exc = 50 # Excess of potential fruits
            self.p_up1 = 0.14 # Upper soil water depletion threshold for water stress effects on affect canopy expansion
            self.p_up2 = 0.69 # Upper soil water depletion threshold for water stress effects on canopy stomatal control
            self.p_up3 = 0.69 # Upper soil water depletion threshold for water stress effects on canopy senescence
            self.p_up4 = 0.8 # Upper soil water depletion threshold for water stress effects on canopy pollination
            self.p_lo1 = 0.72 # Lower soil water depletion threshold for water stress effects on canopy expansion
            self.p_lo2 = 1 # Lower soil water depletion threshold for water stress effects on canopy stomatal control
            self.p_lo3 = 1 #  Lower soil water depletion threshold for water stress effects on canopy senescence
            self.p_lo4 = 1 # Lower soil water depletion threshold for water stress effects on canopy pollination
            self.fshape_w1 = 2.9 # Shape factor describing water stress effects on canopy expansion
            self.fshape_w2 = 6 # Shape factor describing water stress effects on stomatal control
            self.fshape_w3 = 2.7 # Shape factor describing water stress effects on canopy senescence
            self.fshape_w4 = 1 # Shape factor describing water stress effects on pollination
        
        elif crop_code == 15:
            self.Name = 'Wheat'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 0; self.Tupp = 26
            
            if tmax >= 35: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 35; self.Tmax_lo = 40
                
            if tmin <= 5: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 5; self.Tmin_lo = 0
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 14; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = 0.028; self.SxBotQ = 0.008
            self.SeedSize = 1.5; self.PlantPop = 4_500_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.96; self.Kcb = 1.1; self.fage = 0.15
            self.WP = 15; self.WPy = 100; self.fsink = 0.5
            self.HI0 = 0.48; self.dHI_pre = 5; self.a_HI = 10; self.b_HI = 7; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.2; self.p_up2 = 0.65; self.p_up3 = 0.7;  self.p_up4 = 0.8
            self.p_lo1 = 0.65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 5.; self.fshape_w2 = 2.5; self.fshape_w3 = 2.5; self.fshape_w4 = 1
            
        elif crop_code == 27:
            self.Name = 'Rice'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 8; self.Tupp = 30
            
            if tmax >= 35: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 35; self.Tmax_lo = 40
            if tmin <= 5: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 10; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = .5; self.fshape_r = 2.5
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 6; self.PlantPop = 1_000_000 if PlantPop == 0 else PlantPop
            self.CCx = .95; self.Kcb = 1.2; self.fage = 0.15
            self.WP = 19; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .43; self.dHI_pre = 0; self.a_HI = 10; self.b_HI = 7; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0; self.p_up2 = 0.5; self.p_up3 = 0.55;  self.p_up4 = 0.75
            self.p_lo1 = 0.4; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
            
            # ONLY FOR RICE
            self.Aer = 0 # Vol (%) below saturation at which stress begins to occur due to deficient aeration
            self.LagAer = 999 # Number of days lag before aeration stress affects crop growth
            
        elif crop_code == 44:
            self.Name = 'Barley'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 2; self.Tupp = 28
            
            if tmax >= 35: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 35; self.Tmax_lo = 40
                
            if tmin <= 5: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 5; self.Tmin_lo = 0
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 14; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 1.3; self.fshape_r = 1.5
            self.SxTopQ = 0.019; self.SxBotQ = 0.006
            self.SeedSize = 1.5; self.PlantPop = 1_500_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.8; self.Kcb = 1.1; self.fage = .15
            self.WP = 15; self.WPy = 100; self.fsink = 0.5
            self.HI0 = 0.4; self.dHI_pre = 5; self.a_HI = 10; self.b_HI = 5; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.2; self.p_up2 = 0.6; self.p_up3 = 0.55;  self.p_up4 = 0.85
            self.p_lo1 = 0.65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 75:
            self.Name = 'Oats'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 0; self.Tupp = 30
            
            if tmax >= 35: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 35; self.Tmax_lo = 40
                
            if tmin <= 5: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 5; self.Tmin_lo = 0
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 14; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 1.3; self.fshape_r = 1.5
            self.SxTopQ = .019; self.SxBotQ = .006
            self.SeedSize = 1.5; self.PlantPop = 1_600_000 if PlantPop == 0 else PlantPop
            self.CCx = .98; self.Kcb = 1.1; self.fage = .15
            self.WP = 12.4; self.WPy = 100; self.fsink = .5
            self.HI0 = .5; self.dHI_pre = 5; self.a_HI = 10; self.b_HI = 5; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = .45; self.p_up2 = .8; self.p_up3 = .75;  self.p_up4 = .9
            self.p_lo1 = .8; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 79:
            self.Name = 'Millet'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 33
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5

            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 12; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.8; self.fshape_r = 1.3
            self.SxTopQ = .051; self.SxBotQ = .013
            self.SeedSize = 5; self.PlantPop = 180_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = 1.07; self.fage = 0.3
            self.WP = 33.7; self.WPy = 100; self.fsink = 0.5
            self.HI0 = 0.52; self.dHI_pre = 4; self.a_HI = 1; self.b_HI = 3; self.dHI0 = 25
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.35; self.p_up2 = 0.75; self.p_up3 = 0.8;  self.p_up4 = 0.8
            self.p_lo1 = 0.7; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 83:
            self.Name = 'Sorghum'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 8; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
                
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 12; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 1.8; self.fshape_r = 1.3
            self.SxTopQ = 0.016; self.SxBotQ = 0.004
            self.SeedSize = 3; self.PlantPop = 200_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.9; self.Kcb = 1.07; self.fage = 0.3
            self.WP = 33.7; self.WPy = 100; self.fsink = 0.5
            self.HI0 = 0.45; self.dHI_pre = 4; self.a_HI = 1; self.b_HI = 3; self.dHI0 = 25
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.15; self.p_up2 = 0.7; self.p_up3 = 0.7;  self.p_up4 = 0.8
            self.p_lo1 = 0.7; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 6; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 108:
            self.Name = 'Cereals'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5

            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 12; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.8; self.fshape_r = 1.3
            self.SxTopQ = .051; self.SxBotQ = .013
            self.SeedSize = 5; self.PlantPop = 180_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 15.; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .5; self.dHI_pre = 4; self.a_HI = 1; self.b_HI = 3; self.dHI0 = 25
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.35; self.p_up2 = 0.75; self.p_up3 = 0.8;  self.p_up4 = 0.8
            self.p_lo1 = 0.7; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 116:
            self.Name = 'Potato'; self.CropType = 2
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 2; self.Tupp = 26
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 1
            self.GDD_up = 7; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = 0.016; self.SxBotQ = 0.004
            self.SeedSize = 10; self.PlantPop = 58_000 if PlantPop == 0 else PlantPop
            self.CCx = .96; self.Kcb = 1.1; self.fage = 0.15
            self.WP = 18; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .82; self.dHI_pre = 2; self.a_HI = 0; self.b_HI = 10; self.dHI0 = 5
            self.Determinant = 0; self.exc = 100
            self.p_up1 = 0.2; self.p_up2 = 0.55; self.p_up3 = 0.7;  self.p_up4 = 0.9
            self.p_lo1 = 0.6; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 122:
            self.Name = 'Sweet potato'; self.CropType = 2
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 35
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 1
            self.GDD_up = 11; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 1.6; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .013
            self.SeedSize = 15; self.PlantPop = 40_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.94; self.Kcb = 1.1; self.fage = 0.15
            self.WP = 20.; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .55; self.dHI_pre = 2; self.a_HI = 0; self.b_HI = 10; self.dHI0 = 5
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .26; self.p_up2 = .65; self.p_up3 = .69;  self.p_up4 = .85
            self.p_lo1 = .66; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.3; self.fshape_w2 = 3.4; self.fshape_w3 = 2.7; self.fshape_w4 = 1.
        elif crop_code == 125:
            self.Name = 'Cassava'; self.CropType = 2
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 1
            self.GDD_up = 11; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .013
            self.SeedSize = 10; self.PlantPop = 12_500 if PlantPop == 0 else PlantPop
            self.CCx = .88; self.Kcb = .85; self.fage = .5
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .6; self.dHI_pre = 4; self.a_HI = 10; self.b_HI = 4; self.dHI0 = 15
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .2; self.p_up2 = .5; self.p_up3 = .5;  self.p_up4 = .85
            self.p_lo1 = .6; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 136:
            self.Name = 'Taro'; self.CropType = 2
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 1; self.Tupp = 35
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 1
            self.GDD_up = 11; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .013
            self.SeedSize = 25; self.PlantPop = 27_750 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = 1.1; self.fage = .5
            self.WP = 15; self.WPy = 100; self.fsink = .5
            self.HI0 = .8; self.dHI_pre = 4; self.a_HI = 10; self.b_HI = 4; self.dHI0 = 15
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .1; self.p_up2 = .45; self.p_up3 = .45;  self.p_up4 = .85
            self.p_lo1 = .45; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 137:
            self.Name = 'Yams'; self.CropType = 2
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 1
            self.GDD_up = 11; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .013
            self.SeedSize = 10; self.PlantPop = 8_800 if PlantPop == 0 else PlantPop
            self.CCx = .91; self.Kcb = .85; self.fage = .5
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .95; self.dHI_pre = 4; self.a_HI = 10; self.b_HI = 4; self.dHI0 = 15
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .2; self.p_up2 = .5; self.p_up3 = .5;  self.p_up4 = .85
            self.p_lo1 = .6; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 156:
            self.Name = 'Sugarcane'; self.CropType = 1
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 12; self.Tupp = 32
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = -999; self.GDD_lo = -999
            self.Zmin = .3; self.Zmax = 1.8; self.fshape_r = 1.3
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 6.5; self.PlantPop = 140_000 if PlantPop == 0 else PlantPop
            self.CCx = .95; self.Kcb = 1.1; self.fage = .15
            self.WP = 30; self.WPy = 100; self.fsink = .5
            self.HI0 = .35; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .6;  self.p_up4 = .9
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
            self.Aer = 0 # Vol (%) below saturation at which stress begins to occur due to deficient aeration
            self.LagAer = 999 # Number of days lag before aeration stress affects crop growth
        elif crop_code == 157:
            self.Name = 'Sugar beet'; self.CropType = 2
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 3; self.Tupp = 25
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 1
            self.GDD_up = 9; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 1; self.PlantPop = 100_000 if PlantPop == 0 else PlantPop
            self.CCx = .98; self.Kcb = 1.15; self.fage = .15
            self.WP = 18; self.WPy = 100; self.fsink = .5
            self.HI0 = .75; self.dHI_pre = 0; self.a_HI = 6; self.b_HI = 0; self.dHI0 = 20
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .2; self.p_up2 = .65; self.p_up3 = .75;  self.p_up4 = 1
            self.p_lo1 = .6; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 176:
            self.Name = 'Beans'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 9; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 0.3; self.Zmax = 1.7; self.fshape_r = 1.5
            self.SxTopQ = 0.048; self.SxBotQ = 0.012
            self.SeedSize = 10.; self.PlantPop = 131_579 if PlantPop == 0 else PlantPop
            self.CCx = 0.99; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 15; self.WPy = 90; self.fsink = 0.5
            self.HI0 = 0.4; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 1; self.dHI0 = 10
            self.Determinant = 0; self.exc = 50
            self.p_up1 = 0.15; self.p_up2 = 0.6; self.p_up3 = 0.7;  self.p_up4 = 0.88
            self.p_lo1 = 0.65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 2.5; self.fshape_w2 = 3; self.fshape_w3 = 2.5; self.fshape_w4 = 1.
        elif crop_code == 187:
            self.Name = 'Pea'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 27
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = 0.048; self.SxBotQ = 0.012
            self.SeedSize = 5.; self.PlantPop = 810_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.95; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 12; self.WPy = 90; self.fsink = .5
            self.HI0 = .42; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 1; self.dHI0 = 10
            self.Determinant = 0; self.exc = 50
            self.p_up1 = .1; self.p_up2 = .45; self.p_up3 = .45;  self.p_up4 = .8
            self.p_lo1 = .65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3.; self.fshape_w4 = 1.
        elif crop_code == 191:
            self.Name = 'Chickpea'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 8; self.Tupp = 32
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 5.; self.PlantPop = 330_000 if PlantPop == 0 else PlantPop
            self.CCx = .93; self.Kcb = 1.1; self.fage = .3
            self.WP = 16; self.WPy = 90; self.fsink = .5
            self.HI0 = .36; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 1; self.dHI0 = 10
            self.Determinant = 0; self.exc = 50
            self.p_up1 = .3; self.p_up2 = .5; self.p_up3 = .86;  self.p_up4 = 1.
            self.p_lo1 = .65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3.; self.fshape_w4 = 1.
        elif crop_code == 195:
            self.Name = 'Cowpea'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 36
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = 0.048; self.SxBotQ = 0.012
            self.SeedSize = 10.; self.PlantPop = 200_000 if PlantPop == 0 else PlantPop
            self.CCx = .93; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 16; self.WPy = 90; self.fsink = .5
            self.HI0 = .29; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 1; self.dHI0 = 10
            self.Determinant = 0; self.exc = 50
            self.p_up1 = .3; self.p_up2 = .4; self.p_up3 = .6;  self.p_up4 = 0.88
            self.p_lo1 = .65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 2.5; self.fshape_w2 = 3; self.fshape_w3 = 2.5; self.fshape_w4 = 1.
        elif crop_code == 197:
            self.Name = 'Pigeon pea'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = 0.048; self.SxBotQ = 0.012
            self.SeedSize = 5.; self.PlantPop = 810_000 if PlantPop == 0 else PlantPop
            self.CCx = .95; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 15; self.WPy = 90; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 1; self.dHI0 = 10
            self.Determinant = 0; self.exc = 50
            self.p_up1 = .1; self.p_up2 = .45; self.p_up3 = .45;  self.p_up4 = .8
            self.p_lo1 = .65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3.; self.fshape_w4 = 1.
        elif crop_code == 201:
            self.Name = 'Lentils'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = .3; self.Zmax = 1.7; self.fshape_r = 1.5
            self.SxTopQ = 0.048; self.SxBotQ = 0.012
            self.SeedSize = 10.; self.PlantPop = 130_000 if PlantPop == 0 else PlantPop
            self.CCx = .95; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 15; self.WPy = 90; self.fsink = 0.5
            self.HI0 = .4; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 1; self.dHI0 = 10
            self.Determinant = 0; self.exc = 50
            self.p_up1 = 0.15; self.p_up2 = 0.6; self.p_up3 = 0.7;  self.p_up4 = 0.88
            self.p_lo1 = 0.65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 2.5; self.fshape_w2 = 3; self.fshape_w3 = 2.5; self.fshape_w4 = 1.
        elif crop_code == 211:
            self.Name = 'Pulses'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = .3; self.Zmax = 1.7; self.fshape_r = 1.5
            self.SxTopQ = 0.048; self.SxBotQ = 0.012
            self.SeedSize = 10.; self.PlantPop = 130_000 if PlantPop == 0 else PlantPop
            self.CCx = .95; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 15; self.WPy = 90; self.fsink = 0.5
            self.HI0 = .4; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 1; self.dHI0 = 10
            self.Determinant = 0; self.exc = 50
            self.p_up1 = 0.15; self.p_up2 = 0.6; self.p_up3 = 0.7;  self.p_up4 = 0.88
            self.p_lo1 = 0.65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 2.5; self.fshape_w2 = 3; self.fshape_w3 = 2.5; self.fshape_w4 = 1.
        elif crop_code == 217:
            self.Name = 'Cashew'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 35
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 2.; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 230_000 if PlantPop == 0 else PlantPop
            self.CCx = .65; self.Kcb = .8; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .1; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 236:
            self.Name = 'Soybean'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 10; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = 0.012; self.SxBotQ = 0.003
            self.SeedSize = 5.; self.PlantPop = 330_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.98; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 15; self.WPy = 60; self.fsink = 0.5
            self.HI0 = 0.4; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 10
            self.Determinant = 1; self.exc = 50
            self.p_up1 = 0.15; self.p_up2 = 0.6; self.p_up3 = 0.7;  self.p_up4 = 0.85
            self.p_lo1 = 0.65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 242:
            self.Name = 'Groundnut'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 9; self.Tupp = 30
            
            if tmax >= 35: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 35; self.Tmax_lo = 40
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 10; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 3.; self.PlantPop = 88_889 if PlantPop == 0 else PlantPop
            self.CCx = .7; self.Kcb = 1.1; self.fage = .3
            self.WP = 15; self.WPy = 100; self.fsink = .5
            self.HI0 = .24; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 10
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .5; self.p_up2 = .8; self.p_up3 = .9;  self.p_up4 = .95
            self.p_lo1 = .8; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 1; self.fshape_w2 = 2; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 249:
            self.Name = 'Coconut'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 40
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.5; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 390_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = 1.; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 254:
            self.Name = 'Oil palm'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 40
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.5; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 390_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = 1.; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 260:
            self.Name = 'Olives'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 6; self.Tupp = 31
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.7; self.Zmax = 1.7; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 250_000 if PlantPop == 0 else PlantPop
            self.CCx = .5; self.Kcb = .7; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 263:
            self.Name = 'Sheanut'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 40
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.5; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 390_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = 1.; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 267:
            self.Name = 'Sunflower'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 4; self.Tupp = 30
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 12; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 2.; self.fshape_r = 1.3
            self.SxTopQ = 0.05; self.SxBotQ = 0.015
            self.SeedSize = 5; self.PlantPop = 57_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.98; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 18; self.WPy = 60; self.fsink = 0.5
            self.HI0 = 0.35; self.dHI_pre = 5; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 10
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.15; self.p_up2 = 0.6; self.p_up3 = 0.7;  self.p_up4 = 0.85
            self.p_lo1 = 0.65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 2.5; self.fshape_w2 = 2.5; self.fshape_w3 = 2.5; self.fshape_w4 = 1
        elif crop_code == 270:
            self.Name = 'Rapeseed'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 0; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 12; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 1.; self.fshape_r = 1.8
            self.SxTopQ = 0.05; self.SxBotQ = 0.015
            self.SeedSize = 5; self.PlantPop = 440_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.8; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 18.6; self.WPy = 100; self.fsink = 0.5
            self.HI0 = 0.25; self.dHI_pre = 5; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 10
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.2; self.p_up2 = 0.6; self.p_up3 = 0.7;  self.p_up4 = 0.85
            self.p_lo1 = 0.55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.5; self.fshape_w2 = 5; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 289:
            self.Name = 'Sesame'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 35
            
            if tmax >= 35: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 35; self.Tmax_lo = 40
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 10; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 4.5; self.PlantPop = 300_000 if PlantPop == 0 else PlantPop
            self.CCx = .95; self.Kcb = 1.1; self.fage = .3
            self.WP = 18; self.WPy = 100; self.fsink = .5
            self.HI0 = .23; self.dHI_pre = 3; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 10
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .5; self.p_up2 = .8; self.p_up3 = .9;  self.p_up4 = .95
            self.p_lo1 = .8; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 1; self.fshape_w2 = 2; self.fshape_w3 = 3; self.fshape_w4 = 1.
        elif crop_code == 328:
            self.Name = 'Cotton'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 12; self.Tupp = 35
            
            if tmax >= 43: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 43; self.Tmax_lo = 48
                
            if tmin <= 15: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 15; self.Tmin_lo = 10
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 0
            self.GDD_up = 0; self.GDD_lo = 0
            self.Zmin = 0.3; self.Zmax = 2; self.fshape_r =1.5
            self.SxTopQ = 0.052; self.SxBotQ = 0.015
            self.SeedSize = 6; self.PlantPop = 120_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.98; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 15; self.WPy = 70; self.fsink = 0.5
            self.HI0 = 0.35; self.dHI_pre = 5; self.a_HI = 2; self.b_HI = 10; self.dHI0 = 30
            self.Determinant = 0; self.exc = 200
            self.p_up1 = 0.2; self.p_up2 = 0.65; self.p_up3 = 0.75;  self.p_up4 = 0.85
            self.p_lo1 = 0.7; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3; self.fshape_w2 = 2.5; self.fshape_w3 = 2.5; self.fshape_w4 = 1
        elif crop_code == 339:
            self.Name = 'Oilseeds'; self.CropType = 3
            self.PlantMethod = 1; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 0; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 12; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.8
            self.SxTopQ = .05; self.SxBotQ = .015
            self.SeedSize = 5; self.PlantPop = 440_000 if PlantPop == 0 else PlantPop
            self.CCx = .8; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 18; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .25; self.dHI_pre = 5; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 10
            self.Determinant = 1; self.exc = 100
            self.p_up1 = 0.2; self.p_up2 = 0.6; self.p_up3 = 0.7;  self.p_up4 = 0.85
            self.p_lo1 = 0.55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.5; self.fshape_w2 = 5; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 388:
            self.Name = 'Tomato'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 7; self.Tupp = 28
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 0
            self.GDD_up = 0; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = 0.024; self.SxBotQ = 0.006
            self.SeedSize = 20; self.PlantPop = 33_000 if PlantPop == 0 else PlantPop
            self.CCx = 0.75; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 18.; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .63; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = .15; self.p_up2 = .5; self.p_up3 = .7;  self.p_up4 = .92
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 401:
            self.Name = 'Pepper'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 7; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 0
            self.GDD_up = 0; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 15.; self.PlantPop = 65_000 if PlantPop == 0 else PlantPop
            self.CCx = .75; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 20.; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .5; self.dHI_pre = 4; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = .15; self.p_up2 = .45; self.p_up3 = .45;  self.p_up4 = 1
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 403:
            self.Name = 'Onion'; self.CropType = 2
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 0; self.GDD_lo = 0
            self.Zmin = .15; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 5.; self.PlantPop = 330_000 if PlantPop == 0 else PlantPop
            self.CCx = .65; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 19.; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .8; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = .3; self.p_up2 = .5; self.p_up3 = .92;  self.p_up4 = 1
            self.p_lo1 = .65; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 463:
            self.Name = 'Vegetables'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 7; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 0
            self.GDD_up = 0; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .024; self.SxBotQ = .006
            self.SeedSize = 20; self.PlantPop = 50_000 if PlantPop == 0 else PlantPop
            self.CCx = .8; self.Kcb = 1.1; self.fage = .3
            self.WP = 18.; self.WPy = 100; self.fsink = .5
            self.HI0 = .65; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = .15; self.p_up2 = .5; self.p_up3 = .7;  self.p_up4 = .92
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 10000:
            self.Name = 'Vegetables, leafy'; self.CropType = 1
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 7; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = -999; self.GDD_lo = -999
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .024; self.SxBotQ = .006
            self.SeedSize = 7.; self.PlantPop = 80_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = 1.1; self.fage = .3
            self.WP = 15.; self.WPy = 100; self.fsink = .5
            self.HI0 = .65; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .5; self.p_up2 = .5; self.p_up3 = .7;  self.p_up4 = .9
            self.p_lo1 = .8; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 10001:
            self.Name = 'Vegetables, tuber'; self.CropType = 2
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = -999; self.GDD_lo = -999
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 1.; self.PlantPop = 100_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = 1.1; self.fage = .3
            self.WP = 18.; self.WPy = 100; self.fsink = .5
            self.HI0 = .75; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .2; self.p_up2 = .65; self.p_up3 = .75;  self.p_up4 = 1
            self.p_lo1 = .6; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 486:
            self.Name = 'Banana'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 40
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 410_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = 1.1; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 489:
            self.Name = 'Plantains'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 40
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 410_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = 1.1; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 490:
            self.Name = 'Oranges'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 35
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 2.; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 390_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = .8; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 515:
            self.Name = 'Apples'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.5; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 280_000 if PlantPop == 0 else PlantPop
            self.CCx = .8; self.Kcb = .95; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 558:
            self.Name = 'Berries'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .06
            self.SeedSize = 15.; self.PlantPop = 403_333 if PlantPop == 0 else PlantPop
            self.CCx = .7; self.Kcb = .8; self.fage = .3
            self.WP = 15.; self.WPy = 100; self.fsink = .5
            self.HI0 = .5; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .21; self.p_up2 = .5; self.p_up3 = .85;  self.p_up4 = .9
            self.p_lo1 = .51; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 560:
            self.Name = 'Grapes'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 8: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 8; self.Tmin_lo = 3
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 1
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.2; self.Zmax = 1.2; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .06
            self.SeedSize = 15.; self.PlantPop = 403_333 if PlantPop == 0 else PlantPop
            self.CCx = .65; self.Kcb = .85; self.fage = 0.3
            self.WP = 15; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .5; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .21; self.p_up2 = .5; self.p_up3 = .85;  self.p_up4 = .9
            self.p_lo1 = .51; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 572:
            self.Name = 'Avocado'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 35
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.5; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 390_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = .85; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .2; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 603:
            self.Name = 'Tropical fruit'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 35
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 2.; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 410_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = .8; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .1; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 641:
            self.Name = 'Alfalfa'; self.CropType = 1
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = -999; self.GDD_lo = -999
            self.Zmin = 2.; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .02; self.SxBotQ = .01
            self.SeedSize = 40.; self.PlantPop = 2_000_000 if PlantPop == 0 else PlantPop
            self.CCx = .97; self.Kcb = 1.1; self.fage = .3
            self.WP = 15.; self.WPy = 100; self.fsink = .5
            self.HI0 = 1.; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .15; self.p_up2 = .6; self.p_up3 = .98;  self.p_up4 = .9
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 656:
            self.Name = 'Coffee'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 35
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 2.; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 390_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = .8; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .1; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 661:
            self.Name = 'Cocoa'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 10; self.Tupp = 35
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.5; self.Zmax = 1.5; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 390_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = 1.; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .1; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 667:
            self.Name = 'Tea'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 8; self.Tupp = 32
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 2.; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 220_000 if PlantPop == 0 else PlantPop
            self.CCx = .9; self.Kcb = 1.1; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .1; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 821:
            self.Name = 'Fibres'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 5; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 1.7; self.Zmax = 1.7; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 250_000 if PlantPop == 0 else PlantPop
            self.CCx = .5; self.Kcb = .7; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .4; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        elif crop_code == 826:
            self.Name = 'Tobacco'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 7; self.Tupp = 30
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 1; self.PolColdStress = 1; self.TrColdStress = 0
            self.GDD_up = 0; self.GDD_lo = 0
            self.Zmin = .3; self.Zmax = 1.; self.fshape_r = 1.5
            self.SxTopQ = .048; self.SxBotQ = .012
            self.SeedSize = 40.; self.PlantPop = 20_000 if PlantPop == 0 else PlantPop
            self.CCx = .75; self.Kcb = 1.1; self.fage = 0.3
            self.WP = 20.; self.WPy = 100; self.fsink = 0.5
            self.HI0 = .2; self.dHI_pre = 4; self.a_HI = 0; self.b_HI = 3; self.dHI0 = 15
            self.Determinant = 1; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .7;  self.p_up4 = .92
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        
        elif crop_code == 836:
            self.Name = 'Rubber'; self.CropType = 3
            self.PlantMethod = 0; self.CalendarType = 2; self.SwitchGDD = 0
            self.Tbase = 15; self.Tupp = 40
            
            if tmax >= 40: self.Tmax_up = tmax*1; self.Tmax_lo = tmax + 5
            else: self.Tmax_up = 40; self.Tmax_lo = 45
            if tmin <= 10: self.Tmin_up = tmin*1; self.Tmin_lo = tmin - 5
            else: self.Tmin_up = 10; self.Tmin_lo = 5
            
            self.PolHeatStress = 0; self.PolColdStress = 0; self.TrColdStress = 0
            self.GDD_up = 8; self.GDD_lo = 3
            self.Zmin = 2.; self.Zmax = 2.; self.fshape_r = 1.5
            self.SxTopQ = .025; self.SxBotQ = .006
            self.SeedSize = 200.; self.PlantPop = 250_000 if PlantPop == 0 else PlantPop
            self.CCx = .85; self.Kcb = 1.25; self.fage = .1
            self.WP = 17; self.WPy = 100; self.fsink = .5
            self.HI0 = .1; self.dHI_pre = 0; self.a_HI = 0; self.b_HI = 0; self.dHI0 = 0
            self.Determinant = 0; self.exc = 100
            self.p_up1 = .25; self.p_up2 = .5; self.p_up3 = .65;  self.p_up4 = .85
            self.p_lo1 = .55; self.p_lo2 = 1; self.p_lo3 = 1; self.p_lo4 = 1 
            self.fshape_w1 = 3.; self.fshape_w2 = 3; self.fshape_w3 = 3; self.fshape_w4 = 1
        else:
            raise ValueError('Wrong crop name')

        # set any paramaters specified by user
        allowed_keys = {'Ksccx','Ksexpf','Kswp','fcdecline','sfertstress','need_calib','RelativeBio','Ksccx_in',\
                        'fcdecline_in','Ksccx_es','Ksexpf_es','Kswp_es','fcdecline_es','sf_es','relbio_es','HI0','CCx'}#allow to initialize soil fertility stress coefficients

        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allowed_keys)
        
        self.calculate_additional_params()
    
    def calibrate_phenology(self, cell, conf):
        "Calibrate phenology in steps"

        import modules.acea.acea_core as ac, pickle
        cal_step = 15
        # Read inputs
        with open(ac.GetPath(['data','acea','climate', f'{conf.climate_name}_{int(cell["id30"])}.pckl']), 'rb') as f:
            _tmax, _tmin,_,_ = pickle.load(f)
            # Calculate the adjusted Tmin and Tmax based on Tbase and Tupper
            _tmax[_tmax > self.Tupp] = self.Tupp
            _tmax[_tmax < self.Tbase] = self.Tbase
            _tmin[_tmin > self.Tupp] = self.Tupp
            _tavg = np.round((_tmax+_tmin)/2, 1)
            _tavg[_tavg < self.Tbase] = self.Tbase
            
            data = {'GDDcum': np.cumsum(_tavg - self.Tbase),
                    'Date': pd.date_range(start=np.datetime64(str(conf.climate_start)+'-01-01'),\
                                          end=np.datetime64(str(conf.climate_end)+'-12-31'))}
            wdf = pd.DataFrame(data)
       

        cal_step = 15
        cal_steps = [*range(conf.climate_start+cal_step,conf.climate_end+1,cal_step)]
        cal_steps_end = [*range(conf.climate_start+cal_step*2,conf.climate_end+1,cal_step),conf.climate_end]
        
        count = 0; avg_ggds = []
        for mid_year in cal_steps:
            start_year = mid_year-cal_step
            end_year = cal_steps_end[count]; gdds_cum = []
            for cur_year in range(start_year,end_year+1):
                _planting = pd.to_datetime(f'{self.PlantingDate}/{cur_year}'); _harvest = pd.to_datetime(f'{self.HarvestDate}/{cur_year}'); 
                _ggd_start = wdf[wdf['Date'] == _planting]['GDDcum'].values[0]
                
                if _harvest < _planting:
                    _harvest = pd.to_datetime(f'{self.HarvestDate}/{cur_year+1}')
                    if len(wdf[wdf['Date'] == _harvest])==0: continue
                        
                _ggd_end = wdf[wdf['Date'] == _harvest]['GDDcum'].values[0]
                gdds_cum.append(int(_ggd_end-_ggd_start))
                
            # Average growing season GGDs
            gdds_cum = np.array(gdds_cum)
            if any(gdds_cum!=0): avg_GDD = int(np.nanmean(gdds_cum[gdds_cum!=0]))
            else: avg_GDD = 300 # If the environment is very cold, take 300 units
            avg_ggds.append(avg_GDD)

            count +=1
        # Interpolate 
        line = np.poly1d(np.polyfit(cal_steps, avg_ggds, 1))
        
        # Store adjustment coefficients
        self.phenology_calibration = True
        self.coef_years = np.array([*range(conf.climate_start,conf.climate_end+1)], dtype=('int64'))
        if conf.crop_phenology == 'average': self.coefs = np.ones(np.shape(self.coef_years)) * np.mean(line(self.coef_years)/self._Maturity)
        else: self.coefs = line(self.coef_years)/self._Maturity  
               
    def calculate_additional_params(self,): # Copied from classes.py
        # Calculate additional parameters for all self types in mix
        
        # Fractional canopy cover size at emergence
        self.CC0 = self.PlantPop*self.SeedSize*1e-8
        # Root extraction terms
        SxTopQ = self.SxTopQ; SxBotQ = self.SxBotQ
        S1 = self.SxTopQ; S2 = self.SxBotQ
        if S1 == S2:
            SxTop = S1; SxBot = S2
        else:
            if SxTopQ < SxBotQ:
                S1 = SxBotQ; S2 = SxTopQ
            xx = 3*(S2/(S1-S2))
            if xx < 0.5:
                SS1 = (4/3.5)*S1; SS2 = 0
            else:
                SS1 = (xx+3.5)*(S1/(xx+3)); SS2 = (xx-0.5)*(S2/xx)

            if SxTopQ > SxBotQ:
                SxTop = SS1; SxBot = SS2
            else:
                SxTop = SS2; SxBot = SS1

        self.SxTop = SxTop; self.SxBot = SxBot

        # Water stress thresholds
        self.p_up = np.array([self.p_up1,self.p_up2,self.p_up3,self.p_up4])
        self.p_lo = np.array([self.p_lo1,self.p_lo2,self.p_lo3,self.p_lo4])
        self.fshape_w = np.array([self.fshape_w1,self.fshape_w2,self.fshape_w3,self.fshape_w4])