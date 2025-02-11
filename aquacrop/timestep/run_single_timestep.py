import os
import numpy as np
import pandas as pd

from aquacrop.entities.output import Output


from ..entities.totalAvailableWater import TAW
from ..entities.moistureDepletion import Dr
from ..entities.crop import CropStructNT

from ..solution.pre_irrigation import pre_irrigation
from ..solution.irrigation import irrigation
from ..solution.capillary_rise import capillary_rise
from ..solution.germination import germination
from ..solution.growth_stage import growth_stage
from ..solution.canopy_cover import canopy_cover
from ..solution.transpiration import transpiration
from ..solution.groundwater_inflow import groundwater_inflow
from ..solution.harvest_index import harvest_index

if os.getenv("DEVELOPMENT"):
    from ..solution.growing_degree_day import growing_degree_day
    from ..solution.drainage import drainage
    from ..solution.root_zone_water import root_zone_water
    from ..solution.rainfall_partition import rainfall_partition
    from ..solution.check_groundwater_table import check_groundwater_table
    from ..solution.soil_evaporation import soil_evaporation
    from ..solution.root_development import root_development
    from ..solution.infiltration import infiltration
    from ..solution.HIref_current_day import HIref_current_day
    from ..solution.biomass_accumulation import biomass_accumulation
else:
    from ..solution.solution_growing_degree_day import growing_degree_day
    from ..solution.solution_drainage import drainage
    from ..solution.solution_root_zone_water import root_zone_water
    from ..solution.solution_rainfall_partition import rainfall_partition
    from ..solution.solution_check_groundwater_table import check_groundwater_table
    from ..solution.solution_soil_evaporation import soil_evaporation
    from ..solution.solution_root_development import root_development
    from ..solution.solution_infiltration import infiltration
    from ..solution.solution_HIref_current_day import HIref_current_day
    from ..solution.solution_biomass_accumulation import biomass_accumulation

from typing import Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    # Important: classes are only imported when types are checked, not in production.
    from numpy import ndarray
    from aquacrop.entities.clockStruct import ClockStruct
    from aquacrop.entities.initParamVariables import InitialCondition
    from aquacrop.entities.paramStruct import ParamStruct
    from aquacrop.entities.output import Output

def solution_single_time_step(
    init_cond: "InitialCondition",
    param_struct: "ParamStruct",
    clock_struct: "ClockStruct",
    weather_step: "ndarray",
    outputs:"Output",
) ->  Tuple["InitialCondition", "ParamStruct","Output"]:
    """
    Function to perform AquaCrop solution for a single time step

    Arguments:

        init_cond (InitialCondition):  containing current variables+counters

        param_struct (ParamStruct):  contains model paramaters

        clock_struct (ClockStruct):  model time paramaters

        weather_step (numpy.ndarray):  contains precipitation,ET,temp_max,temp_min for current day

        outputs (Output):  object to store outputs

    Returns:

        NewCond (InitialCondition):  containing updated simulation variables+counters

        param_struct (ParamStruct):  contains model paramaters

        outputs (Output):  object to store outputs

    """

    # Unpack structures
    Soil = param_struct.Soil
    CO2 = param_struct.CO2
    precipitation = weather_step[2]
    temp_max = weather_step[1]
    temp_min = weather_step[0]
    et0 = weather_step[3]

    # Store initial conditions in structure for updating %%
    NewCond = init_cond
    if param_struct.water_table == 1:
        GroundWater = param_struct.z_gw[clock_struct.time_step_counter]
    else:
        GroundWater = 0

    # Check if growing season is active on current time step %%
    if clock_struct.season_counter >= 0:
        # Check if in growing season
        CurrentDate = clock_struct.step_start_time
        planting_date = clock_struct.planting_dates[clock_struct.season_counter]
        harvest_date = clock_struct.harvest_dates[clock_struct.season_counter]
        harvest_date = harvest_date + ((harvest_date-planting_date) * param_struct.Seasonal_Crop_List[clock_struct.season_counter].crop_gs_increase/100).round('D')

        growing_season = False
        if param_struct.Seasonal_Crop_List[clock_struct.season_counter].crop_perennial == False: # Annual crops
            if (planting_date <= CurrentDate) and (harvest_date >= CurrentDate) and \
               (NewCond.crop_mature == False) and (NewCond.crop_dead == False):
                growing_season = True
        elif planting_date <= CurrentDate: # Perennial crops
            if clock_struct.season_counter+1 == clock_struct.n_seasons:
                if CurrentDate < (clock_struct.simulation_end_date - pd.Timedelta(1, unit="d")):
                    growing_season = True
            elif CurrentDate < (clock_struct.planting_dates[clock_struct.season_counter+1] - pd.Timedelta(1, unit="d")):
                growing_season = True

        # Assign crop, irrigation management, and field management structures
        Crop_ = param_struct.Seasonal_Crop_List[clock_struct.season_counter]
        Crop_Name = param_struct.CropChoices[clock_struct.season_counter]
        IrrMngt = param_struct.IrrMngt

        if growing_season is True:
            FieldMngt = param_struct.FieldMngt
        else:
            FieldMngt = param_struct.FallowFieldMngt

    else:
        # Not yet reached start of first growing season
        growing_season = False
        # Assign crop, irrigation management, and field management structures
        # Assign first crop as filler crop
        Crop_ = param_struct.Fallow_Crop
        Crop_Name = "fallow"

        Crop_.Aer = 5
        Crop_.Zmin = 0.3
        IrrMngt = param_struct.FallowIrrMngt
        FieldMngt = param_struct.FallowFieldMngt

    # Increment time counters %%
    if growing_season is True:
        # Calendar days after planting
        NewCond.dap = NewCond.dap + 1
        # Growing degree days after planting

        gdd = growing_degree_day(
            Crop_.GDDmethod, Crop_.Tupp, Crop_.Tbase, temp_max, temp_min
        )

        # Update cumulative gdd counter
        NewCond.gdd = gdd
        NewCond.gdd_cum = NewCond.gdd_cum + gdd

        NewCond.growing_season = True
        
        if NewCond.Wr0 == 0: 
            NewCond.Wr0,_,_,_,_,_,_,_,_,_,_ = root_zone_water(Soil.Profile,Soil.Profile.dzsum[-1],NewCond.th,Soil.z_top,float(Crop_.Zmin),Crop_.Aer)
        
        
    else:
        NewCond.growing_season = False

        # Calendar days after planting
        NewCond.dap = 0
        # Growing degree days after planting
        gdd = 0.3
        NewCond.gdd_cum = 0

    # save current timestep counter
    NewCond.time_step_counter = clock_struct.time_step_counter
    NewCond.precipitation = weather_step[2]
    NewCond.temp_max = weather_step[1]
    NewCond.temp_min = weather_step[0]
    NewCond.et0 = weather_step[3]


    class_args = {
        key: value
        for key, value in Crop_.__dict__.items()
        if not key.startswith("__") and not callable(key)
    }
    Crop = CropStructNT(**class_args)
    
    # Run simulations %%
    # 1. Check for groundwater table
    NewCond.th_fc_Adj, NewCond.wt_in_soil, NewCond.S_cr = check_groundwater_table(
        Soil.Profile,
        NewCond.z_gw,
        NewCond.th,
        NewCond.th_fc_Adj,
        NewCond.S_cr,
        param_struct.water_table,
        GroundWater
    )

    # 2. Root development
    NewCond.z_root = root_development(
        Crop,
        Soil.Profile,
        NewCond.dap,
        NewCond.z_root,
        NewCond.delayed_cds,
        NewCond.gdd_cum,
        NewCond.delayed_gdds,
        NewCond.tr_ratio,
        NewCond.th,
        NewCond.canopy_cover,
        NewCond.canopy_cover_ns,
        NewCond.germination,
        NewCond.r_cor,
        NewCond.t_pot,
        NewCond.z_gw,
        gdd,
        growing_season,
        param_struct.water_table,
    )

    # 3. Pre-irrigation
    NewCond, PreIrr = pre_irrigation(
        Soil.Profile, Crop, NewCond, growing_season, IrrMngt
    )

    # 4. Drainage

    NewCond.th, DeepPerc, FluxOut, NewCond.S_rain,NewCond.S_irr,NewCond.S_cr = drainage(
        Soil.Profile,
        NewCond.th,
        NewCond.th_fc_Adj,
        NewCond.S_rain,
        NewCond.S_irr,
        NewCond.S_cr,
    )

    # 5. Surface runoff
    Runoff, Infl, NewCond.day_submerged = rainfall_partition(
        precipitation,
        NewCond.th,
        NewCond.day_submerged,
        FieldMngt.sr_inhb,
        FieldMngt.bunds,
        FieldMngt.z_bund,
        FieldMngt.curve_number_adj_pct,
        Soil.cn,
        Soil.adj_cn,
        Soil.z_cn,
        Soil.nComp,
        Soil.Profile,
    )

    # 6. Irrigation
    NewCond.depletion, NewCond.taw, NewCond.irr_cum, Irr = irrigation(
        IrrMngt.irrigation_method,
        IrrMngt.SMT,
        IrrMngt.AppEff,
        IrrMngt.MaxIrr,
        IrrMngt.IrrInterval,
        IrrMngt.Schedule,
        IrrMngt.depth,
        IrrMngt.MaxIrrSeason,
        NewCond.growth_stage,
        NewCond.irr_cum,
        NewCond.e_pot,
        NewCond.t_pot,
        NewCond.z_root,
        NewCond.th,
        NewCond.dap,
        NewCond.time_step_counter,
        Crop,
        Soil.Profile,
        Soil.z_top,
        growing_season,
        precipitation,
        Runoff,
    )

    # 7. Infiltration
    (
        NewCond.th,
        NewCond.surface_storage,
        DeepPerc,
        Runoff,
        Infl,
        FluxOut,
        NewCond.S_rain, 
        NewCond.S_irr,
    ) = infiltration(
        Soil.Profile,
        NewCond.surface_storage,
        NewCond.th_fc_Adj,
        NewCond.th,
        Infl,
        Irr,
        IrrMngt.AppEff,
        FieldMngt.bunds,
        FieldMngt.z_bund,
        FluxOut,
        DeepPerc,
        Runoff,
        growing_season,
        NewCond.S_rain, 
        NewCond.S_irr,
    )
    # 8. Capillary Rise
    NewCond, CR = capillary_rise(
        Soil.Profile,
        Soil.nLayer,
        Soil.fshape_cr,
        NewCond,
        FluxOut,
        param_struct.water_table,
    )

    # 9. Check germination
    NewCond = germination(
        NewCond,
        Soil.z_germ,
        Soil.Profile,
        Crop.GermThr,
        Crop.PlantMethod,
        gdd,
        growing_season,
    )

    # 10. Update growth stage
    NewCond = growth_stage(Crop, NewCond, growing_season)

    # 11. Canopy cover development
    NewCond = canopy_cover(
        Crop, Soil.Profile, Soil.z_top, NewCond, gdd, et0, growing_season
    )

    # 12. Soil evaporation
    (
        NewCond.e_pot,
        NewCond.th,
        NewCond.stage2,
        NewCond.w_stage_2,
        NewCond.w_surf,
        NewCond.surface_storage,
        NewCond.evap_z,
        Es,
        EsPot,
        NewCond.ET_rain,
        NewCond.ET_irr,
        NewCond.ET_cr,
        NewCond.S_rain,
        NewCond.S_irr,
        NewCond.S_cr,
    ) = soil_evaporation(
        clock_struct.evap_time_steps,
        clock_struct.sim_off_season,
        clock_struct.time_step_counter,
        Soil.Profile,
        Soil.evap_z_min,
        Soil.evap_z_max,
        Soil.rew,
        Soil.kex,
        Soil.fwcc,
        Soil.f_wrel_exp,
        Soil.f_evap,
        Crop.CalendarType,
        Crop.Senescence,
        IrrMngt.irrigation_method,
        IrrMngt.WetSurf,
        FieldMngt.mulches,
        FieldMngt.f_mulch,
        FieldMngt.mulch_pct,
        NewCond.dap,
        NewCond.w_surf,
        NewCond.evap_z,
        NewCond.stage2,
        NewCond.th,
        NewCond.delayed_cds,
        NewCond.gdd_cum,
        NewCond.delayed_gdds,
        NewCond.ccx_w,
        NewCond.canopy_cover_adj,
        NewCond.ccx_act,
        NewCond.canopy_cover,
        NewCond.premat_senes,
        NewCond.surface_storage,
        NewCond.w_stage_2,
        NewCond.e_pot,
        et0,
        Infl,
        precipitation,
        Irr,
        growing_season,
        NewCond.Wetted_area,
        NewCond.ET_rain,
        NewCond.ET_irr,
        NewCond.ET_cr,
        NewCond.S_rain,
        NewCond.S_irr,
        NewCond.S_cr,
    )

    # 13. Crop transpiration
    Tr, TrPot_NS, TrPot, NewCond, IrrNet = transpiration(
        Soil.Profile,
        Soil.nComp,
        Soil.z_top,
        Crop,
        IrrMngt.irrigation_method,
        IrrMngt.NetIrrSMT,
        NewCond,
        et0,
        CO2,
        growing_season,
        gdd,
    )

    # 14. Groundwater inflow
    NewCond, GwIn = groundwater_inflow(Soil.Profile, NewCond)

    # 15. Reference harvest index
    (NewCond.hi_ref, NewCond.yield_form, NewCond.pct_lag_phase) = HIref_current_day( # ,NewCond.HIfinal
        NewCond.hi_ref,
        NewCond.HIfinal,
        NewCond.dap,
        NewCond.delayed_cds,
        NewCond.yield_form,
        NewCond.pct_lag_phase,
        #NewCond.cc_prev,
        NewCond.canopy_cover,
        NewCond.ccx_w,
        Crop,
        growing_season,
    )

    # 16. Biomass accumulation
    (NewCond.biomass, NewCond.biomass_ns,NewCond.StressSFadjNEW,NewCond.StressSFadjpre,NewCond.Tr_ET0_accum,NewCond.WPadj,Crop,) = biomass_accumulation(
        Crop,
        NewCond.dap,
        NewCond.delayed_cds,
        NewCond.hi_ref,
        NewCond.pct_lag_phase,
        NewCond.biomass,
        NewCond.biomass_ns,
        Tr,
        TrPot_NS,
        et0,
        growing_season,
        NewCond.StressSFadjNEW,
        NewCond.StressSFadjpre,
        NewCond.Tr_ET0_accum,
        NewCond.WPadj,
    )
    
    # Update global variables
    param_struct.Seasonal_Crop_List[clock_struct.season_counter].Ksccx=Crop.Ksccx
    param_struct.Seasonal_Crop_List[clock_struct.season_counter].Ksexpf=Crop.Ksexpf
    param_struct.Seasonal_Crop_List[clock_struct.season_counter].Kswp=Crop.Kswp
    param_struct.Seasonal_Crop_List[clock_struct.season_counter].fcdecline=Crop.fcdecline
    param_struct.Seasonal_Crop_List[clock_struct.season_counter].MaxCanopyCD=Crop.MaxCanopyCD

    # 17. Harvest index
    NewCond = harvest_index(
        Soil.Profile, Soil.z_top, Crop, NewCond, et0, temp_max, temp_min, growing_season
    )

    # 18. Yield potential
    NewCond.YieldPot = (NewCond.biomass_ns / 100) * NewCond.harvest_index

    # 19. Crop yield_ (dry and fresh)
    if growing_season is True:
        # Calculate crop yield_ (tonne/ha)
        NewCond.DryYield = round((NewCond.biomass / 100) * NewCond.harvest_index_adj,3)
        NewCond.FreshYield = NewCond.DryYield / (1 - (Crop.YldWC / 100))
        # print( clock_struct.time_step_counter,(NewCond.biomass/100),NewCond.harvest_index_adj)
        # Check if crop has reached maturity
        if ((Crop.CalendarType == 1) and (NewCond.dap >= Crop.Maturity)) or (
            (Crop.CalendarType == 2) and (NewCond.gdd_cum >= Crop.Maturity)
        ):
            # Crop has reached maturity
            NewCond.crop_mature = True

    elif growing_season is False:
        # Crop yield_ is zero outside of growing season
        NewCond.DryYield = 0
        NewCond.FreshYield = 0

    # 20. Root zone water
    _TAW = TAW()
    _water_root_depletion = Dr()
    # thRZ = RootZoneWater()

    Wr, _water_root_depletion.Zt, _water_root_depletion.Rz, _TAW.Zt, _TAW.Rz, _, _, _, _, _, _ = root_zone_water(
        Soil.Profile,
        float(NewCond.z_root),
        NewCond.th,
        Soil.z_top,
        float(Crop.Zmin),
        Crop.Aer,
    )

    # 21. Update net irrigation to add any pre irrigation
    IrrNet = IrrNet + PreIrr
    NewCond.irr_net_cum = NewCond.irr_net_cum + PreIrr

    # Update model outputs %%
    row_day = clock_struct.time_step_counter
    row_gs = clock_struct.season_counter

    # Irrigation
    if growing_season is True:

        NewCond.Pcum = NewCond.Pcum + precipitation
        NewCond.Ecum = NewCond.Ecum + Es
        NewCond.Tcum = NewCond.Tcum + Tr
        NewCond.Rcum = NewCond.Rcum + Runoff
        NewCond.Percum = NewCond.Percum + DeepPerc
        NewCond.CRcum = NewCond.CRcum + CR
        NewCond.GWcum = NewCond.GWcum + GwIn
    
        if IrrMngt.irrigation_method == 4:
            # Net irrigation
            IrrDay = IrrNet
            IrrTot = NewCond.irr_net_cum
        else:
            # Irrigation
            IrrDay = Irr
            IrrTot = NewCond.irr_cum

    else:
        IrrDay = 0
        IrrTot = 0

        NewCond.depletion = _water_root_depletion.Rz
        NewCond.taw = _TAW.Rz

    # Water contents
    outputs.water_storage[row_day, :3] = np.array(
        [clock_struct.time_step_counter, growing_season, NewCond.dap]
    )
    outputs.water_storage[row_day, 3:] = NewCond.th

    # Water fluxes
    # print(f'Saving NewCond.z_gw to outputs: {NewCond.z_gw}')
    outputs.water_flux[row_day, :] = [
        clock_struct.time_step_counter,
        clock_struct.season_counter,
        NewCond.dap,
        Wr,
        NewCond.z_gw,
        NewCond.surface_storage,
        IrrDay,
        Infl,
        Runoff,
        DeepPerc,
        CR,
        GwIn,
        Es,
        EsPot,
        Tr,
        TrPot,
    ]

    # Crop growth
    outputs.crop_growth[row_day, :] = [
        clock_struct.time_step_counter,
        clock_struct.season_counter,
        NewCond.dap,
        gdd,
        NewCond.gdd_cum,
        NewCond.z_root,
        NewCond.canopy_cover,
        NewCond.canopy_cover_ns,
        NewCond.biomass,
        NewCond.biomass_ns,
        NewCond.harvest_index,
        NewCond.harvest_index_adj,
        NewCond.DryYield,
        NewCond.FreshYield,
        NewCond.YieldPot,
        Tr,
        TrPot_NS,
        TrPot,
        Tr/et0,
        NewCond.WPadj
    ]
    
    outputs.ET_color[row_day, :] = [
        clock_struct.time_step_counter,
        clock_struct.season_counter,
        NewCond.dap,
        sum(NewCond.ET_rain),
        sum(NewCond.ET_irr),
        sum(NewCond.ET_cr),
    ]
    
    outputs.S_color[row_day, :] = [
        clock_struct.time_step_counter,
        clock_struct.season_counter,
        NewCond.dap,
        sum(NewCond.S_rain),
        sum(NewCond.S_irr),
        sum(NewCond.S_cr),
    ]

    # Final output (if at end of growing season)
    
    hrv_date = clock_struct.harvest_dates[clock_struct.season_counter]
    plnt_date = clock_struct.planting_dates[clock_struct.season_counter]
    hrv_date_max = hrv_date + ((hrv_date-plnt_date) * Crop.crop_gs_increase/100).round('D')
    if clock_struct.season_counter+1 == clock_struct.n_seasons: 
        hrv_date_crit = clock_struct.simulation_end_date-pd.Timedelta(1, unit="d")
    else: 
        hrv_date_crit = clock_struct.planting_dates[clock_struct.season_counter+1]-pd.Timedelta(14, unit="d")
    if hrv_date_max >= hrv_date_crit: 
        hrv_date_max = hrv_date_crit
    
    
    if clock_struct.season_counter > -1:
        if ((NewCond.crop_mature == True) or (NewCond.crop_dead == True) or (clock_struct.step_end_time >= hrv_date_max)) \
        and (NewCond.harvest_flag == False):
            
            if NewCond.Wr_end == 0: 
                NewCond.Wr_end,_,_,_,_,_,_,_,_,_,_ = root_zone_water(Soil.Profile,Soil.Profile.dzsum[-1],NewCond.th,Soil.z_top,float(Crop.Zmin),Crop.Aer)
            
            # Store final outputs
            anthesis_date = clock_struct.planting_dates[clock_struct.season_counter] + pd.Timedelta(Crop.HIstartCD, unit="d")
            if not clock_struct.season_counter in outputs.final_stats.index:
                # Store final outputs
                outputs.final_stats.loc[row_gs] = [
                    clock_struct.season_counter,
                    Crop_Name,
                    plnt_date,
                    anthesis_date,
                    clock_struct.step_end_time,
                    clock_struct.time_step_counter,
                    NewCond.DryYield,
                    NewCond.FreshYield,
                    NewCond.YieldPot,
                    IrrTot,
                    NewCond.gdd_cum,
                    NewCond.Wr0, 
                    NewCond.Pcum, 
                    NewCond.CRcum, 
                    NewCond.GWcum,
                    NewCond.Wr_end, 
                    NewCond.Ecum, 
                    NewCond.Tcum, 
                    NewCond.Rcum, 
                    NewCond.Percum,
                    NewCond.biomass/100,
                ]           

            if Crop.crop_perennial == False: # Annual crops
            
                outputs.final_watercolor.loc[clock_struct.season_counter] = [clock_struct.season_counter,
                    sum(NewCond.ET_rain),
                    sum(NewCond.ET_irr), 
                    sum(NewCond.ET_cr),
                    sum(NewCond.S_rain), 
                    sum(NewCond.S_irr), 
                    sum(NewCond.S_cr),]

                NewCond.harvest_flag = True # Set harvest flag
            else: # Perennial crops
                temp_flag_=False
                if clock_struct.season_counter+1 == clock_struct.n_seasons:
                    if clock_struct.step_end_time == (clock_struct.simulation_end_date-pd.Timedelta(1, unit="d")):
                        temp_flag_=True

                elif clock_struct.step_end_time == (clock_struct.planting_dates[clock_struct.season_counter+1] - pd.Timedelta(1, unit="d")):
                    temp_flag_=True
                    
                if temp_flag_:
                    outputs.final_watercolor.loc[clock_struct.season_counter] = [clock_struct.season_counter,
                            sum(NewCond.ET_rain),
                            sum(NewCond.ET_irr), 
                            sum(NewCond.ET_cr),
                            sum(NewCond.S_rain), 
                            sum(NewCond.S_irr), 
                            sum(NewCond.S_cr),]
                    NewCond.harvest_flag = True # Set harvest flag
                    outputs.final_stats.iloc[clock_struct.season_counter,11:] = [
                            NewCond.Wr0, 
                            NewCond.Pcum,  
                            NewCond.CRcum, 
                            NewCond.GWcum,
                            NewCond.Wr_end, 
                            NewCond.Ecum, 
                            NewCond.Tcum, 
                            NewCond.Rcum, 
                            NewCond.Percum,
                            NewCond.biomass/100,
                        ]

    return NewCond, param_struct, outputs
