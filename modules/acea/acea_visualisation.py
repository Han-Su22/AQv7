# -*- coding: utf-8 -*-
"""
Created on Fri May 28 10:00:26 2021

@author: MialykO
"""
#%% Import packages
import numpy as np
from matplotlib.cm import get_cmap
import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mtick
import os
import modules.acea.acea_core as ac

#%% Classes
class acea_visualiser:
    
    def __init__(self, conf):
        self.conf = conf
        self.data_folder = ac.GetPath(["data","acea"])
        self.raster_folder_path = ac.GetPath(["outputs", f"{self.conf.project_name}",'rasters'])
        self.clock_start_year = int(self.conf.clock_start.split('/')[0])
        self.clock_end_year = int(self.conf.clock_end.split('/')[0])
        self.max_gs = (self.clock_end_year - self.clock_start_year)-1; self.first_year = self.clock_end_year-self.max_gs+1
        self.lats, self.lons, self.res5 = ac.GetResolution(1)
        self.text_years = f'{self.first_year}_{self.clock_end_year}'
        
        # Get main ACEA outputs
        self.GetMainOutputs()
        
    def GetMainOutputs(self):
        
        # Get production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_rf_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        self.prod_rf = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Rainfed production
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_prod_ir_global_annual_{self.first_year}_{self.clock_end_year}.nc'])
        self.prod_ir = np.ma.filled(ac.ReadNC(_raster_path, f'production-{self.conf.crop_name_short}'),0) # Irrigated production
        
        self.prod_sum = self.prod_rf + self.prod_ir; self.prod_sum = np.ma.masked_where(self.prod_sum == 0, self.prod_sum) # Total production
        prod_rf_frac = self.prod_rf/self.prod_sum; prod_ir_frac = self.prod_ir/self.prod_sum
        
        # Get WFs
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_rf_global_annual_{self.text_years}.nc'])
        wf_green_rf = ac.ReadNC(_raster_path, f'wf_green_rf-{self.conf.crop_name_short}')
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_rf_global_annual_{self.text_years}.nc'])
        wf_blue_rf = ac.ReadNC(_raster_path, f'wf_blue_rf-{self.conf.crop_name_short}')
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_green_ir_global_annual_{self.text_years}.nc'])
        wf_green_ir = np.ma.filled(ac.ReadNC(_raster_path, f'wf_green_ir-{self.conf.crop_name_short}'),0)
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_ir_global_annual_{self.text_years}.nc'])
        wf_blue_ir = np.ma.filled(ac.ReadNC(_raster_path, f'wf_blue_ir-{self.conf.crop_name_short}'),0)
        
        self.wf_gw = wf_blue_rf + wf_green_rf # Only the cells with shallow GW (blue masks green values)
        wf_green_rf = np.ma.filled(wf_green_rf,0); wf_blue_rf = np.ma.filled(wf_blue_rf,0)
        
        # Get total rainfed and irr WF
        self.wf_rf = wf_green_rf + wf_blue_rf; self.wf_ir = wf_green_ir + wf_blue_ir
        self.wf_rf = np.ma.array(self.wf_rf, mask=(self.wf_rf==0)); self.wf_ir = np.ma.array(self.wf_ir, mask=(self.wf_ir==0))
        
        # Get total WFs
        wf_green_rf = wf_green_rf * prod_rf_frac; wf_blue_rf = wf_blue_rf * prod_rf_frac
        wf_green_ir = wf_green_ir * prod_ir_frac; wf_blue_ir = wf_blue_ir * prod_ir_frac
        
        self.wf_green = wf_green_rf + wf_green_ir # Add up rainfed and irr
        self.wf_blue = wf_blue_rf + wf_blue_ir # Add up rainfed and irr
        self.wf_tot = self.wf_green + self.wf_blue
        self.wf_tot = np.ma.array(self.wf_tot, mask=(self.wf_tot==0))
        self.wf_green = np.ma.array(self.wf_green, mask=(self.wf_green==0)); self.wf_blue = np.ma.array(self.wf_blue, mask=(self.wf_blue==0))
    
    def PlotWFMaps(self, save_raster = True):
        "Plot WF maps of different periods"
        
        # Unit WF
        wf_rf5_last = np.ma.mean(self.wf_rf[-5:,:,:], axis=0) # 2012-2016
        wf_ir5_last = np.ma.mean(self.wf_ir[-5:,:,:], axis=0) # 2012-2016
        wf_green5_last = np.ma.mean(self.wf_green[-5:,:,:], axis=0) # 2012-2016
        wf_blue5_last = np.ma.mean(self.wf_blue[-5:,:,:], axis=0) # 2012-2016
        
        wf_tot5_first = np.ma.mean(self.wf_tot[:5,:,:], axis=0) # 1986-1990
        wf_tot5_last = np.ma.mean(self.wf_tot[-5:,:,:], axis=0) # 2012-2016
        wf_tot5_change = (wf_tot5_last/wf_tot5_first - 1) * 100
        
        # WF of production
        wc_tot5_first = np.ma.mean(self.wf_tot[:5,:,:]*self.prod_sum[:5,:,:], axis=0) # 1986-1990
        wc_tot5_last = np.ma.mean(self.wf_tot[-5:,:,:]*self.prod_sum[-5:,:,:], axis=0) # 2012-2016
        wc_tot5_change = (wc_tot5_last/wc_tot5_first - 1) * 100
        
        self.PlotMapCustomLevels(wf_tot5_last, " ", np.around(self.CalculateBars(wf_tot5_last)/50)*50, False, 'Unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)','RdYlBu_r',True, 200)
        self.PlotMapCustomLevels(wf_green5_last, " ", np.around(self.CalculateBars(wf_green5_last)/50)*50, False, 'Green unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)','RdYlBu_r',True, 200)
        self.PlotMapCustomLevels(wf_blue5_last, " ", np.around(self.CalculateBars(wf_blue5_last)/10)*10, False, 'Blue unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)','RdYlBu_r',True, 200)
        
        # to remove
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_blue_ir_global_annual_{self.text_years}.nc'])
        wf_blue_ir = ac.ReadNC(_raster_path, f'wf_blue_ir-{self.conf.crop_name_short}')
        wf_blue_ir_last = np.ma.mean(wf_blue_ir[-5:,:,:], axis=0) # 2012-2016
        self.PlotMapCustomLevels(wf_blue_ir_last, " ", np.around(self.CalculateBars(wf_blue_ir_last)/10)*10, False, 'Blue unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)','RdYlBu_r',True, 200)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_sc3_cot_cwu_blue_ir_global_annual_{self.text_years}.nc'])
        wf_blue_ir = ac.ReadNC(_raster_path, f'cwu_blue_ir-{self.conf.crop_name_short}')
        wf_blue_ir_last = np.ma.mean(wf_blue_ir[-5:,:,:], axis=0) # 2012-2016
        self.PlotMapCustomLevels(wf_blue_ir_last, " ", np.around(self.CalculateBars(wf_blue_ir_last)/10)*10, False, 'Blue irrigated CWU ($mm$ $y^{-1}$)','RdYlBu_r',True, 200)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_sc3_cot_yield_global_annual_{self.text_years}.nc'])
        wf_blue_ir = ac.ReadNC(_raster_path, f'yield-{self.conf.crop_name_short}')
        wf_blue_ir_last = np.ma.mean(wf_blue_ir[-5:,:,:], axis=0) # 2012-2016
        self.PlotMapCustomLevels(wf_blue_ir_last, " ", np.around(self.CalculateBars(wf_blue_ir_last*10))/10, False, 'Irrigated yield ($t$ $ha^{-1}$ $y^{-1}$)','RdYlBu_r',True, 200)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_sc3_cot_matyday_global_annual_{self.text_years}.nc'])
        wf_blue_ir = ac.ReadNC(_raster_path, f'matyday-{self.conf.crop_name_short}')
        wf_blue_ir_last = np.ma.mean(wf_blue_ir[-5:,:,:], axis=0) # 2012-2016
        self.PlotMapCustomLevels(wf_blue_ir_last, " ", np.around(self.CalculateBars(wf_blue_ir_last*10))/10, False, 'Growing season($d$ $y^{-1}$)','RdYlBu_r',True, 200)
        
        # ------
        self.PlotMapCustomLevels(wc_tot5_last/10**3, " ", np.around(self.CalculateBars(wf_tot5_last)/50)*50, False, 'Water footprint of production ($10^3$ $m^3$ $y^{-1}$)','RdYlBu_r',True, 200)
        
        
        # Rainfed and irrigated unit WF
        self.PlotMapWithDist(wf_rf5_last, " ", self.CalculateBars(wf_rf5_last), 4000, 'Rainfed unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)','RdYlBu_r',0, True, True)
        self.PlotMapWithDist(wf_ir5_last, " ", self.CalculateBars(wf_ir5_last), 2000, 'Irrigated unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)','RdYlBu_r',0, True, True)
        
        # Total 2012-2016
        self.PlotMapWithDist(wf_tot5_last, " ", self.CalculateBars(wf_tot5_last, True, [30, 40, 50, 65, 95]), 4000, 'Unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)','RdYlBu_r',0, True, True)
        self.PlotMapWithDist(wc_tot5_last/10**3, "Average water footprint of production 2012-2016", self.CalculateBars(wc_tot5_last/10**3, True, [40, 55, 65, 85, 95]), \
                                 2000, 'Water footprint of production ($10^3$ $m^3$ $y^{-1}$)','RdYlBu_r', 0, True, True)
        
        # Change from 2012-2016 to 1986-1990
        self.PlotMapCustomLevels(wf_tot5_change, "Change in unit water footprint", np.arange(-80, 81, 10), True, 'Change in unit water footprint', 'RdBu_r')
        self.PlotMapCustomLevels(wc_tot5_change, "Change in water footprint of production", self.CalculateBars(wc_tot5_change, False, [40, 50, 65, 80, 95]), True, \
                                 'Change in water footprint of production', 'RdYlBu_r', True)

        # Save rasters
        if save_raster:
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_tot_change.nc'])
            val_name = f'wf_total-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}'
            ac.CreateRaster5(wf_tot5_change, _raster_path, val_name, 'm3 t-1 y-1', title)
            
            _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wc_tot_change.nc'])
            val_name = f'wf_total-{self.conf.crop_name_short}'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}'
            ac.CreateRaster5(wc_tot5_change, _raster_path, val_name, 'm3 y-1', title)
        
    def PlotOtherMaps(self):
        "Plot ET, CWU, and Y maps of different periods"
        
        # Get harvest areas
        harv_area_rf = np.ma.filled(ac.GetHistHarvestedAreas5(self.conf.crop_fao, False),0); harv_area_ir = np.ma.filled(ac.GetHistHarvestedAreas5(self.conf.crop_fao, True),0)
        harv_area_tot = harv_area_ir + harv_area_rf; harv_area_tot = np.ma.masked_where(harv_area_tot == 0, harv_area_tot)
        harv_area_rf_frac = harv_area_rf/harv_area_tot; harv_area_ir_frac = harv_area_ir/harv_area_tot
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_sc1_{self.conf.crop_name_short}_matyday_global_annual_{self.text_years}.nc'])
        growing_season_rf = np.ma.filled(ac.ReadNC(_raster_path, f'matyday-{self.conf.crop_name_short}'),0)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_sc2_{self.conf.crop_name_short}_matyday_global_annual_{self.text_years}.nc'])
        growing_season_rf = growing_season_rf + np.ma.filled(ac.ReadNC(_raster_path, f'matyday-{self.conf.crop_name_short}'),0)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_sc3_{self.conf.crop_name_short}_matyday_global_annual_{self.text_years}.nc'])
        growing_season_ir = np.ma.filled(ac.ReadNC(_raster_path, f'matyday-{self.conf.crop_name_short}'),0)
        
        growing_season = growing_season_rf*harv_area_rf_frac  + growing_season_ir*harv_area_ir_frac
        
        yields = self.prod_sum/harv_area_tot
        yields5_last = np.ma.mean(yields[-5:,:,:], axis=0)
        
        cwu = self.wf_tot*yields/10
        cwu5_last = np.ma.mean(cwu[-5:,:,:], axis=0)
        
        et = cwu/growing_season
        et5_last = np.ma.mean(et[-5:,:,:], axis=0)
        
        self.PlotMapWithDist(et5_last, "Average maize ET 2012-2016", self.CalculateBars(et5_last*20)/20, 5, 'Evapotranspiration ($mm$ $day^{-1}$)','RdYlBu')
        self.PlotMapWithDist(yields5_last, "Average maize yields 2012-2016", self.CalculateBars(yields5_last*10, True, [40, 65, 75, 85, 99])/10, 12, 'Crop yield ($t$ $ha^{-1}$ $y^{-1}$)','RdYlBu')
        self.PlotMapWithDist(cwu5_last, "Average maize CWU 2012-2016", self.CalculateBars(cwu5_last), 700, 'Crop water use ($mm$ $y^{-1}$)','RdYlBu')
        
    def PlotCVMap(self):
        "Plot the coefficient of variation"
        
        # Get & detrend data
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_tot_trended_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_tot_trended = ac.ReadNC(_raster_path, f'wf_total-{self.conf.crop_name_short}')
        else: wf_tot_trended = self.TrendRasterData(self.wf_tot, 1, _raster_path)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_rf_trended_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_rf_trended = ac.ReadNC(_raster_path, f'wf_total-{self.conf.crop_name_short}')
        else: wf_rf_trended = self.TrendRasterData(self.wf_rf, 1, _raster_path)
       
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_wf_ir_trended_global_annual_{self.text_years}.nc'])
        if os.path.exists(_raster_path): wf_ir_trended = ac.ReadNC(_raster_path, f'wf_total-{self.conf.crop_name_short}')
        else: wf_ir_trended = self.TrendRasterData(self.wf_ir, 1, _raster_path)
        
        wf_tot_detrended = self.wf_tot[0,:,:] + (self.wf_tot - wf_tot_trended)
        wf_rf_detrended = self.wf_rf[0,:,:] + (self.wf_rf - wf_rf_trended) 
        wf_ir_detrended = self.wf_ir[0,:,:] + (self.wf_ir - wf_ir_trended)
        wf_tot_detrended[wf_tot_detrended<0] = 0; wf_rf_detrended[wf_rf_detrended<0] = 0; wf_ir_detrended[wf_ir_detrended<0] = 0
        
        # Analyse variability
        std_wf_tot = np.ma.std(wf_tot_detrended, axis=0); mean_wf_tot = np.ma.mean(wf_tot_detrended, axis = 0); cv_wf_tot = std_wf_tot/mean_wf_tot*100
        std_wf_rf = np.ma.std(wf_rf_detrended, axis=0); mean_wf_rf = np.ma.mean(wf_rf_detrended, axis = 0); cv_wf_rf = std_wf_rf/mean_wf_rf*100
        std_wf_ir = np.ma.std(wf_ir_detrended, axis=0); mean_wf_ir = np.ma.mean(wf_ir_detrended, axis = 0); cv_wf_ir = std_wf_ir/mean_wf_ir*100
        
        self.PlotMapWithDist(cv_wf_tot, "Coeffient of variation of detrended maize water footprints 1986-2016", self.CalculateBars(cv_wf_tot), 40, 'Coefficient of variation','RdYlBu_r', 1, True, True)
        # PlotMapWithDist(cv_wf_ir, self.lats, self.lons, "Coeffient of variation of detrended maize water footprints 1986-2016", CalculateBars(cv_wf_ir*10)/10, 50, '','RdYlBu_r', 1)
        
        _raster_path = ac.GetPath([self.raster_folder_path, f'acea_5arc_{self.conf.crop_name_short}_cv_wftot_global_annual_{self.text_years}.nc'])
        val_name = f'cv-{self.conf.crop_name_short}'
        title = f'CV of detrended {self.conf.crop_name}, Climate: {self.conf.climate_name}'
        if not os.path.exists(_raster_path): ac.CreateRaster5(cv_wf_tot, _raster_path, val_name, '%', title)
            
        # GW analysis
        wf_gw_detrended = self.wf_gw[0,:,:] + (self.wf_gw - wf_rf_trended); wf_gw_detrended[wf_gw_detrended<0] = 0
        std_wf_gw = np.ma.std(wf_gw_detrended, axis=0); mean_wf_gw = np.ma.mean(wf_gw_detrended, axis = 0); cv_wf_gw = std_wf_gw/mean_wf_gw*100
        
        # Global summary
        print("Tot: CV {:.1f}%, WF {:.1f} m3 t-1 y-1".format(np.ma.average(cv_wf_tot, weights = np.ma.mean(self.prod_sum[-5:,:,:], axis=0)), \
                                                  np.ma.average(np.ma.average(np.ma.reshape(self.wf_tot[-5:,:,:], (5,-1)), axis = 1, weights = np.ma.reshape(self.prod_sum[-5:,:,:], (5,-1))))
                                                  ))
        
        print("Rainfed: CV {:.1f}%, WF {:.1f} m3 t-1 y-1".format(np.ma.average(cv_wf_rf, weights = np.ma.mean(self.prod_rf[-5:,:,:], axis=0)), \
                                                  np.ma.average(np.ma.average(np.ma.reshape(self.wf_rf[-5:,:,:], (5,-1)), axis = 1, weights = np.ma.reshape(self.prod_rf[-5:,:,:], (5,-1))))
                                                  ))
        print("Irrigated: CV {:.1f}%, WF {:.1f} m3 t-1 y-1".format(np.ma.average(cv_wf_ir, weights = np.ma.mean(self.prod_ir[-5:,:,:], axis=0)), \
                                                  np.ma.average(np.ma.average(np.ma.reshape(self.wf_ir[-5:,:,:], (5,-1)), axis = 1, weights = np.ma.reshape(self.prod_ir[-5:,:,:], (5,-1))))
                                                  ))
        print("GW: CV {:.1f}%, WF {:.1f} m3 t-1 y-1".format(np.ma.average(cv_wf_gw, weights = np.ma.mean(self.prod_rf[-5:,:,:], axis=0)), \
                                                  np.ma.average(np.ma.average(np.ma.reshape(self.wf_gw[-5:,:,:], (5,-1)), axis = 1, weights = np.ma.reshape(self.prod_rf[-5:,:,:], (5,-1))))
                                                  ))
        
    #%% General functions
    def CalculateBars(self, data, start_0=True, quantiles = [30, 40, 50, 65, 95]):
        step = 8
        bar_l1 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[0]))
        bar_l2 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[1]))
        bar_l3 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[2]))
        bar_l4 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[3]))
        bar_l5 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[4]))
        
        start_val = 0 if start_0 else np.floor(np.nanpercentile(np.ma.filled(data, np.nan), 5))
        
        color_lvls = np.array([*np.arange(start_val, bar_l1*.9, np.floor(bar_l1/step)),
                 *np.arange(bar_l1, bar_l2*.9, np.floor(bar_l2/step)),
                 *np.arange(bar_l2, bar_l3*.9, np.floor(bar_l3/step)),
                 *np.arange(bar_l3, bar_l4*.9, np.floor(bar_l4/step)),
                 *np.arange(bar_l4, bar_l5, np.floor(bar_l5/step)),
                 ], dtype=int)
        
        return color_lvls

    def PlotMapWithDist(self, data, title, color_lvls, side_max, plot_label, colourmap='rainbow', plot_type=0, cut_low_lat = True, normilised = False):
        "Map with a latitudinal distribution side chart and a colorbar"
        
        # Process the data
        data_lat_med = np.ma.median(data, axis=1)
        data_lat_q10 = np.nanpercentile(np.ma.filled(data, np.nan), 10, axis=1)
        if cut_low_lat:
            data_lat_med[int(-0.29*self.res5[0]):] = np.ma.masked
            data_lat_q10[int(-0.29*self.res5[0]):] = np.ma.masked
        
        # Set up colours
        if max(color_lvls) > 250:
            bar_ticks = np.round(color_lvls[np.arange(0, len(color_lvls), 3)]/10)*10
        else:
            bar_ticks = np.round(color_lvls[np.arange(0, len(color_lvls), 3)]*10)/10
        
        cur_cmap = get_cmap(colourmap)
        if normilised:
            normi = mpl_colors.BoundaryNorm(color_lvls, cur_cmap.N, extend='both')
        
        # Create figure
        fig = plt.figure(figsize=(20,5), dpi=300)
        gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1], wspace=0.0, hspace=0.0, right=0.757)
        
        # Subplot 1 # cax = plt.axes(projection=ccrs.PlateCarree())
        ax1 = plt.subplot(gs[0], projection=ccrs.PlateCarree())
        if normilised:
            plot1 = plt.contourf(self.lons, self.lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, norm = normi, transform=ccrs.PlateCarree(), extend = 'max')
                # Colours: https://matplotlib.org/stable/tutorials/colors/colormaps.html
        else:
            plot1 = plt.contourf(self.lons, self.lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, transform=ccrs.PlateCarree(), extend = 'max')
                # Colours: https://matplotlib.org/stable/tutorials/colors/colormaps.html        
        
        ax1.set_xlim([-179, 179]); ax1.set_ylim([-60, 68.5])
        ax1.coastlines(resolution='110m') # Display the coastlines ("110m", "50m", and "10m")
        ax1.add_feature(cfeature.BORDERS, edgecolor='grey', alpha=0.4)
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='grey', alpha=0.3, linestyle='--', draw_labels=True)
        gl.top_labels = False; gl.right_labels=False
        # plt.title(title)
        
        # Subplot 2
        ax2 = plt.subplot(gs[1], sharey = ax1)
        plt.plot(data_lat_q10, self.lats, color='black', linewidth=1)
        plt.grid(linewidth=1, color='grey', alpha=0.3, linestyle='--')
        ax2.set_yticklabels([])
            
        ax2.fill_betweenx(self.lats, 0, data_lat_med, color='grey', alpha=0.5)
        # plt.title('Median latitudinal CV')
        ax2.set_xlim([0, side_max])
        
        # Final additions, show figure
        if plot_type == 1:
            ax2.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
            plt.colorbar(plot1, ax=ax2, shrink=.9, ticks = bar_ticks, format='%2.0f %%', label = plot_label) # https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.colorbar.html
        else:
            if max(bar_ticks) > 30:
                plt.colorbar(plot1, ax=ax2, shrink=.95, ticks = bar_ticks, format='%2.0f', label = plot_label)
            else:
                plt.colorbar(plot1, ax=ax2, shrink=.95, ticks = bar_ticks, format='%2.1f', label = plot_label)
        # cbar.ax.minorticks_on()
        fig.subplots_adjust(wspace=0, hspace=0)
        # fig.suptitle(title)
        plt.show()
        
    def PlotMapCustomLevels(self, data, title, color_lvls, percentage,plot_label, colourmap='rainbow', normilised = False, dpi_val=300):
        "General map with a colorbar"
        # Set up colours
        cur_cmap = get_cmap(colourmap)
        
        # Create figure
        plt.figure(figsize=(20,5), dpi=dpi_val)
        cur_cmap = get_cmap(colourmap)
        if normilised:
            normi = mpl_colors.BoundaryNorm(color_lvls, cur_cmap.N, extend='both')
        
        # Subplot 1
        cax = plt.axes(projection=ccrs.PlateCarree())
        if normilised:
            plot1 = plt.contourf(self.lons, self.lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, norm = normi, transform=ccrs.PlateCarree(), extend = 'both')
            # Colours: https://matplotlib.org/stable/tutorials/colors/colormaps.html
        else:
            plot1 = plt.contourf(self.lons, self.lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, transform=ccrs.PlateCarree(), extend = 'both')
            
        cax.set_xlim([-179, 179]); cax.set_ylim([-60, 68.5])
        cax.coastlines(resolution='110m') # Display the coastlines ("110m", "50m", and "10m")
        cax.add_feature(cfeature.BORDERS, edgecolor='grey', alpha=0.4)
        gl = cax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='grey', alpha=0.3, linestyle='--', draw_labels=True)
        gl.top_labels = False; gl.right_labels=False
        # plt.title(title)
        
        # Final additions, show figure
        if percentage:
            plt.colorbar(plot1, ax=cax, shrink=.9, pad=.015, format='%2.0f %%', label = plot_label) # https://matplotlib.org/devdocs/api/_as_gen/matplotlib.pyplot.colorbar.html
        else:
            plt.colorbar(plot1, ax=cax, shrink=.9, pad=.015, label = plot_label)
            
        plt.show()
        
    def TrendRasterData(self, data_ori, save=False, fpath=False, limit=True):
        "Detrend historical raster"
        
        cols = data_ori.shape[2]; years = np.arange(self.first_year, self.clock_end_year+1)
        empty_arr = np.zeros((data_ori.shape))
        data_trended = np.ma.masked_where(empty_arr == 0, empty_arr)
            
        for row in range(0, int(cols/2)):
            print('Detrending row %i' % row)
            for col in range(0, cols):
                data_temp = data_ori[:,row,col]
                if all(np.ma.getmask(data_temp)): continue
                if limit:
                    limit_up = np.percentile(data_temp,90); limit_lo = np.percentile(data_temp,10)
                    med_val = np.percentile(data_temp,50)
                    if limit_up > med_val*2: data_temp[data_temp>limit_up] = limit_up
                    if limit_lo < med_val/2: data_temp[data_temp<limit_lo] = limit_lo
                    
                coef, intercept = np.ma.polyfit(years, data_temp, 1)
                data_trended[:,row,col] = years*coef + intercept
        
        data_trended[data_trended<0] = 0
        
        if save:
            val_name = f'wf-{self.conf.crop_name_short}'; title = 'Trended WF total'
            title = f'Crop: {self.conf.crop_name}, Climate: {self.conf.climate_name}, Variable: crop yield scaled'
            ac.CreateHistRaster5(data_trended, fpath, val_name, 'm3 t-1 y-1', title)
            
        return data_trended