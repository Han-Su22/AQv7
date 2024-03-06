# -*- coding: utf-8 -*-
"""
Created on Fri May 28 10:00:26 2021

@author: MialykO
"""
#%% Import packages
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from modules.met_brewer import met_brew
import modules.acea.acea_core as ac
import importlib
import modules.footprints.footprints_core as fc
from scipy.stats import gaussian_kde
import scipy.stats as stats
from cartopy.mpl.gridliner import LATITUDE_FORMATTER
from matplotlib.cm import get_cmap
import matplotlib.colors as mpl_colors
from matplotlib import gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mtick

#%% Functions
def PlotMapWithDist(data, title, color_lvls, side_max, plot_label, colourmap='rainbow', plot_type=0, cut_low_lat = True, normilised = False):
    "Map with a latitudinal distribution side chart and a colorbar"
    
    # Process the data
    lats, lons, res5 = ac.GetResolution(1)
    data_lat_med = np.ma.median(data, axis=1)
    data_lat_q10 = np.nanpercentile(np.ma.filled(data, np.nan), 10, axis=1)
    if cut_low_lat:
        data_lat_med[int(-0.29*res5[0]):] = np.ma.masked
        data_lat_q10[int(-0.29*res5[0]):] = np.ma.masked
    
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
        plot1 = plt.contourf(lons, lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, norm = normi, transform=ccrs.PlateCarree(), extend = 'max')
            # Colours: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    else:
        plot1 = plt.contourf(lons, lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, transform=ccrs.PlateCarree(), extend = 'max')
            # Colours: https://matplotlib.org/stable/tutorials/colors/colormaps.html        
    
    ax1.set_xlim([-179, 179]); ax1.set_ylim([-60, 68.5])
    ax1.coastlines(resolution='110m') # Display the coastlines ("110m", "50m", and "10m")
    ax1.add_feature(cfeature.BORDERS, edgecolor='grey', alpha=0.4)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='grey', alpha=0.3, linestyle='--', draw_labels=True)
    gl.top_labels = False; gl.right_labels=False
    # plt.title(title)
    
    # Subplot 2
    ax2 = plt.subplot(gs[1], sharey = ax1)
    plt.plot(data_lat_q10, lats, color='black', linewidth=1)
    plt.grid(linewidth=1, color='grey', alpha=0.3, linestyle='--')
    ax2.set_yticklabels([])
        
    ax2.fill_betweenx(lats, 0, data_lat_med, color='grey', alpha=0.5)
    # plt.title('Median latitudinal CV')
    ax2.set_xlim([0, side_max])
    
    # Final additions, show figure
    if plot_type == 1: # Percentage
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
    
def PlotMapCustomLevels(data, title, color_lvls, percentage, plot_label, colourmap='rainbow', normilised = False, dpi_val=300):
    "General map with a colorbar"
    # Set up colours
    cur_cmap = get_cmap(colourmap)
    lats, lons, res5 = ac.GetResolution(1)
    
    # Create figure
    plt.figure(figsize=(20,5), dpi=dpi_val)
    cur_cmap = get_cmap(colourmap)
    if normilised:
        normi = mpl_colors.BoundaryNorm(color_lvls, cur_cmap.N, extend='both')
    
    # Subplot 1
    cax = plt.axes(projection=ccrs.PlateCarree())
    if normilised:
        plot1 = plt.contourf(lons, lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, norm = normi, transform=ccrs.PlateCarree(), extend = 'both')
        # Colours: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    else:
        plot1 = plt.contourf(lons, lats, data, levels = color_lvls, alpha = .9, cmap=cur_cmap, transform=ccrs.PlateCarree(), extend = 'both')
        
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
    
def PlotLatDist_unitChange(crop_num, crop_names, harv_area_avg, f_avg_new, f_avg_old, lats, title, save_path=False, step=10):
    "Latitudinal distribution of footprints for multiple crops"
    
    fig, axs = plt.subplots(1, crop_num, sharey=True, figsize=(2*crop_num, 6), dpi=300, )
    fig.subplots_adjust(top=0.85)
    
    c = 0
    for ax in axs.flatten():
        # Calculate numbers
        _threshold = np.nanpercentile(np.ma.filled(harv_area_avg[c], np.nan), 10)
        _f_change = np.around((f_avg_new[c]/f_avg_old[c]-1)*100,1)
        _f_change[harv_area_avg[c]<=_threshold] = np.ma.masked; _f_change = np.ma.median(_f_change, axis=1)
        _ha_avg = np.ma.sum(harv_area_avg[c], axis=1)
        
        # Primary axis for 
        data_lat_med = np.array([np.mean(np.ma.filled(_f_change[i:i+step],0)) for i in range(0, len(_f_change), step)])
        data_lat_med_pos = data_lat_med*1; data_lat_med_pos[data_lat_med_pos > 0] = 0
        data_lat_med_neg = data_lat_med*1; data_lat_med_neg[data_lat_med_neg < 0] = 0

        ax.plot([0,0], [lats[0],lats[-1]], c='silver', ls="--", lw=1)
        ax.fill_betweenx(lats[::step], 0, data_lat_med_pos, color='teal', alpha=0.6, linewidth=0)
        ax.fill_betweenx(lats[::step], data_lat_med_neg, 0, color='tomato', alpha=0.6, linewidth=0)
        
        # Secondary axis for HA
        data = np.array([np.nansum(np.ma.filled(_ha_avg[i:i+step],0))/10**6 for i in range(0, len(_ha_avg), step)]) # Sum every X elements
        ax_top = ax.twiny() # add axis
        ax_top.fill_betweenx(lats[::step], 0, data, color='silver', alpha=0.5, linewidth=0)
        ax_top.fill_betweenx(lats[::step], -data, 0, color='silver', alpha=0.5, linewidth=0)
        ax_top.spines['right'].set_visible(False); ax_top.spines['left'].set_visible(False)
        ax_top.set_xlabel('Harvested area (Mha)'); _max=max(data)
        ax_top.set_xlim(-_max, _max)
        ax_top.set_xticks([-int(_max/1.5), 0, int(_max/1.5)])
        ax_top.set_xticklabels([int(_max/1.5), 0, int(_max/1.5)])
        
        # Styling
        ax.set_xlim([-80,80]); ax.grid(axis = 'y', alpha=0.5)
        ax.spines['right'].set_visible(False);ax.spines['top'].set_visible(False)
        if c == 0:
            ax.set_ylabel('Latitude')
            ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
            ax.spines.left.set_position(("axes", -0.1))
        else:
            ax.tick_params(left = False)
            ax.spines['left'].set_visible(False)
            
        ax.set_title(f'{crop_names[c]}')
        ax.set_xlabel('Change (%)'); # ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
        c += 1
    
    plt.ylim([-50, 80]); fig.suptitle(title)
    
    if save_path: plt.savefig(save_path, bbox_inches ="tight") #transparent = True,
    plt.show()

def PlotDens_unitChange(globalF, xmax, crop_num, crop_names, f_avg_new, f_wavg_new, f_avg_old, \
                        f_wavg_old, prodF_new, colours, xlabel, save_path=False):
    "Change of density distribution of footprints for multiple crops"
    
    x_vals =  [*range(0,xmax,int(xmax/200))] # Specifying the limits of our data
    
    fig, axs = plt.subplots(crop_num, 1, sharex=True, figsize=(10, 1.4*crop_num), dpi=300, )
    c = 0
    for ax in axs.flatten():
        # 2015-2019
        _f_avg_new = f_avg_new[c,:,:]; _f_avg_old = f_avg_old[c,:,:];
        F_prod = np.around(np.ma.sum(prodF_new[c,:,:])/10**9 / globalF*100,1)
        _threshold = np.nanpercentile(np.ma.filled(_f_avg_new, np.nan), 95); _f_avg_new[_f_avg_new>=_threshold] = np.ma.masked
        _f_avg_new = np.array(_f_avg_new[_f_avg_new.mask==False].flatten()); kde = stats.gaussian_kde(_f_avg_new)
        _f_avg_new = gaussian_kde(_f_avg_new, bw_method=kde.factor / 2.5); _f_avg_new = _f_avg_new(x_vals)
        ax.fill_between(x_vals, _f_avg_new, color=colours[c], alpha = 0.8, linewidth=0) # Add distribution
        ax.plot(f_wavg_new[c], np.max(_f_avg_new)/2, 'wo', mec='black', ms=10) # Add global average
        
        # 1990-1994
        _threshold = np.nanpercentile(np.ma.filled(_f_avg_old, np.nan), 95); _f_avg_old[_f_avg_old>=_threshold] = np.ma.masked
        _f_avg_old = np.array(_f_avg_old[_f_avg_old.mask==False].flatten()); kde = stats.gaussian_kde(_f_avg_old)
        _f_avg_old = gaussian_kde(_f_avg_old, bw_method=kde.factor / 2.5); _f_avg_old = _f_avg_old(x_vals)
        ax.fill_between(x_vals, _f_avg_old, color=colours[c], alpha = 0.4, linewidth=0) # Add distribution
        ax.plot(f_wavg_old[c], np.max(_f_avg_new)/2, 'o', c='black', ms=7) # Add global average
        
        # Add text
        unitF_change = np.around((f_wavg_new[c]/f_wavg_old[c]-1)*100,1)
        unitF_change = f'+{unitF_change}' if unitF_change > 0 else unitF_change
        
        if c == 0:
            ax.annotate('2017-2019', xy=(f_wavg_new[c], np.max(_f_avg_new)/2), xytext=(f_wavg_new[c]*1.4, np.max(_f_avg_new)/1.1), fontsize=10, color='grey',
                        arrowprops=dict(arrowstyle="->", linestyle='--', color='silver', shrinkA=2, shrinkB=5, \
                                lw=1, connectionstyle="arc3,rad=0.3"))
            ax.annotate('1990-1992', xy=(f_wavg_old[c], np.max(_f_avg_new)/2), xytext=(f_wavg_old[c]*1.3, np.max(_f_avg_new)/1.8), fontsize=10, color='grey',
                        arrowprops=dict(arrowstyle="->", linestyle='--', color='silver', shrinkA=2, shrinkB=5, \
                                lw=1, connectionstyle="arc3,rad=-0.3"))
                
            ax.annotate('Historical change (%)', xy=(xmax*0.88, np.max(_f_avg_new)/1.3), xytext=(xmax*0.62, np.max(_f_avg_new)/1.3), fontsize=10, color='grey',
                        arrowprops=dict(arrowstyle="->", linestyle='--', color='silver', shrinkA=5, shrinkB=15, \
                                lw=1, connectionstyle="arc3,rad=-0.3"))
            ax.annotate('Global share (%)', xy=(xmax*0.95, np.max(_f_avg_new)/1.3), xytext=(xmax*0.8, np.max(_f_avg_new)/4.2), fontsize=10, color='grey',
                        arrowprops=dict(arrowstyle="->", linestyle='--', color='silver', shrinkA=5, shrinkB=15, \
                                lw=1, connectionstyle="arc3,rad=0.3"))
            
        ax.annotate('', xy=(f_wavg_new[c], np.max(_f_avg_new)/2), xytext=(f_wavg_old[c], np.max(_f_avg_new)/2),
                        arrowprops=dict(arrowstyle="-|>", color='black', shrinkB=5, lw=1.5))
        
        ax.scatter(xmax*0.88, np.max(_f_avg_new)/1.3, s=750, facecolors='none', edgecolors='grey')
        ax.scatter(xmax*0.95, np.max(_f_avg_new)/1.3, s=750, facecolors='none', edgecolors='silver')
        ax.text(xmax*0.88, np.max(_f_avg_new)/1.3, f'{unitF_change}', color='grey', fontsize=8, ha='center', va='center')
        ax.text(xmax*0.95, np.max(_f_avg_new)/1.3, f'{F_prod}', color='grey', fontsize=8, ha='center', va='center')
        
        # Styling
        ax.grid(axis = 'x', alpha=0.5); ax.set_yticks([])
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.set_ylabel(f'{crop_names[c]}')
        
        c += 1
        
    ax.set_xlabel(xlabel); plt.xlim([0, xmax])
    
    if save_path: plt.savefig(save_path, bbox_inches ="tight") #transparent = True,
    plt.show()
    
def PlotTimeseries_WF(crops):
    "Historical average footprints for multiple crops"
    
    crop_list = fc.GetSimulatedCrops()
    folder = ac.GetPath(["outputs", 'main'])
    global_data= ac.GetPath([folder, f'GlobalStats_All_1990_2019.txt'])
    global_unitwf = pd.read_csv(global_data, skiprows=34, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)
    global_ha_2019 = pd.read_csv(global_data, skiprows=313, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.) + \
                    pd.read_csv(global_data, skiprows=344, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)
    global_prod_2019 = pd.read_csv(global_data, skiprows=3, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)
    global_ha_2019_top50 = global_ha_2019.tail(1).T.sort_values(by=[2019], ascending=False).head(50)
    global_unitwf_change = global_unitwf*1
    
    for name, values in global_unitwf.iteritems():
        if global_unitwf_change[name][1990] == 0: global_unitwf_change[name] = 0.
        else:
            avg_unitwf = global_unitwf_change[name][[1990,1991,1992]].sum()/3
            global_unitwf_change[name] = (global_unitwf_change[name]/avg_unitwf-1)
    
    
    colours = met_brew(name="Cross", n=len(crop_list['Group'].unique()), brew_type="continuous")
    
    plt.figure(figsize=(10,8), dpi=300)
    count=0; group_names = []
    for group in crop_list['Group'].unique():
        crops = crop_list[crop_list['Group'] == group]['ACEA'].to_list()

        crops_max = global_unitwf_change[crops].T.quantile(.75, axis = 0).to_list()
        crops_avg_w = ((global_unitwf_change[crops].T * global_prod_2019[crops].T)).sum() / global_prod_2019[crops].T.sum().to_list()
        crops_min = global_unitwf_change[crops].T.quantile(.25, axis = 0).to_list()
        # crops = list(set(crops) & set(global_ha_2019_top20.index.to_list()))
        if len(crops)>0:
            plt.fill_between(global_unitwf_change.index, crops_max, crops_min, color=colours[count], alpha = .4, linewidth=0)
            plt.plot(global_unitwf_change.index,crops_avg_w, color=colours[count])
            group_names.append(group)
        count+=1; 
    
    plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0))
    plt.legend(group_names, loc='upper left')
    plt.grid(alpha=0.5); plt.show()
    
    a=1
    
    # fig, (ax0,ax1,ax2) = plt.subplots(3, 1, sharex=True, figsize=(10, 12), dpi=300)
    # plt_range = [*range(clock_start_year, clock_end_year+1)]
    
    # # Plot unit WFs
    # for c in range(crop_num):
    #     data = np.ma.average(unitF_tot[c,:,:,:], axis=(1,2),  weights = prod_tot[c,:,:,:])
    #     ax0.plot(plt_range, data, c=colours[c], lw=2)
    
    # ax0.set_ylabel(unitlabel)
    # ax0.grid(alpha=0.5); ax0.spines['top'].set_visible(False); ax0.spines['right'].set_visible(False)
    
    # # Plot production WFs
    # data = np.zeros((crop_num,len(plt_range)))
    # for c in range(crop_num):
    #     data[c,:] = np.ma.sum(unitF_tot[c,:,:,:]*prod_tot[c,:,:,:], axis=(1,2))
    #     ax1.plot(plt_range, data[c,:]/10**9, c=colours[c], lw=2)
        
    # ax1.set_ylabel(prodlabel)
    # ax1.grid(alpha=0.5); ax1.spines['top'].set_visible(False); ax1.spines['right'].set_visible(False)
    
    # # Plot cummulative production WFs
    # data_cum = np.cumsum(data, axis=0)
    # for c in range(crop_num):
    #     if c == 0: ax2.fill_between(plt_range, data_cum[c,:]/10**9, color=colours[c])
    #     else: ax2.fill_between(plt_range, data_cum[c,:]/10**9, data_cum[c-1,:]/10**9, color=colours[c])
    
    # ax2.set_ylabel(prodlabel)
    # ax2.grid(alpha=0.5); ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
    
    # fig.legend(crop_names, loc='center right', ncol=1)
    # plt.subplots_adjust(right=0.85)
    
    # if save_path: plt.savefig(save_path, transparent = True,bbox_inches ="tight")
    # plt.show()

#%% Classes
class wf_visualiser:
    data_folder = ac.GetPath(["data","acea"])
    lats, _, _ = ac.GetResolution(1)
    def __init__(self, crops, save_name_add, colour_palette, years=[1990,2019], years_avg=5):
        self.crop_num = len(crops); self.crop_names = crops['Name'].to_list()
        self.colours = met_brew(name=colour_palette, n=self.crop_num, brew_type="continuous")
        self.crops = list(crops.index.to_list())
        self.save_name_add = save_name_add; self.years_avg = years_avg
        self.clock_start_year = years[0]; self.clock_end_year = years[1]
        self.max_gs = (self.clock_end_year - self.clock_start_year)+1
        
        print("============Visualisation of crop WFs============")
        print("Getting main outputs...")
        self.GetMainOutputs() # Get main ACEA outputs
        
    def GetMainOutputs(self):
        self.prodF_new = np.zeros((self.crop_num, 2160, 4320))*np.ma.masked
        self.unitWF_avg_old = self.prodF_new*1; self.unitWF_avg_new =  self.prodF_new*1; self.harv_area_avg =  self.prodF_new*1
        self.unitWF_wavg_old = np.zeros(self.crop_num); self.unitWF_wavg_new = self.unitWF_wavg_old*1
        
        c_count = 0
        for c in self.crops:
            project = f"{c}_phd_mialyk_2022"
            raster_folder_path =  ac.GetPath(["outputs", project,'rasters'])
            
            # Get production
            _raster_path = ac.GetPath([raster_folder_path, f'acea_5arc_{c}_prod_ir_global_annual_{self.clock_start_year}_{self.clock_end_year}.nc'])
            _prod_tot = np.ma.filled(ac.ReadNC(_raster_path, f'production-{c}'), 0)
            _raster_path = ac.GetPath([raster_folder_path, f'acea_5arc_{c}_prod_rf_global_annual_{self.clock_start_year}_{self.clock_end_year}.nc'])
            _prod_tot = _prod_tot + np.ma.filled(ac.ReadNC(_raster_path, f'production-{c}'), 0)
            _prod_tot = np.ma.masked_where(_prod_tot == 0, _prod_tot)
                        
            # Get WF
            _raster_path = ac.GetPath([raster_folder_path, f'acea_5arc_{c}_wf_tot_global_annual_{self.clock_start_year}_{self.clock_end_year}.nc'])
            _data_WF = ac.ReadNC(_raster_path, f'wf_total-{c}')
            _data_WF[(_data_WF < np.nanpercentile(np.ma.filled(_data_WF, np.nan), 1)) & (_data_WF > np.nanpercentile(np.ma.filled(_data_WF, np.nan), 99))] = np.ma.masked
            
            
            self.unitWF_wavg_old[c_count] = np.ma.mean(np.ma.average(_data_WF[:self.years_avg,:,:], axis=(1,2), weights = _prod_tot[:self.years_avg,:,:])) # Old weighted avg
            self.unitWF_wavg_new[c_count] = np.ma.mean(np.ma.average(_data_WF[-self.years_avg:,:,:], axis=(1,2), weights = _prod_tot[-self.years_avg:,:,:])) # Recent weighted avg
            self.unitWF_avg_old[c_count, :, :] = np.ma.median(_data_WF[:self.years_avg,:,:], axis=0) # Old avg
            self.unitWF_avg_new[c_count, :, :] = np.ma.median(_data_WF[-self.years_avg:,:,:], axis=0) # Recent avg
            self.prodF_new[c_count, :, :] = _prod_tot[-1,:,:] * _data_WF[-1,:,:] # Recent prod avg
            
            # Get harvested areas
            crop_fao = getattr(importlib.import_module(f'projects.{project}'), 'project_conf')# Import project data
            crop_fao = crop_fao.crop_fao
            harv_area_rf = np.ma.filled(fc.GetHistHarvestedAreas5(crop_fao, False),0); harv_area_ir = np.ma.filled(fc.GetHistHarvestedAreas5(crop_fao, True),0)
            harv_area_tot = harv_area_ir + harv_area_rf; _harv_area_tot = np.ma.masked_where(harv_area_tot == 0, harv_area_tot)
            self.harv_area_avg[c_count, :, :] = np.ma.mean(_harv_area_tot[-self.years_avg:,:,:], axis=0)
            
            c_count +=1
            
    # def PlotTimeseries_WF(self):
    #     print('- Plotting historical timeseries')
    #     fc.PlotTimeseries_prod_unit(self.clock_start_year, self.clock_end_year, self.crop_num, self.crop_names, self.unitWF_tot, self.prod_tot, self.colours, \
    #                       'Unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)', 'Production water footprint ($10^9 m^3$ $y^{-1}$)')
    # def PlotTimeseries_unitWF(self):
       
    #     fig, ax1 = plt.subplots(1, 1,figsize=(10, 4), dpi=300)
    #     plt_range = [*range(self.clock_start_year,self.clock_end_year+1)]
        
    #     for c in range(self.crop_num):
    #         data = np.ma.average(self.unitWF_tot[c,:,:,:], axis=(1,2),  weights = self.prod_tot[c,:,:,:])
    #         ax1.plot(plt_range, data)
        
    #     ax1.set_ylabel('Unit water footprint ($m^3$ $t^{-1}$ $y^{-1}$)')
    #     fig.legend(self.crop_names, loc='center right', ncol=1)
    #     ax1.grid()
    #     plt.show()
    def CalculateBars(self, data, start_0=True, quantiles = [30, 40, 50, 65, 95]):
        step = 8
        bar_l1 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[0]))
        bar_l2 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[1]))
        bar_l3 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[2]))
        bar_l4 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[3]))
        bar_l5 = np.floor(np.nanpercentile(np.ma.filled(data, np.nan), quantiles[4]))
        
        start_val = 0 if start_0 else np.floor(np.nanpercentile(np.ma.filled(data, np.nan), 2))
        
        color_lvls = np.array([*np.arange(start_val, bar_l1*.9, np.floor(bar_l1/step)),
                 *np.arange(bar_l1, bar_l2*.9, np.floor(bar_l2/step)),
                 *np.arange(bar_l2, bar_l3*.9, np.floor(bar_l3/step)),
                 *np.arange(bar_l3, bar_l4*.9, np.floor(bar_l4/step)),
                 *np.arange(bar_l4, bar_l5, np.floor(bar_l5/step)),
                 ], dtype=int)
        
        return color_lvls
    
    def PlotMaps_WF(self):
        print('- Plotting maps')
        for c in range(self.crop_num):
            data = self.unitWF_avg_new[c,:,:]*1
            PlotMapWithDist(data, " ", self.CalculateBars(data, False), int(np.nanpercentile(np.ma.filled(data, np.nan), 90)), f'Unit water footprint of {self.crop_names[c]} ($m^3$ $t^{-1}$ $y^{-1}$)','Spectral_r',0, True, True)
            PlotMapCustomLevels((data/self.unitWF_avg_old[c,:,:] - 1) * 100, " ", np.arange(-80, 81, 10), True, f"Change in unit water footprint of {self.crop_names[c]}", 'Spectral_r')
            
    def PlotLatDist_unitWF(self):
        print('- Plotting latitudinal distribution')
        save_path = ac.GetPath(["outputs", 'main', 'figures', f'{self.save_name_add}_unitWF_lat.png'])
        PlotLatDist_unitChange(self.crop_num, self.crop_names, self.harv_area_avg, self.unitWF_avg_new, self.unitWF_avg_old, self.lats, 'Unit water footprint', save_path)
        
    def PlotDens_unitWF(self, globalF, xmax):
        print('- Plotting density distribution')
        save_path = ac.GetPath(["outputs", 'main', 'figures', f'{self.save_name_add}_unitWF_dens.png'])
        PlotDens_unitChange(globalF, xmax, self.crop_num, self.crop_names, self.unitWF_avg_new, self.unitWF_wavg_new, \
                               self.unitWF_avg_old, self.unitWF_wavg_old, self.prodF_new, self.colours, 'Unit water footprint ($m^3$ $t^{-1}$ $yr^{-1}$)', save_path)

class lf_visualiser:
    data_folder = ac.GetPath(["data","acea"])
    lats, _, _ = ac.GetResolution(1)
    
    def __init__(self, crops, years, years_avg, save_name_add):
        self.colours = list(np.array(list(crops.values()))[:,1])
        self.crop_names = list(np.array(list(crops.values()))[:,0])
        self.crops = list(crops.keys()); self.crop_num = len(crops)
        self.save_name_add = save_name_add
        self.clock_start_year = years[0]; self.clock_end_year = years[1]
        self.max_gs = (self.clock_end_year - self.clock_start_year)+1
        self.years_avg = years_avg
        print("============Visualisation of crop LFs============")
        print("Getting main outputs...")
        self.GetMainOutputs() # Get main ACEA outputs
        
    def GetMainOutputs(self):
        self.prodF_new = np.zeros((self.crop_num, 2160, 4320))*np.ma.masked
        self.unitLF_avg_old = self.prodF_new*1; self.unitLF_avg_new =  self.prodF_new*1; self.harv_area_avg =  self.prodF_new*1
        self.unitLF_wavg_old = np.zeros(self.crop_num); self.unitLF_wavg_new = self.unitLF_wavg_old*1
        
        c_count = 0
        for c in self.crops:
            project = f"{c}_phd_mialyk_2022"
            raster_folder_path =  ac.GetPath(["outputs", project,'rasters'])
            
            # Get production
            _raster_path = ac.GetPath([raster_folder_path, f'acea_5arc_{c}_prod_ir_global_annual_{self.clock_start_year}_{self.clock_end_year}.nc'])
            _prod_tot = np.ma.filled(ac.ReadNC(_raster_path, f'production-{c}'), 0)
            _raster_path = ac.GetPath([raster_folder_path, f'acea_5arc_{c}_prod_rf_global_annual_{self.clock_start_year}_{self.clock_end_year}.nc'])
            _prod_tot = _prod_tot + np.ma.filled(ac.ReadNC(_raster_path, f'production-{c}'), 0)
            _prod_tot = np.ma.masked_where(_prod_tot == 0, _prod_tot)
            
            # Get WF
            _raster_path = ac.GetPath([raster_folder_path, f'acea_5arc_{c}_lf_tot_global_annual_{self.clock_start_year}_{self.clock_end_year}.nc'])
            _data_LF = ac.ReadNC(_raster_path, f'lf_total-{c}') 
            
            self.unitLF_wavg_old[c_count] = np.ma.mean(np.ma.average(_data_LF[:self.years_avg,:,:], axis=(1,2), weights = _prod_tot[:self.years_avg,:,:])) # Old weighted avg
            self.unitLF_wavg_new[c_count] = np.ma.mean(np.ma.average(_data_LF[-self.years_avg:,:,:], axis=(1,2), weights = _prod_tot[-self.years_avg:,:,:])) # Recent weighted avg
            self.unitLF_avg_old[c_count, :, :] = np.ma.mean(_data_LF[:self.years_avg,:,:], axis=0) # Old avg
            self.unitLF_avg_new[c_count, :, :] = np.ma.mean(_data_LF[-self.years_avg:,:,:], axis=0) # Recent avg
            self.prodF_new[c_count, :, :] = _prod_tot[-1,:,:] * _data_LF[-1,:,:] # Recent prod avg
            
            # Get harvested areas
            crop_fao = getattr(importlib.import_module(f'projects.{project}'), 'project_conf')# Import project data
            crop_fao = crop_fao.crop_fao
            harv_area_rf = np.ma.filled(fc.GetHistHarvestedAreas5(crop_fao, False),0); harv_area_ir = np.ma.filled(fc.GetHistHarvestedAreas5(crop_fao, True),0)
            harv_area_tot = harv_area_ir + harv_area_rf; _harv_area_tot = np.ma.masked_where(harv_area_tot == 0, harv_area_tot)
            self.harv_area_avg[c_count, :, :] = np.ma.mean(_harv_area_tot[-self.years_avg:,:,:], axis=0)
            
            c_count +=1
    
    # def PlotTimeseries_LF(self):
    #     print('- Plotting historical timeseries')
    #     fc.PlotTimeseries_prod_unit(self.clock_start_year, self.clock_end_year, self.crop_num, self.crop_names, self.unitLF_tot, self.prod_tot, self.colours, \
    #                       'Unit land footprint ($m^2$ $t^{-1}$ $y^{-1}$)', 'Production land footprint ($10^9 m^2$ $y^{-1}$)')
        
    def PlotLatDist_unitLF(self):
        print('- Plotting latitudinal distribution')
        save_path = ac.GetPath(["outputs", 'main', 'figures', f'{self.save_name_add}_unitLF_lat.png'])
        fc.PlotLatDist_unitChange(self.crop_num, self.crop_names, self.harv_area_avg, self.unitLF_avg_new, self.unitLF_avg_old, self.lats, 'Unit land footprint', save_path)
        
    def PlotDens_unitLF(self, globalF, xmax):
        print('- Plotting density distribution')
        save_path = ac.GetPath(["outputs", 'main', 'figures', f'{self.save_name_add}_unitLF_dens.png'])
        fc.PlotDens_unitChange(globalF, xmax, self.crop_num, self.crop_names,self.unitLF_avg_new, self.unitLF_wavg_new, \
                               self.unitLF_avg_old, self.unitLF_wavg_old, self.prodF_new, self.colours, 'Unit land footprint ($m^2$ $t^{-1}$ $yr^{-1}$)', save_path)