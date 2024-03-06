# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 08:43:42 2022

@author: MialykO
"""

#%% Import packages
import numpy as np, os, importlib, pandas as pd
from  datetime import datetime
import modules.acea.acea_core as ac, modules.footprints.footprints_core as fc, pandas as pd

#%% Classes
class footprint_reporter:
    "Analyse multiple crops at once"
    
    def __init__(self, crops, years, title):
        self.crops_name = crops['Name'].values.tolist()
        self.crops_acea = crops['ACEA'].values.tolist()
        self.crops_fao = crops['FAO'].values.tolist()
        self.crops_groups = crops['Group'].values.tolist()
        
        self.crop_num = len(crops)
        self.clock_start_year = years[0]; self.clock_end_year = years[1]
        self.gs_num = (self.clock_end_year - self.clock_start_year)+1
        self.title = title
        self.text_years = f'{self.clock_start_year}_{self.clock_end_year}'
        
    def GlobalAnalysis_unitLF(self):
        "Txt annual report of LF analysis"
        
        # 1. Create report file
        out_report = ac.GetPath(["outputs", 'general'])
        if not os.path.exists(out_report): os.makedirs(out_report)
        out_report = ac.GetPath([out_report, f'GlobalAnalysis_LF_{self.text_years}.txt'])
        if os.path.exists(out_report): os.remove(out_report)
        ac.printLog(out_report, datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        
        # 2. Collect data
        LF = np.zeros((self.crop_num, self.gs_num, 8))
        c_count = 0
        for c in self.crops_acea:
            project = f"{c}_phd_mialyk_2022"
            rep_folder_path =  ac.GetPath(["outputs", project,'reports'])
            crop_fao = getattr(importlib.import_module(f'projects.{project}'), 'project_conf')# Import project data
            crop_fao = crop_fao.crop_fao
            
            # Get LF
            _raster_path = ac.GetPath([rep_folder_path, f'{crop_fao}_LF_GlobalStats_{self.text_years}.txt'])
            _data = pd.read_csv(_raster_path, skiprows=65, sep='\t', nrows=30)
            _harv_area = np.around(_data['Harvested area[ha]'].values/10**5,2) # Harvested area [Mha]
            
            LF[c_count,:,0] = _data['Unit LF[m2 t-1]'].values # Unit LF avg [m2 t-1]
            LF[c_count,:,1] = _data['AvgS[%]'].values # Weighted suitability
            LF[c_count,:,2] = np.around(_data['HA_HS'].values*_harv_area/100,2) # HS area[Mha]
            LF[c_count,:,3] = np.around(_data['HA_MS'].values*_harv_area/100,2) # MS area[Mha]
            LF[c_count,:,4] = np.around(_data['HA_PS'].values*_harv_area/100,2) # PS area[Mha]
            LF[c_count,:,5] = _data['LF_HS'].values # Unit LF avg for HS [m2 t-1]
            LF[c_count,:,6] = _data['LF_MS'].values # Unit LF avg for MS [m2 t-1]
            LF[c_count,:,7] = _data['LF_PS'].values # Unit LF avg for PS [m2 t-1]
            
            c_count +=1
        
        # 3. Write report
        # LF HS
        ac.printLog(out_report, '\tLF HS'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,5])))

        # LF MS
        ac.printLog(out_report, '\tLF MS'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,6])))

        # LF PS
        ac.printLog(out_report, '\tLF PS')
        ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,7])))

        # LF
        ac.printLog(out_report, '\tUnit LF[m2 t-1]'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,0])))

        # Suitability HS
        ac.printLog(out_report, '\tHS[Mha]'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num):  ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,2])))
        
        # Suitability MS
        ac.printLog(out_report, '\tMS[Mha]'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,3])))

        # Suitability PS
        ac.printLog(out_report, '\tPS[Mha]'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,4])))
        
        # Suitability avg abs
        ac.printLog(out_report, '\tAvg Suitability'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,1])))
        
    def GlobalReport_CropSummary(self, avg_year=3):
        "Summary report of global values"
        
        # Extract data
        folder = ac.GetPath(["outputs", 'general'])
        _crop_report = ac.GetPath([folder, f'GlobalStats_All_{self.text_years}.txt'])
        WFu = pd.read_csv(_crop_report, skiprows=34, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        WFu_g = pd.read_csv(_crop_report, skiprows=65, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        crop_yield = pd.read_csv(_crop_report, skiprows=189, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        ha_rf =  pd.read_csv(_crop_report, skiprows=313, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)*10**4
        ha_ir = pd.read_csv(_crop_report, skiprows=344, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)*10**4
        prod = pd.read_csv(_crop_report, skiprows=3, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)
        wf_prod = WFu * prod
        LFu = pd.read_csv(_crop_report, skiprows=158, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        lf_prod = ha_rf+ha_ir
        y_rf = pd.read_csv(_crop_report, skiprows=375, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        y_ir = pd.read_csv(_crop_report, skiprows=406, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        
        # Cooking data
        avg_wf_prod = wf_prod.tail(avg_year).mean(axis=0).values
        avg_lf_prod = lf_prod.tail(avg_year).mean(axis=0).values
        lfu_change = (LFu.tail(avg_year).mean(axis=0)/LFu.head(avg_year).mean(axis=0)).values-1; lfu_change[lfu_change==np.inf] = 0
        world_ir_lf = ha_ir.tail(avg_year).sum().sum()/(ha_rf+ha_ir).tail(avg_year).sum().sum()
        world_wfg_pct = (WFu_g * prod).tail(avg_year).sum().sum()/wf_prod.tail(avg_year).sum().sum()
        wfu_change = (WFu.tail(avg_year).mean(axis=0)/WFu.head(avg_year).mean(axis=0)).values-1; wfu_change[wfu_change==np.inf] = 0
        world_wfu_change = (wfu_change*prod.tail(avg_year).mean(axis=0).values).sum()/prod.tail(avg_year).mean(axis=0).sum()
        world_lfu_change = (lfu_change*prod.tail(avg_year).mean(axis=0).values).sum()/prod.tail(avg_year).mean(axis=0).sum()
        wft_change = avg_wf_prod/wf_prod.head(avg_year).mean(axis=0).values-1; wft_change[wft_change==np.inf] = 0
        lft_change = avg_lf_prod/lf_prod.head(avg_year).mean(axis=0).values-1; lft_change[lft_change==np.inf] = 0
        world_wft_change = avg_wf_prod.sum()/wf_prod.head(avg_year).mean(axis=0).sum()-1
        world_lft_change = avg_lf_prod.sum()/lf_prod.head(avg_year).mean(axis=0).sum()-1
        
        ha_rf_change = ha_rf.tail(avg_year).mean(axis=0).values/ha_rf.head(avg_year).mean(axis=0).values-1; ha_rf_change[ha_rf_change==np.inf] = 0
        ha_ir_change = ha_ir.tail(avg_year).mean(axis=0).values/ha_ir.head(avg_year).mean(axis=0).values-1; ha_ir_change[ha_ir_change==np.inf] = 0
        world_ha_rf_change = ha_rf.tail(avg_year).mean(axis=0).sum()/ha_rf.head(avg_year).mean(axis=0).sum()-1
        world_ha_ir_change = ha_ir.tail(avg_year).mean(axis=0).sum()/ha_ir.head(avg_year).mean(axis=0).sum()-1
         
        y_rf_change = (y_rf.tail(avg_year).mean(axis=0)/y_rf.head(avg_year).mean(axis=0)).values-1; y_rf_change[y_rf_change==np.inf] = 0
        y_ir_change = (y_ir.tail(avg_year).mean(axis=0)/y_ir.head(avg_year).mean(axis=0)).values-1; y_ir_change[y_ir_change==np.inf] = 0
        
        # Final dataframe
        col_names = ['Crop', 'Group', '% of total WF', '% of total LF', 'Irrigated area','Crop yield[t ha-1]', 'WFu[m3 t-1]', 'Green', 'Blue', 
                     'WFu change','WFt change','LFu change', 'LFt change', 'HArf change','HAir change', 'Yrf change','Yir change']
        final_df = pd.DataFrame(columns = col_names, index=self.crops_acea)
        
        # Data per crop
        final_df['Crop'] = self.crops_name; final_df['Group'] = self.crops_groups
        final_df['% of total WF'] = avg_wf_prod/avg_wf_prod.sum()
        final_df['% of total LF'] = avg_lf_prod/avg_lf_prod.sum()
        final_df['Irrigated area'] = (ha_ir/(ha_rf+ha_ir)).tail(avg_year).mean(axis=0).values
        final_df['Crop yield[t ha-1]'] = crop_yield.tail(avg_year).mean(axis=0).values
        final_df['WFu[m3 t-1]'] = WFu.tail(avg_year).mean(axis=0).values
        final_df['Green'] = (WFu_g/WFu).tail(avg_year).mean(axis=0).values
        final_df['Blue'] = 1-final_df['Green']
        final_df['WFu change'] = wfu_change
        final_df['WFt change'] = wft_change
        final_df['LFu change'] = lfu_change
        final_df['LFt change'] = lft_change
        final_df['HArf change'] = ha_rf_change
        final_df['HAir change'] = ha_ir_change
        final_df['Yrf change'] = y_rf_change
        final_df['Yir change'] = y_ir_change
        
        # Data per crop group
        for group in set(self.crops_groups):
            # crops = list(np.array(self.crops_acea)[np.array(self.crops_groups) == group])
            
            df_crops = final_df[final_df['Group']==group]
            crops = df_crops.index.to_list()
            _wfu = wf_prod[crops].tail(avg_year).sum().sum()/prod[crops].tail(avg_year).sum().sum()
            _wfu_change = _wfu/(wf_prod[crops].head(avg_year).sum().sum()/prod[crops].head(avg_year).sum().sum())-1
            _wft_change = wf_prod[crops].tail(avg_year).sum().sum()/wf_prod[crops].head(avg_year).sum().sum()-1
            
            _lfu = lf_prod[crops].tail(avg_year).sum().sum()/prod[crops].tail(avg_year).sum().sum()
            _lfu_change = _lfu/(lf_prod[crops].head(avg_year).sum().sum()/prod[crops].head(avg_year).sum().sum())-1
            _lft_change = lf_prod[crops].tail(avg_year).sum().sum()/lf_prod[crops].head(avg_year).sum().sum()-1
            
            _ha_rf_change = ha_rf[crops].tail(avg_year).sum().sum()/ha_rf[crops].head(avg_year).sum().sum()-1
            _ha_ir_change = ha_ir[crops].tail(avg_year).sum().sum()/ha_ir[crops].head(avg_year).sum().sum()-1
            
            _y_rf_change = ((y_rf[crops].tail(avg_year)*ha_rf[crops].tail(avg_year)).sum().sum()/ha_rf[crops].tail(avg_year).sum().sum()) \
                            /((y_rf[crops].head(avg_year)*ha_rf[crops].head(avg_year)).sum().sum()/ha_rf[crops].head(avg_year).sum().sum())-1
            _y_ir_change = ((y_ir[crops].tail(avg_year)*ha_ir[crops].tail(avg_year)).sum().sum()/ha_ir[crops].tail(avg_year).sum().sum()) \
                            /((y_ir[crops].head(avg_year)*ha_ir[crops].head(avg_year)).sum().sum()/ha_ir[crops].head(avg_year).sum().sum())-1
            
            crop_group = [group, '-', df_crops['% of total WF'].sum(), df_crops['% of total LF'].sum(),
                          ha_ir[crops].tail(avg_year).sum().sum()/(ha_rf[crops]+ha_ir[crops]).tail(avg_year).sum().sum(),
                          prod[crops].tail(avg_year).sum().sum()/(ha_rf[crops]+ha_ir[crops]).tail(avg_year).sum().sum(),
                          _wfu,
                          (WFu_g[crops].tail(avg_year)*prod[crops].tail(avg_year)).sum().sum()/wf_prod[crops].tail(avg_year).sum().sum(),
                          1-(WFu_g[crops].tail(avg_year)*prod[crops].tail(avg_year)).sum().sum()/wf_prod[crops].tail(avg_year).sum().sum(),
                          _wfu_change,_wft_change,
                          _lfu_change,_lft_change,
                          _ha_rf_change,_ha_ir_change,
                          _y_rf_change,_y_ir_change
                          ]
            final_df= final_df.append(pd.Series(crop_group, index=final_df.columns, name=group[:3].lower()), ignore_index=False)
            
        # World average
        world = ['World', '-','1','1', world_ir_lf, '-', '-', world_wfg_pct, 1-world_wfg_pct, world_wfu_change, world_wft_change,
                 world_lfu_change, world_lft_change, world_ha_rf_change, world_ha_ir_change, '-', '-']
        final_df= final_df.append(pd.Series(world, index=final_df.columns, name='wrl'), ignore_index=False)
        final_df.to_csv(ac.GetPath([folder, f'GlobalAnalysisChange_{self.text_years}.txt']), sep='\t')
    
    
    def GlobalTotals(self):
        "Calculate the main numbers"
        
        # Extract data
        folder = ac.GetPath(["outputs", 'main', 'summary'])
        _crop_report = ac.GetPath([folder, f'Global_stats_all_crops_All_{self.text_years}.txt'])
        prod = pd.read_csv(_crop_report, skiprows=3, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)
        crop_groups = pd.read_csv(_crop_report, skiprows=1, sep='\t', nrows=1, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)
        
        # WF [m3 t-1]
        WFu_g = pd.read_csv(_crop_report, skiprows=65, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        WFu_b_cr = pd.read_csv(_crop_report, skiprows=96, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        WFu_b_ir = pd.read_csv(_crop_report, skiprows=127, sep='\t', nrows=30, encoding='latin1',index_col=0).replace('--', 0.).replace(np.nan, 0.)
        
        # HA [m2]
        ha_rf =  pd.read_csv(_crop_report, skiprows=313, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)*10**4
        ha_ir = pd.read_csv(_crop_report, skiprows=344, sep='\t', nrows=30, index_col=0, encoding='latin1').replace('--', 0.).replace(np.nan, 0.)*10**4
        
        # Final dataframe
        # col_names = ['Year', 'WF total [Bm3]', 'WF green','WF blue ir', 'WF blue cr', 'Harvested area [Bm2]', 'Harvested area rf', 'Harvested area ir']
        years = list(range(1990,2020))
        final_df = pd.DataFrame()                       
        final_df['Year'] = years
        final_df['WF green'] = np.sum((WFu_g*prod).to_numpy(), axis=(1)) / 10**9
        final_df['WF blue cr'] = np.sum((WFu_b_cr*prod).to_numpy(), axis=(1))/ 10**9
        final_df['WF blue ir'] = np.sum((WFu_b_ir*prod).to_numpy(), axis=(1))/ 10**9
        final_df['WF total [Bm3]'] = final_df['WF green'] + final_df['WF blue cr'] + final_df['WF blue ir']
        final_df['Harvested area rf'] = np.sum(ha_rf.to_numpy(), axis=(1)) / 10**9
        final_df['Harvested area ir'] = np.sum(ha_ir.to_numpy(), axis=(1)) / 10**9
        final_df['Harvested area [Bm2]'] = final_df['Harvested area rf'] + final_df['Harvested area ir']
        
        for c_group in sorted(set(crop_groups.values.tolist()[0])):
            crop_list = np.array(self.crops_acea)[np.array(self.crops_groups)==c_group]
            final_df[f'{c_group} harvested area'] = (np.sum(ha_rf[crop_list], axis=(1)) + np.sum(ha_ir[crop_list], axis=(1))).to_numpy() / 10**9 
            final_df[f'{c_group} WF total'] = (np.sum((WFu_g*prod)[crop_list], axis=(1)) + np.sum((WFu_b_cr*prod)[crop_list], axis=(1)) + np.sum((WFu_b_ir*prod)[crop_list], axis=(1))).to_numpy() / 10**9 
        
        final_df.to_csv(ac.GetPath([folder, f'Global_totals_{self.text_years}.txt']), index=False, sep='\t')
        print(final_df)
    
    def NationalReport_AllCrops_And_Countries(self):
        "Create a table with avg hist national WF"
        
        folder_path = ac.GetPath(["outputs", 'main', 'national_wf_stats_hist_borders'])
        crops = fc.GetSimulatedCrops()
        cols = ['crop_code', 'country_code', 'year', 'harvarea_ha','irrigated_fraction', 'production_t', 'crop_yield_t_ha', 'wfg_m3t', 'wfb_cr_m3t', 'wfb_i_m3t', 'wf_tot_m3t']
        final_table = pd.DataFrame(columns=cols)
        
        crop_count = 1
        for crop_code in crops['FAO']:
            print(f'Crop {crop_count} out of {len(crops)}'); crop_count+=1
            _report_path = ac.GetPath([folder_path, f'{crop_code}_WF_HistNations_annual_1990_2019.csv'])
            np_data = pd.read_csv(_report_path, index_col=0).replace('--', 0.).replace(np.nan, 0.).to_numpy()
            for r in range(len(np_data)):
                final_table.loc[len(final_table)] = [crop_code, np_data[r,0], np_data[r,1], int(np_data[r,7]),
                                                     np_data[r,8], int(np_data[r,6]), np_data[r,10],
                                                     np_data[r,2], np_data[r,3], np_data[r,4], np_data[r,5]]
    
        # Save the final table
        out_report = ac.GetPath(["outputs", 'main', 'summary', f'AllCrops_WF_HistNations_annual_{self.text_years}.csv'])
        if os.path.exists(out_report): os.remove(out_report)
        final_table.to_csv(out_report, index=False)
    
    def GlobalReport_AllCrops(self):
        "Txt annual report of global averages"
        
        # 1. Create report file
        out_folder = ac.GetPath(["outputs", 'main'])
        if not os.path.exists(out_folder): os.makedirs(out_folder)
        out_report = ac.GetPath([out_folder,'summary', f'Global_stats_all_crops_{self.title}_{self.text_years}.txt'])
        if os.path.exists(out_report): os.remove(out_report)
        
        # 2. Collect data
        WF_info = np.zeros((self.crop_num, self.gs_num, 4)); LF = np.zeros((self.crop_num, self.gs_num, 3));
        ir_info =  np.zeros((self.crop_num, self.gs_num, 3)); crop_info = np.zeros((self.crop_num, self.gs_num, 7))
        c_count = 0
        for c in self.crops_acea:
            if os.path.exists(ac.GetPath(["outputs",f"{c}_drv_phd_mialyk_2022"])): project = f"{c}_drv_phd_mialyk_2022"
            elif os.path.exists(ac.GetPath(["outputs",f"{c}_phd_mialyk_2022"])): project = f"{c}_phd_mialyk_2022"
            elif os.path.exists(ac.GetPath(["outputs",f"{c}_phd_mialyk"])): project = f"{c}_phd_mialyk"
            elif os.path.exists(ac.GetPath(["outputs",f"{c}_drv_phd_mialyk"])): project = f"{c}_drv_phd_mialyk"
            else:
                raise Exception("Project is not found") 

            crop_fao = getattr(importlib.import_module(f'projects.{project}'), 'project_conf').crop_fao# Import project data
                
            _raster_path = ac.GetPath([out_folder,'global_wf_stats', f'{crop_fao}_WF_GlobalStats_{self.text_years}.txt'])
            
            # Get irrigation info
            _data = pd.read_csv(_raster_path, skiprows=33, sep='\t', nrows=30).replace('--', 0.).replace(np.nan, 0.)
            ir_info[c_count,:,0] = np.around(_data['Irrigation[mm]'].values.astype(float),2)
            ir_info[c_count,:,1] = np.around(_data['CWU[mm]'].values.astype(float),2)  # irrigated CWU
            crop_info[c_count,:,2] = np.around(_data['Harvested area[ha]'].values.astype(float),2) # irrigated HA
            crop_info[c_count,:,5] = np.around(_data['Yield[t ha-1]'].values.astype(float),2) # irrigated yield
                
            # Get rainfed info
            _data = pd.read_csv(_raster_path, skiprows=1, sep='\t', nrows=30).replace('--', 0.).replace(np.nan, 0.)
            ir_info[c_count,:,2] = np.around(_data['CWU[mm]'].values.astype(float),2) # rainfed info to compare
            crop_info[c_count,:,3] = np.around(_data['Harvested area[ha]'].values.astype(float),2) # rainfed HA
            crop_info[c_count,:,4] = np.around(_data['Yield[t ha-1]'].values.astype(float),2) # rainfed yield
            
            # Get general info
            _data = pd.read_csv(_raster_path, skiprows=65, sep='\t', nrows=30).replace('--', 0.).replace(np.nan, 0.)
            crop_info[c_count,:,0] = np.around(_data['Production[t]'].values.astype(float),2)
            crop_info[c_count,:,1] = np.around(_data['Yield[t ha-1]'].values.astype(float),2)
            crop_info[c_count,:,6] = np.around(_data['Scaling factor'].values.astype(float),2)
            
            # Get WF info
            WF_info[c_count,:,0] = _data['Unit WF[m3 t-1]'].values.astype(float) # Unit WF tot (green+blue) [m3 t-1]
            WF_info[c_count,:,1] = np.around(_data['Green'].values.astype(float),2) # Unit WF green [m3 t-1]
            WF_info[c_count,:,2] = np.around(_data['Blue cr'].values.astype(float),2) # Unit WF blue cr [m3 t-1]
            WF_info[c_count,:,3] = np.around(_data['Blue ir'].values.astype(float),2) # Unit WF blue ir [m3 t-1]
            
            # Get total LF
            # _raster_path = ac.GetPath([out_folder, 'global_lf_stats', f'{crop_fao}_LF_GlobalStats_{self.text_years}.txt'])
            # _data = pd.read_csv(_raster_path, skiprows=65, sep='\t', nrows=30).replace('--', 0.).replace(np.nan, 0.)
            
            # LF[c_count,:,0] = _data['Unit LF[m2 t-1]'].values.astype(float)
            # LF[c_count,:,1] = np.around(_data['Production[t]'].values.astype(float) * LF[c_count,:,0] /10**9,2) # Production[Gm2]

            c_count +=1
        
        # 3. Write report
        ac.printLog(out_report, "Crop name\t" + "\t".join(self.crops_name))
        ac.printLog(out_report, "FAO code\t" + "\t".join(map(str, self.crops_fao)))
        ac.printLog(out_report, "FAO group\t" + "\t".join(map(str, self.crops_groups)))
        
        # Production
        ac.printLog(out_report, 'Production[t]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, crop_info[:,gs,0])))
        
        # Unit WF tot
        ac.printLog(out_report, 'Unit WF tot[m3 t-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, WF_info[:,gs,0])))
        
        # Unit WF green
        ac.printLog(out_report, 'Unit WF green[m3 t-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, WF_info[:,gs,1])))
        
        # Unit WF blue ir
        ac.printLog(out_report, 'Unit WF blue cr[m3 t-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, WF_info[:,gs,2])))
        
        # Unit WF blue ir
        ac.printLog(out_report, 'Unit WF blue ir[m3 t-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, WF_info[:,gs,3])))
        
        # Unit LF
        ac.printLog(out_report, 'Unit LF[m2 t-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,0])))
        
        # # Prod LF tot
        # ac.printLog(out_report, '\tprod LF[Gm2]'); ac.printLog(out_report, "Year\t" + "\t".join(self.crops_name))
        # for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, LF[:,gs,1])))
        
        # Yield 
        ac.printLog(out_report, 'Yield[t ha-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, crop_info[:,gs,1])))
        
        # CWU rf
        ac.printLog(out_report, 'CWU rf[mm]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, ir_info[:,gs,2])))
        
        # CWU ir
        ac.printLog(out_report, 'CWU ir[mm]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, ir_info[:,gs,1])))
        
        # Irrigation
        ac.printLog(out_report, 'Irrigation[mm]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, ir_info[:,gs,0])))
        
        # Harvested area rf
        ac.printLog(out_report, 'Harvested area rf[ha]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, crop_info[:,gs,3])))
        
        # Harvested area ir
        ac.printLog(out_report, 'Harvested area ir[ha]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, crop_info[:,gs,2])))
        
        # Yield rf
        ac.printLog(out_report, 'Yield rf[t ha-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, crop_info[:,gs,4])))
        
        # Yield ir
        ac.printLog(out_report, 'Yield ir[t ha-1]\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, crop_info[:,gs,5])))
        
        # Scaling factor
        ac.printLog(out_report, 'Scaling factor\t' + "\t".join(self.crops_acea))
        for gs in range(self.gs_num): ac.printLog(out_report, f"{self.clock_start_year+gs}\t" + "\t".join(map(str, crop_info[:,gs,6])))