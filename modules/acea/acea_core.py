# -*- coding: utf-8 -*-
"""
Purpose: General functions for ACEA processing scripts

Created: Thu Feb 25 11:03:11 2021

Author: Oleksandr Mialyk (o.mialyk@utwente.nl) at University of Twente (the Netherlands)
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import os
from osgeo import gdal

#%% General
def printLog(file, text):
    print(text)
    with open(file,'a') as file:
        print(text, file=file)
        
def GetPath(items):
    
    if type(items) == list:
        items = os.path.join(*map(str, items))
    if items.find(os.getcwd()) > 0:
        return items
    else:
        return os.path.join(os.getcwd(),items)

def ReadNC(file, var='3-dimensional Scientific Dataset'):
    "Read a NETCDF4 raster"
   
    data = Dataset(file, "r", format="NETCDF4")
    try: data = np.ma.squeeze(data[f"/{var}"][:,:,:]) # check if 3D
    except:
        try: data = np.ma.squeeze(data[f"/{var}"][:,:]) # check if 2D
        except: data = np.ma.squeeze(data[f"/{var}"][:]) # otherwise 1D
    return data

def GetResolution(resolution=0):
    if resolution == 1: # 5 arc minutes
        incrm = round(0.5/6,6)
        lims = [2160, 4320]
        lats = np.arange(90-incrm/2, -lims[0]/24, -incrm)
        lons = np.arange(-180+incrm/2, lims[1]/24, incrm)
    else: # 30 arc minutes
        incrm = 0.5
        lims = [360, 720]
        lats = np.arange(90-incrm/2, -lims[0]/4, -incrm)
        lons = np.arange(-180+incrm/2, lims[1]/4, incrm)
        
    return lats, lons, lims

def GetCountryCode(country_name):
    "Get countries"
    
    with open(GetPath(["data", "acea",'grid','regions_countries.txt']), 'rb') as f:
        file = pd.read_csv(f, sep='\t', lineterminator='\r')
    
    try: country_code = int(file[file['ISO3'] == str(country_name)]['FAOSTAT'])
    except: country_code = False; raise ValueError(f'No such country as {country_name}')
        
    return country_code

def GetAllCells30():
    "Select global land cells at 30 arc minutes"

    gridcells30 = ReadNC(GetPath(["data", "acea",'grid',"id30.nc"]),'id30')
    row30_list, col30_list = np.ma.where(gridcells30>=0); id30_list = gridcells30[row30_list,col30_list]
    gridcells30 = pd.DataFrame(np.array([id30_list,row30_list,col30_list]).astype(int).transpose(), columns=['id30', 'rowy30', 'colx30']).drop_duplicates().astype(int)
    
    return gridcells30

def GetAllCells5():
    "Select global land cells at 5 arc minutes"
    
    gridcells5 = ReadNC(GetPath(["data", "acea",'grid',"id5.nc"]))
    row5_list, col5_list = np.where(gridcells5>=0); id5_list = gridcells5[row5_list, col5_list]
    gridcells5 = pd.DataFrame(np.array([id5_list,row5_list, col5_list]).astype(int).transpose(), columns=['id5', 'rowy5', 'colx5']).drop_duplicates().astype(int)
    
    return gridcells5

def GetSimulationCells30(data):
    "Select grid cells for simulation at 30 arc minutes"
    
    gridcells30 = []
    if type(data) == str:
        if data == 'world': # All cropland of a crop
            gridcells30 = GetAllCells30()
            
            return gridcells30['id30'].tolist()
        else: # Country code
            gridcells30 = GetCountryCells30(data)['id30'].tolist()
    elif type(data) == list:
        if all(isinstance(n, int) for n in data):
            gridcells30 = data[:]
        else:
            for country in data:
                temp = GetCountryCells30(country)
                gridcells30 = gridcells30 + temp['id30'].tolist()
    else:
        raise ValueError('Wrong configuration for gridcells was passed')
        
    return gridcells30

def GetHarvestedAreas5(crop_code, irrigated, source = 'spam2010'):
    "Select harvested areas for simulation at 5 arc minutes"
    
    if source == 'spam2010': #SPAM2010 dataset https://doi.org/10.7910/DVN/PRFF8V
        try:
            if irrigated: dataset = gdal.Open(GetPath(["data", "acea",'harvested_areas',crop_code,f'spam2010V2r0_global_H_{crop_code}_I.tif']))
            else: dataset = gdal.Open(GetPath(["data", "acea",'harvested_areas',crop_code,f'spam2010V2r0_global_H_{crop_code}_R.tif']))
        except: raise ValueError(f'Harvested area files for {crop_code} were not found')
        data = dataset.GetRasterBand(1).ReadAsArray()
        return np.ma.masked_where(data <= 0, data)
    elif source == 'mirca2000':
        try:
            if irrigated: dataset = ReadNC(GetPath(["data", "acea",'harvested_areas',crop_code,f'{crop_code}_ir_MIRCA_data.nc']), 'value')[0,:,:]
            else: dataset = ReadNC(GetPath(["data", "acea",'harvested_areas',crop_code,f'{crop_code}_rf_MIRCA_data.nc']), 'value')[0,:,:]
        except: raise ValueError(f'Harvested area files for {crop_code} were not found')
        return dataset
    elif source == 'gaez2015':
        try:
            if irrigated: dataset =  gdal.Open(GetPath(["data", "acea",'harvested_areas',crop_code,f'GAEZAct2015_HarvArea_Irrigated_{crop_code}.tif']))
            else: dataset =  gdal.Open(GetPath(["data", "acea",'harvested_areas',crop_code,f'GAEZAct2015_HarvArea_Rainfed_{crop_code}.tif']))
        except: raise ValueError(f'Harvested area files for {crop_code} were not found')
        data = dataset.GetRasterBand(1).ReadAsArray()
        return np.ma.masked_where(data <= 0, data)
    else:
        raise ValueError(f'Harvested area files for {crop_code} were not found')
        return

# def GetCountryCells5(country_code):
#     "Select grid cells of a country at 5 arc minutes"
    
#     gridcells5 = np.array([ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'id5').flatten(), \
#                        ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'rowy5').flatten(), \
#                        ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'colx5').flatten(), \
#                        ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'faostat').flatten()])
    
#     pd_gridcells5 = pd.DataFrame(gridcells5.transpose(), columns=['id5', 'rowy5', 'colx5', 'faostat'])
#     country_cells = pd_gridcells5.where(pd_gridcells5['faostat']==country_code).dropna()[['id5', 'rowy5', 'colx5']].drop_duplicates().astype(int)
    
#     if len(country_cells) == 0: # For no longer existing countries (e.g. USSR)
#         hist_id5 = np.load(GetPath(["data", "acea",'grid','faostat_country_id5.npz']), allow_pickle=True)['id5'] # Read the output file
#         cells_id5 = hist_id5[hist_id5[:,1]==country_code][:,0]
#         _country_cells = gridcells5[:,np.isin(gridcells5[0,:], cells_id5)]
#         country_cells = pd.DataFrame(_country_cells[:-1,:].transpose(), columns=['id5', 'rowy5', 'colx5']).drop_duplicates().astype(int)
    
#     return country_cells

def GetCountryCells30(country_name):
    "Select grid cells of a country at 30 arc minutes"
    
    faostat_code = GetCountryCode(country_name)
    
    # For current countries
    m49_countries = ReadNC(GetPath(["data", "acea",'grid',"countries_m49_30arc.nc"]))
    id30 = ReadNC(GetPath(["data", "acea",'grid',"id30.nc"]),'id30')
    
    with open(GetPath(["data", "footprints",'regions',"FAOSTAT_hist_countries.txt"]), 'rb') as f:
        _country_data = pd.read_csv(f, sep='\t', lineterminator='\r')  
    m49_code = _country_data.loc[_country_data['Country Code'] == faostat_code]['M49 Code'].values[0]
    row30_list, col30_list = np.where(m49_countries==m49_code); id30_list = id30[row30_list,col30_list]
    country_cells = pd.DataFrame(np.array([id30_list,row30_list,col30_list]).astype(int).transpose(), columns=['id30', 'rowy30', 'colx30'])
    
    return country_cells
# def GetRegionList():
#     "Get world regions"
    
#     with open(GetPath(["data", "acea",'grid','regions_countries.txt']), 'rb') as f:
#         country_data = pd.read_csv(f, sep='\t', lineterminator='\r')
    
#     regional_codes = country_data[['Region Code','Region Name']].drop_duplicates().dropna().sort_values(by=['Region Name'])
    
#     return regional_codes.to_numpy()

# def GetSubRegionListByID(region_id):
#     "Get subregions by region id"
    
#     with open(GetPath(["data", "acea",'grid','regions_countries.txt']), 'rb') as f:
#         country_data = pd.read_csv(f, sep='\t', lineterminator='\r')
    
#     subregional_codes = country_data.where(country_data['Region Code']==region_id).dropna(how='all')[['Sub-region Code', 'Intermediate Region Code', 'Sub-region Name', 'Intermediate Region Name']].drop_duplicates().sort_values(by=['Sub-region Name', 'Intermediate Region Name'])

#     return subregional_codes.to_numpy()

# def GetSubRegionCells5(sub_region_id, sub_sub_region_id):
#     "Select global cells of a sub-region"
    
#     with open(GetPath(["data", "acea",'grid','regions_countries.txt']), 'rb') as f:
#         country_data = pd.read_csv(f, sep='\t', lineterminator='\r')
        
#     gridcells5 = np.array([ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'id5').flatten(), \
#                        ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'rowy5').flatten(), \
#                        ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'colx5').flatten(), \
#                        ReadNC(GetPath(["data", 'acea','grid','current_countries_5arc.nc']), 'faostat').flatten()])
    
#     pd_gridcells5 = pd.DataFrame(gridcells5.transpose(), columns=['id5', 'rowy5', 'colx5', 'faostat'])
    
    
#     if np.isnan(sub_sub_region_id):countries = country_data.where(country_data['Sub-region Code']==sub_region_id).dropna(how='all')
#     else:countries = country_data.where(country_data['Intermediate Region Code']==sub_sub_region_id).dropna(how='all')
    
#     # List of cells
#     pd_gridcells5 = pd_gridcells5[pd_gridcells5['faostat'].isin(countries['FAOSTAT'].to_numpy())].drop_duplicates()
    
#     # A raster mask
#     cols = 4320; pd_gridcells5_mask = np.zeros((int(cols/2), cols))
#     _pd_gridcells5 = pd_gridcells5.to_numpy()
#     for j in range(len(_pd_gridcells5)):
#         row5 = _pd_gridcells5[j][1]; col5 = _pd_gridcells5[j][2]
#         pd_gridcells5_mask[int(row5), int(col5)] = 1
    
#     return pd_gridcells5, np.ma.masked_where(pd_gridcells5_mask == 0, pd_gridcells5_mask)

def GetListOfVaribales():
    'Get the dictionary with variables calculated by ACEA'
    
    values = {'yield': ['Dry crop yields', 't ha-1 yr-1 (dry matter)'],
            'biom': ['Dry biomass', 't ha-1 yr-1 (dry matter)'],
            'pirnreq': ['Irrigation demand', 'mm'],
            'aet': ['Evapotranspiration', 'mm'],
            'soilevap': ['Evaporation', 'mm'],
            'transp': ['Transpiration', 'mm'],
            'hi': ['Harvest index', 'fraction'],
            'cwu_green': ['CWU green', 'mm'],
            'cwu_blue_cr': ['CWU capillary rise', 'mm'],
            'cwu_blue_ir': ['CWU irrigation', 'mm'],
            'gdd': ['Growing degree days', 'GDDs'],

            'plantday': ['Planting day', 'day of year'],
            'anthday': ['Days until anthesis', 'days from planting'],
            'matyday': ['Days until maturity', 'days from planting'],
            'harvyear': ['Year of harvest', 'calendar year'],
            'plantyear': ['Year of planting', 'calendar year'],
            }
    
    return values
    
def GetListOfScenarios():
    'Get the dictionary with simultion scenarios in ACEA'
    
    values = {1: 'Rainfed, no shallow groundwater',
                2: 'Rainfed considering shallow groundwater',
                3: 'Irrigated with surface irrigation',
                4: 'Irrigated with sprinklers',
                5: 'Irrigated with drip',
                6: 'Irrigated with surface irrigation considering shallow groundwater',
            }

    return values
def GetFAOtoFreshFactor(crop_code):
    "Get official FAO dry to fresh convertion factors by a FAOSTAT crop code"
    
    path = GetPath(["data", "acea",'crop_parameters', 'crop_dry_to_fresh_ratios.npz'])
    data = np.load(path, allow_pickle=True)['crop_ratios'] # Read the output file
    
    return data[data[:,0]==crop_code][0][1]

def GetIrrSystems5(consider_all=True):
    "Get the fraction of irrigation systems per irrigated cell in 5 acr minute arrays"
    
    ir_syst_data = ReadNC(GetPath(["data", 'acea', 'harvested_areas', 'cell_share_irrigation_systems_5arcmin.nc']), '3-dimensional Scientific Dataset')[:,:,:3] # Surface and Sprinkler fractions. Reference: https://doi.org/10.5194/hess-19-3073-2015
    ir_syst_data[ir_syst_data<0] = 0
    
    if consider_all:
        ir_syst_data[np.sum(ir_syst_data,2)==0] = [1, 0, 0]; sys_sum = np.sum(ir_syst_data,2)
        ir_syst_data_sur = ir_syst_data[:,:,0]/sys_sum
        ir_syst_data_spr = ir_syst_data[:,:,1]/sys_sum
        ir_syst_data_drp = ir_syst_data[:,:,2]/sys_sum
            	# Convert to system fraction, assume surface if nothing is reported
    else: # Only surface
        ir_syst_data_sur = np.ones((ir_syst_data.shape[0],ir_syst_data.shape[1]))
        ir_syst_data_spr = ir_syst_data_sur*0
        ir_syst_data_drp = ir_syst_data_sur*0
    
    return ir_syst_data_sur, ir_syst_data_spr, ir_syst_data_drp

#%% Rasters
def CreateHistRaster5(values, path, val_name, val_unit, title, start, max_gs):
    "Create 5 arc minute raster with many bands (years)"
    
    cols5 = 4320; rows5 = int(cols5/2) # Latitude and Longitude for 5 arcmin grid
    
    # Create a dataset
    raster = Dataset(path,'w', format='NETCDF4') #'w' stands for write
    raster.title = f'ACEA outputs for {title}'
    raster.institution = 'University of Twente, Netherlands'; raster.contact = 'Han Su (han_su20@163.com), Mialyk Oleksandr (o.mialyk@utwente.nl)'
    
    # Dimensions
    raster.createDimension('lat', rows5); raster.createDimension('lon', cols5); raster.createDimension('time', max_gs)
    
    # Variables
    lon = raster.createVariable('lon', 'f8', 'lon'); lat = raster.createVariable('lat', 'f8', 'lat'); time = raster.createVariable('time', 'f8', 'time')
    val = raster.createVariable(val_name, 'f4', ('time','lat','lon',), \
                                   zlib = True, complevel = 6, fill_value = 1.e20)
    # Adding values
    incrm = round(0.5/6,6)
    lat[:] = np.arange(90-incrm/2, -rows5/24, -incrm)
    lat.long_name = 'lat'; lat.units = 'degrees_north'
    lon[:] = np.arange(-180+incrm/2, cols5/24, incrm)
    lon.long_name = 'lon'; lon.units = 'degrees_east'
    time[:] = np.arange(start, start+max_gs, 1)
    time.long_name = 'time'; time.units = 'years'; time.calendar = 'standard'
    val[:,:,:] = values; val.units = val_unit
    raster.close()
    
def CreateRaster5(values, path, val_name, val_unit, title):
    "Create 5 arc minute raster with many bands (years)"
    
    cols5 = 4320; rows5 = int(cols5/2) # Latitude and Longitude for 5 arcmin grid
    
    # Create a dataset
    raster = Dataset(path,'w', format='NETCDF4') #'w' stands for write
    raster.title = f'ACEA outputs for {title}'
    raster.institution = 'University of Twente, Netherlands'; raster.contact = 'Han Su (han_su20@163.com), Mialyk Oleksandr (o.mialyk@utwente.nl)'
    
    # Dimensions
    raster.createDimension('lat', rows5); raster.createDimension('lon', cols5)
    
    # Variables
    lon = raster.createVariable('lon', 'f8', 'lon'); lat = raster.createVariable('lat', 'f8', 'lat')
    val = raster.createVariable(val_name, 'f4', ('lat','lon',), \
                                   zlib = True, complevel = 6, fill_value = 1.e20)
    # Adding values
    incrm = round(0.5/6,6)
    lat[:] = np.arange(90-incrm/2, -rows5/24, -incrm)
    lat.long_name = 'lat'; lat.units = 'degrees_north'
    lon[:] = np.arange(-180+incrm/2, cols5/24, incrm)
    lon.long_name = 'lon'; lon.units = 'degrees_east'
    val[:,:] = values; val.units = val_unit
    raster.close()
    
def CreateHistRaster30(values, path, val_name, val_unit, title, start, max_gs):
    "Create 30 arc minute raster with many bands (years)"
    
    cols30 = 720; rows30 = int(cols30/2) # Latitude and Longitude for 30 arcmin grid
    
    # Create a dataset
    raster = Dataset(path,'w', format='NETCDF4') #'w' stands for write
    raster.title = f'ACEA outputs for {title}'
    raster.institution = 'University of Twente, Netherlands'; raster.contact = 'Han Su (han_su20@163.com), Mialyk Oleksandr (o.mialyk@utwente.nl)'
    
    # Dimensions
    raster.createDimension('lat', rows30); raster.createDimension('lon', cols30); raster.createDimension('time', max_gs)
    
    # Variables
    lon = raster.createVariable('lon', 'f8', 'lon'); lat = raster.createVariable('lat', 'f8', 'lat'); time = raster.createVariable('time', 'f8', 'time')
    val = raster.createVariable(val_name, 'f4', ('time','lat','lon',), \
                                   zlib = True, complevel = 6, fill_value = 1.e20)
    # Adding values
    incrm = 0.5
    lat[:] = np.arange(90-incrm/2, -rows30/4, -incrm)
    lat.long_name = 'lat'; lat.units = 'degrees_north'
    lon[:] = np.arange(-180+incrm/2, cols30/4, incrm)
    lon.long_name = 'lon'; lon.units = 'degrees_east'
    time[:] = np.arange(start, start+max_gs, 1)
    time.long_name = 'time'; time.units = 'years'; time.calendar = 'standard'
    val[:,:,:] = values; val.units = val_unit
    raster.close()
    
def CreateRaster30(values, path, val_name, val_unit, title):
    "Create 5 arc minute raster with many bands (years)"
    
    cols30 = 720; rows30 = int(cols30/2) # Latitude and Longitude for 30 arcmin grid
    
    # Create a dataset
    raster = Dataset(path,'w', format='NETCDF4') #'w' stands for write
    raster.title = f'ACEA outputs for {title}'
    raster.institution = 'University of Twente, Netherlands'; raster.contact = 'Han Su (han_su20@163.com), Mialyk Oleksandr (o.mialyk@utwente.nl)'
    
    # Dimensions
    raster.createDimension('lat', rows30); raster.createDimension('lon', cols30)
    
    # Variables
    lon = raster.createVariable('lon', 'f8', 'lon'); lat = raster.createVariable('lat', 'f8', 'lat')
    val = raster.createVariable(val_name, 'f4', ('lat','lon',), \
                                   zlib = True, complevel = 6, fill_value = 1.e20)
    # Adding values
    incrm = 0.5
    lat[:] = np.arange(90-incrm/2, -rows30/4, -incrm)
    lat.long_name = 'lat'; lat.units = 'degrees_north'
    lon[:] = np.arange(-180+incrm/2, cols30/4, incrm)
    lon.long_name = 'lon'; lon.units = 'degrees_east'
    val[:,:] = values; val.units = val_unit
    raster.close()
    
def CreateLayeredRaster(values, layers, path, title):
    "Create raster with many layers"
    
    rows = values.shape[1]; cols = values.shape[2] # Latitude and Longitude
    
    # Create a dataset
    raster = Dataset(path,'w', format='NETCDF4') #'w' stands for write
    raster.title = 'AquaCrop-Earth@lternatives raster: %s' % title
    raster.institution = 'University of Twente, Netherlands'; raster.contact = 'Han Su (han_su20@163.com), o.mialyk@utwente.nl'
    
    # Dimensions
    raster.createDimension('lat', rows); raster.createDimension('lon', cols)
    
    # Variables
    lats = raster.createVariable('lat', 'f4', ('lat',)); lons = raster.createVariable('lon', 'f4', ('lon',))
    
    i = 0
    for layer_name in layers:
        val = raster.createVariable(layer_name, 'f4', ('lat','lon',), \
                                   zlib = True, complevel = 6, fill_value = -999)
        val[:,:] = values[i,:,:]; i += 1
    
    # Adding values
    if rows == 2160:
        incrm = round(0.5/6,6)
        lats[:] = np.arange(90-incrm/2, -rows/24, -incrm)
        lons[:] = np.arange(-180+incrm/2, cols/24, incrm)
    else:
        incrm = 0.5
        lats[:] = np.arange(90-incrm/2, -rows/4, -incrm)
        lons[:] = np.arange(-180+incrm/2, cols/4, incrm)
    
    lats.long_name = 'lat'; lats.units = 'degrees_north'
    lons.long_name = 'lon'; lons.units = 'degrees_east'
    
    raster.close()
    
# Outputs
def GetOutputsOfScenario(conf, scenario, scaling=True, cwu=True):
    "Get main ACEA outputs of a specific scenario"
    
    sim_yields, sim_cwu_green, sim_cwu_blue_ir, sim_cwu_blue_cr = [0,0,0,0]
    dry_to_fresh = GetFAOtoFreshFactor(conf.crop_fao)
    _,_,res5 = GetResolution(1)
    clock_start_year = int(conf.clock_start.split('/')[0]); clock_end_year = int(conf.clock_end.split('/')[0])
    max_gs = (clock_end_year - clock_start_year)-1; first_year = clock_end_year-max_gs+1
    raster_folder_path = GetPath(["outputs", f"{conf.project_name}",'rasters'])
    text_years = f'{first_year}_{clock_end_year}'
    
    # Get scaling factors for yields
    if scaling:
        scaling_factors = GetScalingFactors(conf.crop_fao, conf.crop_name_short, text_years)
    else:
        scaling_factors = np.ones((max_gs,res5[0],res5[1]))
    
    _raster_path = GetPath([raster_folder_path, f'acea_5arc_sc{scenario}_{conf.crop_name_short}_yield_global_annual_{text_years}.nc'])
    
    if os.path.exists(_raster_path): 
        sim_yields = ReadNC(_raster_path, f'yield-{conf.crop_name_short}') / dry_to_fresh * scaling_factors # [fresh t ha-1]
        
        if cwu: # Return CWU as well
            sim_cwu_green = ReadNC(GetPath([raster_folder_path, f'acea_5arc_sc{scenario}_{conf.crop_name_short}_cwu_green_global_annual_{text_years}.nc']), f'cwu_green-{conf.crop_name_short}')
            
            _raster_path = GetPath([raster_folder_path, f'acea_5arc_sc{scenario}_{conf.crop_name_short}_cwu_blue_cr_global_annual_{text_years}.nc'])
            if os.path.exists(_raster_path): sim_cwu_blue_cr = ReadNC(_raster_path, f'cwu_blue_cr-{conf.crop_name_short}')
    
                
            _raster_path = GetPath([raster_folder_path, f'acea_5arc_sc{scenario}_{conf.crop_name_short}_cwu_blue_ir_global_annual_{text_years}.nc'])
            if os.path.exists(_raster_path): sim_cwu_blue_ir = ReadNC(_raster_path, f'cwu_blue_ir-{conf.crop_name_short}')
        return sim_yields, sim_cwu_green, sim_cwu_blue_ir, sim_cwu_blue_cr
    else:
        print('No outputs were found'); return

def GetScalingFactors(crop_fao, crop_name_short, text_years):
    'Get scaling factors for yields with a limitation of 99.5%'
    
    _raster_path = f'scaling_factors_5arc_{crop_fao}_moving_avg_global_annual_{text_years}.nc'
    _raster_path = GetPath(["data","footprints",'scaling_factors', _raster_path])
    if os.path.exists(_raster_path): 
        scaling_factors = ReadNC(_raster_path, f'scaling_factor-{crop_name_short}')
        _threshold = np.percentile(scaling_factors[scaling_factors.mask==False],99) # Limit to the highest 1 % to avoid crazy values
        scaling_factors[(scaling_factors > _threshold) & (scaling_factors.mask==False)] = _threshold
    else: raise ValueError('Moving average of scaling factors is missing')
    
    return scaling_factors