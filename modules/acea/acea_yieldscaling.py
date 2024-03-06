# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:36:18 2022

@author: SuH
"""


import pandas as pd
from sqlalchemy import create_engine
engine = create_engine('postgresql://postgres:0220@localhost/postgres')
import xarray as xr


def read_spam(crop,farmingsystem):
    temp=pd.read_sql(f"select * from spam.spam2010v2r0_global_h_t{farmingsystem}",engine)
    temp=temp[['iso3','x','y',f'{crop}_{farmingsystem}']]
    temp=temp.loc[temp[f'{crop}_{farmingsystem}']>0,:]
    temp=temp.rename(columns={"x": "lon", "y": "lat",f'{crop}_{farmingsystem}':crop})
    
    return temp

def read_gaez(crop,farmingsystem):
    temp=pd.read_sql(f"select * from gaez.gaez_area_{farmingsystem}_location",engine)
    temp=temp[['country','x','y',f'{crop}_{farmingsystem}']]
    temp=temp.loc[temp[f'{crop}_{farmingsystem}']>0,:]
    temp=temp.rename(columns={"x": "lon", "y": "lat",f'{crop}_{farmingsystem}':crop,'country':'iso3'})
    
    return temp


def read_acea(filename,time=2010):
    temp=xr.open_dataset(filename)
    temp=temp.sel(time=time)
    
    return temp

irr_sys=xr.open_dataset(r'data\acea\harvested_areas\cell_share_irrigation_systems_5arcmin_ss.nc')

def combine_spam_acea(spam_,acea1,acea2,added_col,acea_col,\
                      irrigated=0,resol1=0.083334,resol2=0.500001,irr_sys=irr_sys,dry_fresh=1,acea3=[]):

    cols=spam_.columns
    
    if 'yield' in added_col:
        dryfresh=dry_fresh
    else:
        dryfresh=1
 
    if irrigated==1:
        #only two bands here

        tgt_x = xr.DataArray(list(spam_['lat']), dims="points")
        tgt_y = xr.DataArray(list(spam_['lon']), dims="points")
        temp=irr_sys.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol1)
        band1=temp['Band1'].values
        band2=temp['Band2'].values

        temp=acea1.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol2)
        acea_col1=temp[acea_col].values

        temp=acea2.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol2)
        acea_col2=temp[acea_col].values

        spam_['band1']=band1
        spam_['band2']=band2
        spam_['acea_col1']=acea_col1/dryfresh
        spam_['acea_col2']=acea_col2/dryfresh

        spam_=spam_.fillna(0)
        spam_.loc[(spam_['band1']+spam_['band2'])==0,'band1']=1
        sum_=(spam_['band1']+spam_['band2']).values
        spam_['band1']=spam_['band1']/sum_
        spam_['band2']=spam_['band2']/sum_

        spam_[added_col]=spam_['band1']*spam_['acea_col1']+spam_['band2']*spam_['acea_col2']
        
        spam_.loc[spam_['acea_col2']==0,added_col]=spam_['acea_col1']
        spam_.loc[spam_['acea_col1']==0,added_col]=spam_['acea_col2']
        
    elif irrigated==2:
        #only two bands here

        tgt_x = xr.DataArray(list(spam_['lat']), dims="points")
        tgt_y = xr.DataArray(list(spam_['lon']), dims="points")
        temp=irr_sys.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol1)
        band1=temp['Band1'].values
        band2=temp['Band2'].values
        band3=temp['Band3'].values

        temp=acea1.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol2)
        acea_col1=temp[acea_col].values

        temp=acea2.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol2)
        acea_col2=temp[acea_col].values
        
        temp=acea3.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol2)
        acea_col3=temp[acea_col].values

        spam_['band1']=band1
        spam_['band2']=band2
        spam_['band3']=band3
        spam_['acea_col1']=acea_col1/dryfresh
        spam_['acea_col2']=acea_col2/dryfresh
        spam_['acea_col3']=acea_col3/dryfresh

        spam_=spam_.fillna(0)
        spam_.loc[(spam_['band1']+spam_['band2']+spam_['band3'])==0,'band1']=1
        sum_=(spam_['band1']+spam_['band2']+spam_['band3']).values
        spam_['band1']=spam_['band1']/sum_
        spam_['band2']=spam_['band2']/sum_
        spam_['band3']=spam_['band3']/sum_

        spam_[added_col]=spam_['band1']*spam_['acea_col1']+spam_['band2']*spam_['acea_col2']+\
            spam_['band3']*spam_['acea_col3']
        
        #spam_.loc[spam_['acea_col2']==0,added_col]=spam_['acea_col1']
        #spam_.loc[spam_['acea_col1']==0,added_col]=spam_['acea_col2']
        
    elif irrigated==0:
        #acea2 with groundwater

        tgt_x = xr.DataArray(list(spam_['lat']), dims="points")
        tgt_y = xr.DataArray(list(spam_['lon']), dims="points")

        temp=acea1.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol2)
        acea_col1=temp[acea_col].values
        temp=acea2.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol2)
        acea_col2=temp[acea_col].values

        spam_['acea_col1']=acea_col1/dryfresh
        spam_['acea_col2']=acea_col2/dryfresh

        spam_=spam_.fillna(0)

        spam_[added_col]=spam_['acea_col1']
        
        spam_.loc[spam_['acea_col2']>0,added_col]=spam_['acea_col2']
        
    if irrigated==1 or irrigated==2:
        
        return spam_[list(cols)+[added_col]],spam_[['band1','band2']]

    else:
        return spam_[list(cols)+[added_col]]
    
def add_farm_size(farmsize_,t_):
    
    import pandas as pd
    import warnings
    
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
    columns=t_.columns
    resol=0.083333
    
    temp_=t_.set_index(['lat','lon']).to_xarray()

    tgt_x = xr.DataArray(list(farmsize_['lat_y']), dims="points")
    tgt_y = xr.DataArray(list(farmsize_['lon_x']), dims="points")

    temp=temp_.sel(lat=tgt_x,lon=tgt_y,method="nearest", tolerance=resol)
    
    for col in columns:
        if col in ['lon','lat']:
            pass
        else:
            col1=temp[col].values
            farmsize_[col]=col1
    
    farmsize_=farmsize_.drop(['Unnamed: 0'], axis=1)
    return farmsize_