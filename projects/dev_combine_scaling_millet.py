# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 15:34:36 2022

@author: HanS, MialykO
"""
import modules.acea.acea_core as ac
import os
import numpy as np
import pandas as pd

project1 = r'outputs\PearlMillet\pmil_scaling.csv'
project2 = r'outputs\SmallMillet\smil_scaling.csv'

scaling1=pd.read_csv(project1, index_col='Unnamed: 0')
scaling2=pd.read_csv(project2, index_col='Unnamed: 0')

scaling1['mergedscaling']=0
scaling2['mergedscaling']=0

countries=set(list(scaling1.index)+list(scaling2.index))

for cty in countries:
    if cty in scaling1.index and cty in scaling2.index and\
    scaling1.loc[cty,'mergedscaling']==0 and scaling2.loc[cty,'mergedscaling']==0:
        area_spam=scaling1.loc[cty,'area_SPAM']+scaling2.loc[cty,'area_SPAM']
        yield_=(scaling1.loc[cty,'area_SPAM']*scaling1.loc[cty,'yield']+\
            scaling2.loc[cty,'area_SPAM']*scaling2.loc[cty,'yield'])/area_spam
        yieldf3=(scaling1.loc[cty,'area_SPAM']*scaling1.loc[cty,'yieldf3']+\
            scaling2.loc[cty,'area_SPAM']*scaling2.loc[cty,'yieldf3'])//area_spam
            
        scaling1.loc[cty,'scaling']=scaling1.loc[cty,'faoyield']/yield_
        scaling2.loc[cty,'scaling']=scaling2.loc[cty,'faoyield']/yield_
        
        scaling1.loc[cty,'scalingf3']=scaling1.loc[cty,'faoyield']/yieldf3
        scaling2.loc[cty,'scalingf3']=scaling2.loc[cty,'faoyield']/yieldf3
        
        scaling1.loc[cty,'scaling_area']=scaling1.loc[cty,'area_FAO']/area_spam
        scaling2.loc[cty,'scaling_area']=scaling2.loc[cty,'area_FAO']/area_spam
        
        scaling1.loc[cty,'mergedscaling']=1
        scaling2.loc[cty,'mergedscaling']=1
        

scaling1.to_csv(project1)
scaling2.to_csv(project2)
        
        
        
        
        