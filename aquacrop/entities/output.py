import pandas as pd
import numpy as np


class Output:
    """
    Class to hold output data

    During Simulation these are numpy arrays and are converted to pandas dataframes
    at the end of the simulation

    Atributes:
    
        water_flux (pandas.DataFrame, numpy.array): Daily water flux changes

        water_storage (pandas.DataFrame, numpy array): daily water content of each soil compartment

        crop_growth (pandas.DataFrame, numpy array): daily crop growth variables

        final_stats (pandas.DataFrame, numpy array): final stats at end of each season

    """

    def __init__(self, time_span, initial_th):

        self.water_storage = np.zeros((len(time_span)-1, 3 + len(initial_th)))
        self.water_flux = np.zeros((len(time_span)-1, 16))
        self.crop_growth = np.zeros((len(time_span)-1, 20))
        
        self.ET_color = np.zeros((len(time_span)-1, 6))
        self.S_color = np.zeros((len(time_span)-1, 6))
        
        self.final_stats = pd.DataFrame(
            columns=[
                "Season",#0
                "crop Type",#1
                'Planting Date (YYYY/MM/DD)',#2
                'Anthesis Date (YYYY/MM/DD)',#3
                "Harvest Date (YYYY/MM/DD)",#4
                "Harvest Date (Step)",#5
                "Dry yield (tonne/ha)",#6
                "Fresh yield (tonne/ha)",#7
                "Yield potential (tonne/ha)",#8
                "Seasonal irrigation (mm)",#9
                'GDD cummulative',#10
                'Initial Soil Water (mm)',#11
                'Precipitation (mm)',#12
                'Capillary rise (mm)',#13
                'GW inflow (mm)',#14
                'Final Soil Water (mm)',#15
                'Evaporation (mm)',#16
                'Transpiration (mm)',#17
                'Runoff (mm)',#18
                'Percolation (mm)',#19
                'Biomass (dry tonne/ha)',#20
            ])
            
        self.final_watercolor = pd.DataFrame(
            columns=[
                "Season",#0
                'ET green (mm)',#1
                'ET blue irr (mm)',#2
                'ET blue cr (mm)',#3
                'Storage green (mm)',#4
                'Storage blue irr (mm)',#5
                'Storage blue cr (mm)'#6
            ]
        )