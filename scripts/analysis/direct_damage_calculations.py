"""Estimate direct damages to physical assets exposed to hazards

"""
import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np

from analysis_utils import *
from tqdm import tqdm
tqdm.pandas()

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448

    
    damage_curve = pd.read_csv(os.path.join(incoming_data_path,
                            'damage_assessment_information',
                            'damage_functions',
                            'damage_curves_energy_flooding.csv'))
    damage_fraction = curve_interpolation(damage_curve['flood_depth_m'].values,
                                        damage_curve['substation'].values,
                                        0.25)

    damage_uncertainty = 1.0
    damage_fraction = damage_uncertainty*damage_fraction
    if damage_fraction >= 100.0:
        damage_fraction = 100.0

    print (damage_fraction)




if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)