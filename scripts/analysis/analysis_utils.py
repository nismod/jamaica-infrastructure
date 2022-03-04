"""Functions for preprocessing data
"""
import sys
import os
import json

import pandas as pd
import geopandas as gpd
from scipy.interpolate import interp1d
from scipy import integrate
import math
import numpy as np
from tqdm import tqdm
tqdm.pandas()


def load_config():
    """Read config.json"""
    config_path = os.path.join(os.path.dirname(__file__), "..", "..", "config.json")
    with open(config_path, "r") as config_fh:
        config = json.load(config_fh)
    return config


def geopandas_read_file_type(file_path, file_layer, file_database=None):
    if file_database is not None:
        return gpd.read_file(os.path.join(file_path, file_database), layer=file_layer)
    else:
        return gpd.read_file(os.path.join(file_path, file_layer))

def curve_interpolation(x_curve,y_curve,x_new):
    if x_new <= x_curve[0]:
        return y_curve[0]
    elif x_new >= x_curve[-1]:
        return y_curve[-1]
    else:
        interpolate_values = interp1d(x_curve, y_curve)
        return interpolate_values(x_new)


def expected_risks_pivot(v,probabilites,probability_threshold,flood_protection_column):
    """Calculate expected risks
    """
    prob_risk = sorted([(p,getattr(v,str(p))) for p in probabilites],key=lambda x: x[0])
    if probability_threshold != 1:
        probability_threshold = getattr(v,flood_protection_column)
        if probability_threshold > 0:
            prob_risk = [pr for pr in prob_risk if pr[0] <= 1.0/probability_threshold]
    
    if len(prob_risk) > 1:
        risks = integrate.trapz(np.array([x[1] for x in prob_risk]), np.array([x[0] for x in prob_risk]))
    elif len(prob_risk) == 1:
        risks = 0.5*prob_risk[0][0]*prob_risk[0][1]
    else:
        risks = 0
    return risks

def risks_pivot(dataframe,index_columns,probability_column,
            risk_column,flood_protection_column,expected_risk_column,
            flood_protection=None,flood_protection_name=None):
    
    """
    Organise the dataframe to pivot with respect to index columns
    Find the expected risks
    """
    if flood_protection is None:
        # When there is no flood protection at all
        expected_risk_column = f"{expected_risk_column}_undefended"
        probability_threshold = 1 
    else:
        expected_risk_column = f"{expected_risk_column}_{flood_protection_name}"
        probability_threshold = 0 
        
    probabilites = list(set(dataframe[probability_column].values.tolist()))
    df = (dataframe.set_index(index_columns).pivot(
                                    columns=probability_column
                                    )[risk_column].reset_index().rename_axis(None, axis=1)).fillna(0)
    df.columns = df.columns.astype(str)
    df[expected_risk_column] = df.progress_apply(lambda x: expected_risks_pivot(x,probabilites,
                                                        probability_threshold,
                                                        flood_protection_column),axis=1)
    
    return df[index_columns + [expected_risk_column]]