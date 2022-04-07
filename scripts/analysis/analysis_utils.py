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
    interpolate_values = interp1d(x_curve, y_curve,fill_value=(min(y_curve),max(y_curve)),bounds_error=False)
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

def risks(dataframe,index_columns,probabilities,
            expected_risk_column,
            flood_protection_period=0,flood_protection_name=None):
    
    """
    Organise the dataframe to pivot with respect to index columns
    Find the expected risks
    """
    if flood_protection_name is None and flood_protection_period == 0:
        # When there is no flood protection at all
        expected_risk_column = f"{expected_risk_column}_undefended"
        probability_columns = [str(p) for p in probabilities]
        
    elif flood_protection_period > 0:
        if flood_protection_name is None:
            expected_risk_column = f"{expected_risk_column}_{flood_protection_period}_year_protection"
        else:
            expected_risk_column = f"{expected_risk_column}_{flood_protection_name}"
        
        probabilities = [pr for pr in probabilities if pr <= 1.0/flood_protection_period]
        probability_columns = [str(p) for p in probabilities]
    else:
        # When there is no flood protection at all
        expected_risk_column = f"{expected_risk_column}_{flood_protection_name}"
        probability_columns = [str(p) for p in probabilities]
        
    dataframe.columns = dataframe.columns.astype(str)
    dataframe[expected_risk_column] = list(integrate.trapz(dataframe[probability_columns].to_numpy(),
                                            np.array([probabilities*len(dataframe.index)]).reshape(dataframe[probability_columns].shape)))
    
    return dataframe[index_columns + [expected_risk_column]]

def calculate_discounting_arrays(discount_rate=4.5, growth_rate=5.0,
                                start_year=2020,end_year=2050,
                                maintain_period=4):
    """Set discount rates for yearly and period maintenance costs

    Parameters
    ----------
    discount_rate
        yearly discount rate
    growth_rate
        yearly growth rate

    Returns
    -------
    discount_rate_norm
        discount rates to be used for the costs
    discount_rate_growth
        discount rates to be used for the losses
    min_main_dr
        discount rates for 4-year periodic maintenance
    max_main_dr
        discount rates for 8-year periodic maintenance

    """
    discount_rates = []
    growth_rates = []

    for year in range(start_year,end_year+1):
        discount_rates.append(
            1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year))

    if isinstance(growth_rate, float):
        for year in range(start_year,end_year+1):
            growth_rates.append(
                1.0*math.pow(1.0 + 1.0*growth_rate/100.0, year - start_year))
    else:
        for i, (year,rate) in enumerate(growth_rate):
            if year > start_year:
                growth_rates.append(prod([1 + v[1]/100.0 for v in growth_rate[:i]]))
            else:
                growth_rates.append(1)  

    maintain_years = np.arange(start_year, end_year+1,maintain_period)
    maintain_rates = []
    for year in maintain_years[1:]:
        maintain_rates.append(1.0 / math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year))

    return np.array(discount_rates), np.array(growth_rates), np.array(maintain_rates)

# getting Product
def prod(val) :
    res = 1
    for ele in val:
        res *= ele
    return res

def calculate_growth_rate_factor(growth_rate=5.0,
                                start_year=2020,end_year=2050,
                        ):
    """Set discount rates for yearly and period maintenance costs

    Parameters
    ----------
    discount_rate
        yearly discount rate
    growth_rate
        yearly growth rate

    Returns
    -------
    discount_rate_norm
        discount rates to be used for the costs
    discount_rate_growth
        discount rates to be used for the losses
    min_main_dr
        discount rates for 4-year periodic maintenance
    max_main_dr
        discount rates for 8-year periodic maintenance

    """
    growth_rates = []

    if isinstance(growth_rate, float):
        for year in range(start_year,end_year+1):
            growth_rates.append(
                1.0*math.pow(1.0 + 1.0*growth_rate/100.0, year - start_year))
    else:
        for i, (year,rate) in enumerate(growth_rate):
            if year > start_year:
                growth_rates.append(prod([1 + v[1]/100.0 for v in growth_rate[:i]]))
            else:
                growth_rates.append(1)  

    return np.array(growth_rates)


def calculate_discounting_rate_factor(discount_rate=4.5,
                                start_year=2020,end_year=2050,
                                maintain_period=4,skip_year_one=False):
    """Set discount rates for yearly and period maintenance costs

    Parameters
    ----------
    discount_rate
        yearly discount rate
    growth_rate
        yearly growth rate

    Returns
    -------
    discount_rate_norm
        discount rates to be used for the costs
    discount_rate_growth
        discount rates to be used for the losses
    min_main_dr
        discount rates for 4-year periodic maintenance
    max_main_dr
        discount rates for 8-year periodic maintenance

    """
    discount_rates = []
    maintain_years = np.arange(start_year+1, end_year+1,maintain_period)
    for year in range(start_year,end_year+1):
        if year in maintain_years:
            discount_rates.append(
                1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year))
        else:
            if skip_year_one is True:
                discount_rates.append(0)
            else:
                discount_rates.append(1)

    return np.array(discount_rates)

def extract_growth_rate_info(growth_rates,time_column,rate_column,start_year=2000,end_year=2100):
    growth_year_rates = []
    growth_rates_times = list(sorted(growth_rates[time_column].values.tolist()))
    # And create parameter values
    for y in range(start_year,end_year+1):
        if y in growth_rates_times:
            growth_year_rates.append((y,growth_rates.loc[growth_rates[time_column] == y,rate_column].values[0]))
        elif y < growth_rates_times[0]:
            growth_year_rates.append((y,growth_rates.loc[growth_rates[time_column] == growth_rates_times[0],rate_column].values[0]))
        elif y > growth_rates_times[-1]:
            growth_year_rates.append((y,growth_rates.loc[growth_rates[time_column] == growth_rates_times[-1],rate_column].values[0])) 

    return growth_year_rates