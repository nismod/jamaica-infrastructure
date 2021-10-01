"""Estimate direct damages to physical assets exposed to hazards

"""
import sys
import os

import pandas as pd
import geopandas as gpd
from shapely import wkb
import numpy as np
from SALib.sample import morris
import SALib.analyze.morris 

from analysis_utils import *
from tqdm import tqdm
tqdm.pandas()

def quantiles(dataframe,grouping_by_columns,grouped_columns):
    quantiles_list = ['mean','min','max','median','q5','q95']
    df_list = []
    for quant in quantiles_list:
        if quant == 'mean':
            # print (dataframe)
            df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].mean()
        elif quant == 'min':
            df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].min()
        elif quant == 'max':
            df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].max()
        elif quant == 'median':
            df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].quantile(0.5)
        elif quant == 'q5':
            df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].quantile(0.05)
        elif quant == 'q95':
            df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].quantile(0.95)

        df.rename(columns=dict((g,'{}_{}'.format(g,quant)) for g in grouped_columns),inplace=True)
        df_list.append(df)
    return pd.concat(df_list,axis=1).reset_index()

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
    epsg_jamaica = 3448

    direct_damages_results = os.path.join(output_data_path,"direct_damages")

    summary_results = os.path.join(output_data_path,"direct_damages_summary")
    if os.path.exists(summary_results) == False:
        os.mkdir(summary_results)

    asset_data_details = pd.read_csv(os.path.join(processed_data_path,
                        "networks",
                        "network_layers_hazard_intersections_details.csv"))

    param_values = pd.read_csv('parameter_combinations.txt', sep=" ")
    uncertainty_columns = param_values.columns.values.tolist()[1:]
    damage_results_types = ["direct_damages","EAD"]
    for asset_info in asset_data_details.itertuples():
        asset_damages_results = os.path.join(direct_damages_results,f"{asset_info.asset_gpkg}_{asset_info.asset_layer}")
        for damages in damage_results_types:
            damage_files = [os.path.join(
                                asset_damages_results,
                                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_{damages}_parameter_set_{param.parameter_set}.csv"
                                ) for param in param_values.itertuples()]
            damage_results = [pd.read_csv(file) for file in damage_files if os.path.isfile(file) is True]

            if damage_results:
                damage_results = pd.concat(damage_results,axis=0,ignore_index=True)
                index_columns = [c for c in damage_results.columns.values.tolist() if c != "direct_damage_cost" and "EAD" not in c]
                index_columns = [i for i in index_columns if i not in uncertainty_columns]
                damage_columns = [c for c in damage_results.columns.values.tolist() if c == "direct_damage_cost" or "EAD" in c]
                summarised_damages = quantiles(damage_results,index_columns,damage_columns)
                summarised_damages.to_csv(os.path.join(summary_results,
                            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_{damages}.csv"),index=False)

            print (f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer} {damages}")

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)