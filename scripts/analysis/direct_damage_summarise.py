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
    # asset_data_details = asset_data_details[(asset_data_details["sector"] == "energy") & (asset_data_details["asset_layer"] == "edges")]
    # print (asset_data_details)

    param_values = pd.read_csv('parameter_combinations.txt', sep=" ")
    uncertainty_columns = param_values.columns.values.tolist()[1:]
    damage_results_types = ["direct_damages","EAD"]
    include_confidence_values = [True,False]
    for confidence in include_confidence_values:
        for asset_info in asset_data_details.itertuples():
            asset_damages_results = os.path.join(direct_damages_results,f"{asset_info.asset_gpkg}_{asset_info.asset_layer}")
            for damages_type in damage_results_types:
                damage_files = [os.path.join(
                                    asset_damages_results,
                                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_{damages_type}_parameter_set_{param.parameter_set}.csv"
                                    ) for param in param_values.itertuples()]
                damage_results = [pd.read_csv(file) for file in damage_files if os.path.isfile(file) is True]
                print ("* Done with creating list of all dataframes")

                if damage_results:
                    haz_rcp_epochs = list(set(damage_results[0].set_index(["hazard","rcp","epoch"]).index.values.tolist()))
                    print (haz_rcp_epochs)
                    summarised_damages = []
                    for i,(haz,rcp,epoch) in enumerate(haz_rcp_epochs):
                        damages = [df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch)] for df in damage_results]
                        damages = pd.concat(damages,axis=0,ignore_index=True)
                        print ("* Done with concatinating all dataframes")
                        if confidence is False:
                            damages.drop("confidence",axis=1,inplace=True)
                            confidence_string = "_without_confidence_value"
                        else:
                            confidence_string = ""
                        index_columns = [c for c in damages.columns.values.tolist() if (
                                                    c not in ["direct_damage_cost","direct_reopen_cost"]
                                                    ) and ("EAD_" not in c) and ("EAR_" not in c)]
                        index_columns = [i for i in index_columns if i not in uncertainty_columns]
                        damage_columns = [c for c in damages.columns.values.tolist() if (
                                                    c in ["direct_damage_cost","direct_reopen_cost"]
                                                    ) or ("EAD_" in c) or ("EAR_" in c)]
                    
                        if len(damages.index) > 0:
                            summarised_damages.append(quantiles(damages,index_columns,damage_columns))
                        print (f"* Done with {haz} {rcp} {epoch}")
                    summarised_damages = pd.concat(summarised_damages,axis=0,ignore_index=True)
                    print ("* Done with concatinating results")
                    summarised_damages.to_csv(os.path.join(summary_results,
                                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_{damages_type}{confidence_string}.csv"),index=False)

                print (f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer} {damages_type} with confidence value {confidence}")

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)