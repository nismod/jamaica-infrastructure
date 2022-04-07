"""Generate hazard-damage curves
"""
import os
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
from collections import OrderedDict
from summarise_utils import *
from tqdm import tqdm
tqdm.pandas()

def convert_to_usd(x,loss_column):
    if ("$J" in str(x.damage_cost_unit)) or ("J$" in str(x.damage_cost_unit)) or ("JD" in str(x.damage_cost_unit)):
        return jamaica_currency_conversion()*x[loss_column]
    else:
        return x[loss_column]

def convert_to_jd(x,loss_column):
    if ("$US" in str(x.damage_cost_unit)) or ("US$" in str(x.damage_cost_unit)) or ("USD" in str(x.damage_cost_unit)):
        return (1.0/jamaica_currency_conversion())*x[loss_column]
    else:
        return x[loss_column]

def get_electricity_costs(electricity_data_path,
                            electricity_damage_df,
                            electricity_id_column,
                            electricity_damage_column,
                            electricity_capacity_column,
                            layer_type
                        ):
    electricity_cap_df = gpd.read_file(os.path.join(electricity_data_path,"electricity_network_v1.1.gpkg"),layer=f"{layer_type}_v1.1")
    electricity_damage_df = pd.merge(electricity_damage_df,
                                electricity_cap_df[[electricity_id_column,electricity_capacity_column]],
                                how="left",on=[electricity_id_column]).fillna(0)

    electricity_damage_df[
                    electricity_damage_column
                    ] = electricity_damage_df[
                                                electricity_capacity_column
                                                ]*electricity_damage_df[
                                                            electricity_damage_column
                                                            ]
    return electricity_damage_df

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']

    asset_summary_results = os.path.join(output_data_path,
                            "asset_attributes_summarised")
    if os.path.exists(asset_summary_results) == False:
        os.mkdir(asset_summary_results)

    asset_groupby_column = "asset_grouping_column"
    asset_cost_unit_column = "asset_cost_unit_column"
    asset_group_columns = [
                    "asset_capacity_column",
                    "asset_min_cost_column",
                    "asset_max_cost_column",
                    "asset_min_reopen_cost_column",
                    "asset_max_reopen_cost_column",
                    "asset_population_column",
                    "asset_usage_column"
                ]
    asset_data_details = pd.read_csv(os.path.join(processed_data_path,
                        "networks",
                        "network_layers_details.csv"))

    output_excel = os.path.join(asset_summary_results,'asset_summary_results.xlsx')
    output_wrtr = pd.ExcelWriter(output_excel) 
    for asset_info in asset_data_details.itertuples():
        asset_group_indexes = getattr(asset_info,asset_groupby_column).split(",")
        for group_index in asset_group_indexes:
            group_indexes = [group_index] + getattr(asset_info,asset_cost_unit_column).split(",")
            asset_groups = OrderedDict()
            for gr in asset_group_columns:
                asset_groups[gr]= getattr(asset_info,gr)
            # print (asset_groups)
            asset_groups = dict([(k,gr) for k,gr in asset_groups.items() if gr != "none"])
            # print (asset_groups)
            asset_df = gpd.read_file(os.path.join(processed_data_path,asset_info.path),layer=asset_info.asset_layer).to_crs(epsg=3448)
            if asset_info.asset_gpkg == "rail" and asset_info.asset_layer == "nodes":
                asset_df = asset_df[asset_df["asset_type"] == "station"]
            asset_df["count"] = 1
            asset_sums = ["count"] 
            if asset_info.asset_layer == "edges" or asset_info.asset_layer == "edges_v1.1":
                asset_df["length"] = 0.001*asset_df.geometry.length
                asset_sums += ["length"] 
            elif asset_info.asset_layer == "areas":
                asset_df["area"] = asset_df.geometry.area
                asset_sums += ["area"] 
            
            min_groups = [item for sublist in [gr.split(",") for v,gr in asset_groups.items() if "max_" not in v] for item in sublist]
            max_groups = [item for sublist in [gr.split(",") for v,gr in asset_groups.items() if "min_" not in v] for item in sublist]
            # print (min_groups,max_groups)

            sum_df = asset_df.groupby(group_indexes)[asset_sums].sum().reset_index()
            min_df = asset_df.groupby(group_indexes)[min_groups].min().reset_index()
            min_df.rename(columns=dict([(v,f"min_{v}") for v in min_groups if "min_" not in v]) ,inplace=True)
            # min_groups = [f"min_{v}" for v in min_df. if "min_" not in v]
            max_df = asset_df.groupby(group_indexes)[max_groups].max().reset_index()
            max_df.rename(columns=dict([(v,f"max_{v}") for v in max_groups if "max_" not in v]) ,inplace=True)
            # max_groups = [f"max_{v}" for v in max_groups if "max_" not in v]

            min_max_df = pd.merge(sum_df,min_df,how="left",on=group_indexes).fillna(0)
            min_max_df = pd.merge(min_max_df,max_df,how="left",on=group_indexes).fillna(0)
            # print (min_max_df.columns)
            if "length" in min_max_df.columns.values.tolist():
                print (f"Total kms of {asset_info.asset_gpkg} {asset_info.asset_layer} {min_max_df.length.values.sum()}")
            for c in [v for v in min_max_df.columns.values.tolist() if v not in group_indexes]:
                # print (c)
                min_max_df[c] = min_max_df[c].apply(lambda v:f"{int(v):,}" if v > 10 else f"{round(v,2):,}")
            
            min_max_df.to_excel(output_wrtr,
                            sheet_name=f"{asset_info.asset_gpkg.split('_')[0]}_{asset_info.asset_layer.split('_')[0]}_{group_index}",
                            index=False)
            print (f"Done with layer {asset_info.asset_gpkg} {asset_info.asset_layer}")

    output_wrtr.save()

        


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
