"""Generate hazard-damage curves
"""
import os
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
from collections import OrderedDict
from summarise_utils import *
from jamaica_sector_attributes import (jamaica_sector_attributes,jamaica_currency_conversion)
from tqdm import tqdm
tqdm.pandas()

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

    if isinstance(electricity_damage_column,list):
        for ed in electricity_damage_column:
            electricity_damage_df[
                        ed
                        ] = electricity_damage_df[
                                                    electricity_capacity_column
                                                    ]*electricity_damage_df[
                                                                ed
                                                                ]
    else:
        electricity_damage_df[
                        electricity_damage_column
                        ] = electricity_damage_df[
                                                    electricity_capacity_column
                                                    ]*electricity_damage_df[
                                                                electricity_damage_column
                                                                ]
    electricity_cap_df.drop(electricity_capacity_column,axis=1)
    return electricity_damage_df

def filter_damaged_assets(sector,damage_data_path,
                            damage_string,asset_dataframe,
                            layer_key):
    asset_id_column = sector[f"{layer_key}_id_column"]
    asset_filter_column = sector[f"{layer_key}_damage_filter_column"]
    asset_filter_list = sector[f"{layer_key}_damage_categories"]
    damages = pd.read_csv(
                    os.path.join(
                        damage_data_path,
                        f"{sector['sector_gpkg'].replace('.gpkg','')}_{sector[f'{layer_key}_layer']}_{damage_string}.csv"
                        )
                    )
    
    if asset_filter_column is not None:
        asset_ids = asset_dataframe[asset_dataframe[asset_filter_column].isin(asset_filter_list)][asset_id_column].values.tolist()
        damages = damages[damages[asset_id_column].isin(asset_ids)]
    
    return damages

def get_damage_columns(damage_dictionary,damage_dataframe):
    if damage_dictionary["damage_type"] == "direct_damages":
        damage_columns = ["exposure"]
    else:
        damage_columns = []
    for dg_c in damage_dictionary["damage_columns"]:
        damage_columns += [c for c in damage_dataframe.columns.values.tolist() if dg_c in c]
    return damage_columns

def damages_grouped(damages,damage_groupby,damage_sum_columns,sector_name=None,layer_name=None,current_convert="J$"):
    damages = damages.groupby(
                    damage_groupby,dropna=False
                    ).agg(
                        dict(
                            zip(
                                damage_sum_columns,["sum"]*len(damage_sum_columns)
                                )
                            )
                        ).reset_index() 
    damage_cost_columns = [d for d in damage_sum_columns if d != "exposure"]
    if current_convert == "J$":
        if (
            "$US" in damages["damage_cost_unit"].values[0]
            ) or (
                    "US$" in damages["damage_cost_unit"].values[0]
                ) or (
                        "USD" in damages["damage_cost_unit"].values[0]
                    ):
            for d in damage_cost_columns:
                damages[d] = (1.0/jamaica_currency_conversion())*damages[d]
    damages = damages.drop("damage_cost_unit",axis=1)
    damages["sector"] = sector_name
    damages["subsector"] = layer_name
    return damages

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
    figures_data_path = config['paths']['figures']


    damage_data_path = os.path.join(output_data_path,
                        "direct_damages_summary")

    damage_cases = [
                    {
                        "damage_type":"direct_damages",
                        "damage_columns":["direct_damage_cost","direct_reopen_cost"],
                        "index_columns":["damage_cost_unit","hazard","rp","rcp","epoch"]
                    },
                    {
                        "damage_type": "EAD",
                        "damage_columns":["EAD","EAR"],
                        "index_columns":["damage_cost_unit","hazard","rcp","epoch"]
                    }
                ]
    
    sector_attributes = jamaica_sector_attributes()

    for dg_cs in damage_cases:
        damages_combined = []
        damage_string = f"{dg_cs['damage_type']}_without_confidence_value"
        csv_outputs = os.path.join(damage_data_path,
                            f"hazard_rcp_epoch_{dg_cs['damage_type']}_without_confidence_value_JD.csv")
        for sector in sector_attributes:
            asset_data_path = os.path.join(processed_data_path,
                            "networks",
                            sector["sector"])
            if sector["sector_label"] == "Potable water":
                sector["sector_gpkg"] = "pipelines_NWC.gpkg"
                edges = get_sector_layer(sector,asset_data_path,"edge")
            else:
                edges = get_sector_layer(sector,asset_data_path,"edge")

            if len(edges) > 0:
                edges_damages = filter_damaged_assets(sector,
                                                    damage_data_path,
                                                    damage_string,
                                                    edges,
                                                    "edge")
                
                damage_columns = get_damage_columns(dg_cs,edges_damages)
                if sector["sector"] == "energy":
                    edges_damages = get_electricity_costs(asset_data_path,
                                edges_damages,
                                sector["edge_id_column"],
                                damage_columns,
                                "max",
                                sector["edge_layer"]
                                )
                damages_combined.append(
                                    damages_grouped(
                                            edges_damages,
                                            dg_cs["index_columns"],
                                            damage_columns,
                                            sector["sector_label"],
                                            sector["edge_layer"]
                                            )
                                    )

            if sector["sector_label"] == "Potable water":
                sector["sector_gpkg"] = "potable_facilities_NWC.gpkg"
                nodes = get_sector_layer(sector,asset_data_path,"node")
            else:
                nodes = get_sector_layer(sector,asset_data_path,"node")
            if len(nodes) > 0:
                nodes_damages = filter_damaged_assets(sector,
                                                damage_data_path,
                                                damage_string,
                                                nodes,
                                                "node")
                damage_columns = get_damage_columns(dg_cs,nodes_damages)
                if sector["sector"] == "energy":
                    nodes_damages = get_electricity_costs(asset_data_path,
                                nodes_damages,
                                sector["node_id_column"],
                                damage_columns,
                                "capacity",
                                sector["node_layer"]
                                )
                damages_combined.append(
                                damages_grouped(
                                            nodes_damages,
                                            dg_cs["index_columns"],
                                            damage_columns,
                                            sector["sector_label"],
                                            sector["node_layer"]
                                            )
                                )
            areas = get_sector_layer(sector,asset_data_path,"area")
            if len(areas) > 0:
                areas_damages = filter_damaged_assets(sector,
                                                damage_data_path,
                                                damage_string,
                                                areas,
                                                "area")
                damage_columns = get_damage_columns(dg_cs,areas_damages)
                damages_combined.append(
                                    damages_grouped(
                                        areas_damages,
                                        dg_cs["index_columns"],
                                        damage_columns,
                                        sector["sector_label"],
                                        sector["area_layer"]
                                        )
                                    )
            print (f"* Done with {sector['sector_label']}")
        damages_combined = pd.concat(damages_combined,axis=0,ignore_index=True).fillna(0)
        damages_combined.to_csv(csv_outputs,index=False)

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
