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

def filter_damaged_assets(sector,damage_data_path,
                            damage_string,asset_dataframe,
                            layer_key):
    asset_id_column = sector[f"{layer_key}_id_column"]
    asset_filter_column = sector[f"{layer_key}_damage_filter_column"]
    asset_filter_list = sector[f"{layer_key}_damage_categories"]
    file = os.path.join(
                            damage_data_path,
                            f"{sector['sector_gpkg'].replace('.gpkg','')}_{sector[f'{layer_key}_layer']}_{damage_string}.csv"
                            )
    if os.path.isfile(file) is True:
        damages = pd.read_csv(file)
        if asset_filter_column is not None:
            asset_ids = asset_dataframe[asset_dataframe[asset_filter_column].isin(asset_filter_list)][asset_id_column].values.tolist()
            damages = damages[damages[asset_id_column].isin(asset_ids)]
    else:
        damages = []
    
    return damages

def get_damage_columns(damage_dictionary,damage_dataframe):

    # if damage_dictionary["damage_type"] == "direct_damages":
    #     damage_columns = ["exposure"]
    # else:
    #     damage_columns = []
    # for dg_c in damage_dictionary["damage_columns"]:
    #     damage_columns += [c for c in damage_dataframe.columns.values.tolist() if dg_c in c]
    # return damage_columns

    dc = damage_dictionary["damage_columns"].copy()
    damage_sum_columns = []
    for d in dc:
        damage_sum_columns += [c for c in damage_dataframe.columns.values.tolist() if f"{d}_" in c]

    dc = [c for c in damage_sum_columns if "designed_protection" in c]
    if len(dc) > 0:
        damage_sum_columns = dc

    return damage_sum_columns

def damages_grouped(damages,damage_groupby,damage_sum_columns,sector_name=None,layer_name=None,current_convert="J$"):
    damages.columns = damages.columns.map(str)
    damages["rcp"] = damages["rcp"].map(str)
    damages["sector"] = sector_name
    damages["subsector"] = layer_name
    damages = damages.groupby(
                    damage_groupby + ["sector","subsector"],dropna=False
                    ).agg(
                        dict(
                            zip(
                                damage_sum_columns,["sum"]*len(damage_sum_columns)
                                )
                            )
                        ).reset_index() 
    # damage_cost_columns = [d for d in damage_sum_columns if d != "exposure"]
    # if current_convert == "J$":
    #     if (
    #         "$US" in damages["damage_cost_unit"].values[0]
    #         ) or (
    #                 "US$" in damages["damage_cost_unit"].values[0]
    #             ) or (
    #                     "USD" in damages["damage_cost_unit"].values[0]
    #                 ):
    #         for d in damage_cost_columns:
    #             damages[d] = (1.0/jamaica_currency_conversion())*damages[d]
    # damages = damages.drop("damage_cost_unit",axis=1)
    damages.columns = [str(c).replace("undefended","").replace("designed_protection","") for c in damages.columns.values.tolist()]
    return damages

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']


    damage_data_path = os.path.join(output_data_path,
                        "direct_damages_summary")

    summary_results = os.path.join(output_data_path,"damage_loss_sums")
    if os.path.exists(summary_results) == False:
        os.mkdir(summary_results)
    
    damage_cases = [
                    # {
                    #     "damage_type":"direct_damages",
                    #     "damage_columns":["direct_damage_cost","direct_reopen_cost"],
                    #     "index_columns":["damage_cost_unit","hazard","rp","rcp","epoch"]
                    # },
                    {
                        "damage_type": "EAD_EAEL",
                        "damage_data_path": os.path.join(output_data_path,"direct_damages_summary"),
                        "damage_columns":["EAD","EAEL"],
                        "index_columns":["damage_cost_unit","economic_loss_unit","hazard","rcp","epoch"]
                    }
                ]
    
    sector_attributes = jamaica_sector_attributes()

    for dg_cs in damage_cases:
        damages_combined = []
        damage_string = dg_cs['damage_type']
        csv_outputs = os.path.join(summary_results,
                            f"hazard_rcp_epoch_{dg_cs['damage_type']}_totals.csv")
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
                if len(edges_damages) > 0:
                    damage_columns = get_damage_columns(dg_cs,edges_damages)
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
                if len(nodes_damages) > 0:
                    damage_columns = get_damage_columns(dg_cs,nodes_damages)
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
                if len(areas_damages) > 0:
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

    start_year = 2019
    end_year = 2100
    damage_columns = [str(n) for n in np.arange(start_year,end_year+1,1)]
    damage_data_path = os.path.join(output_data_path,
                        "loss_damage_timeseries")
    csv_outputs = os.path.join(summary_results,
                                "hazard_rcp_timeseries_totals.csv")
    damages_cases = ["EAD","EAEL"]
    damages_types = ["amin","mean","amax"]
    damages_combined = []
    for dg_cs in damages_cases:
        for dt in damages_types:
            damage_string = f"{dg_cs}_timeseries_{dt}"
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
                    if len(edges_damages) > 0:
                        edges_damages["risk_type"] = dg_cs
                        edges_damages["val_type"] = dt
                        dam = damages_grouped(
                                                    edges_damages,
                                                    ["hazard","rcp","risk_type","val_type"],
                                                    damage_columns,
                                                    sector["sector_label"],
                                                    sector["edge_layer"]
                                                    )
                        damages_combined.append(dam)

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
                    if len(nodes_damages) > 0:
                        nodes_damages["risk_type"] = dg_cs
                        nodes_damages["val_type"] = dt
                        dam = damages_grouped(
                                                nodes_damages,
                                                ["hazard","rcp","risk_type","val_type"],
                                                damage_columns,
                                                sector["sector_label"],
                                                sector["node_layer"]
                                                )
                    
                        damages_combined.append(dam)

                areas = get_sector_layer(sector,asset_data_path,"area")
                if len(areas) > 0:
                    areas_damages = filter_damaged_assets(sector,
                                                    damage_data_path,
                                                    damage_string,
                                                    areas,
                                                    "area")
                    if len(areas_damages) > 0:
                        areas_damages["risk_type"] = dg_cs
                        areas_damages["val_type"] = dt
                        dam = damages_grouped(
                                                areas_damages,
                                                ["hazard","rcp","risk_type","val_type"],
                                                damage_columns,
                                                sector["sector_label"],
                                                sector["area_layer"]
                                                )

                        damages_combined.append(dam)
                print (f"* Done with {sector['sector_label']} {dt} {dg_cs}")
    damages_combined = pd.concat(damages_combined,axis=0,ignore_index=True).fillna(0)
    damages_combined.to_csv(csv_outputs,index=False)

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
