"""Generate hazard-damage curves
"""
import os
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from plot_utils import *
from tqdm import tqdm
tqdm.pandas()
from jamaica_sector_plotting_attributes import (
                                                jamaica_sector_attributes, 
                                                jamaica_currency_conversion, 
                                                jamaica_port_and_airport_nodes
                                                )
from collections import namedtuple, OrderedDict

JAMAICA_EXTENT = (598251, 838079, 610353, 714779)
JAMAICA_GRID_EPSG = 3448

def get_asset_total_damage_values(sector,damage_data_path,
                            damage_string,asset_dataframe,
                            damages_filter_columns,damages_filter_values,
                            damage_groupby,
                            damage_sum_columns,layer_key,damage_divisor=1.0):
    asset_id_column = sector[f"{layer_key}_id_column"]
    asset_filter_column = sector[f"{layer_key}_damage_filter_column"]
    asset_filter_list = sector[f"{layer_key}_damage_categories"]
    damages = pd.read_csv(
                    os.path.join(
                        damage_data_path,
                        f"{sector['sector_gpkg'].replace('.gpkg','')}_{sector[f'{layer_key}_layer']}_{damage_string}.csv"
                        )
                    )
    damages = damages.set_index(damages_filter_columns)
    damages = damages[damages.index.isin(damages_filter_values)].reset_index()
    if asset_filter_column is not None:
        asset_ids = asset_dataframe[asset_dataframe[asset_filter_column].isin(asset_filter_list)][asset_id_column].values.tolist()
        damages = damages[damages[asset_id_column].isin(asset_ids)]
    
    eael_values = [c for c in damages.columns.values.tolist() if "EAEL_" in c]
    if len(eael_values) == 0:
        damage_sum_columns = [c for c in damage_sum_columns if c != "EAEL"]

    dc = [c for c in damages.columns.values.tolist() if "designed_protection_mean" in c]
    if len(dc) > 0:
        damage_sum_columns = [f"{d}_designed_protection_mean" for d in damage_sum_columns]
    else:
        damage_sum_columns = [f"{d}_undefended_mean" for d in damage_sum_columns]

    damages = damages.groupby(
                    [asset_id_column] + damage_groupby,dropna=False
                    ).agg(
                        dict(
                            zip(
                                damage_sum_columns,["sum"]*len(damage_sum_columns)
                                )
                            )
                        ).reset_index()
    if damage_divisor != 1:
        for d in damage_sum_columns:
            damages[d] = 1.0*damages[d]/damage_divisor
    return pd.merge(
                    asset_dataframe[[asset_id_column,sector[f"{layer_key}_classify_column"],"geometry"]],
                    damages,how="left",on=[asset_id_column]).fillna(0)

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
    figures_data_path = config['paths']['figures']


    damage_data_path = os.path.join(output_data_path,
                        "direct_damages_summary")
    damage_string = "EAD_EAEL" 
    damage_types = ["EAD","EAEL"]
    damage_legends = ["Expected Annual Damages (J$ million)",
                    "Expected Annual Losses (J$ million/day)"]
    damage_titles = ["Expected Annual Damages",
                    "Expected Annual Losses"]
    damage_groupby = ["damage_cost_unit","economic_loss_unit","epoch"]
    damages_filter_columns = ["epoch"]
    damages_filter_values = [2010]
    sector_attributes = jamaica_sector_attributes()

    no_value_string = "No risk/exposure/operation"
    for sector in sector_attributes:
        print (f"* Starting sector {sector['sector_label']}")
        asset_data_path = os.path.join(processed_data_path,
                        "networks",
                        sector["sector"])  

        if sector["sector_label"] == "Potable water":
            sector["sector_gpkg"] = "pipelines_NWC.gpkg"
            edges = get_sector_layer(sector,asset_data_path,"edge")
        else:
            edges = get_sector_layer(sector,asset_data_path,"edge")

        if len(edges) > 0:
            edges_damages = get_asset_total_damage_values(sector,
                                                damage_data_path,
                                                damage_string,
                                                edges,
                                                damages_filter_columns,
                                                damages_filter_values,
                                                damage_groupby,damage_types,"edge",damage_divisor=1.0e6)

            if len(edges_damages) > 0:
                for idx, (dt,lt,titles) in enumerate(list(zip(damage_types,damage_legends,damage_titles))):
                    damage_column = [c for c in edges_damages if dt in c]
                    if len(damage_column) > 0:
                        damage_column = damage_column[0]
                        if edges_damages[damage_column].sum() > 0:
                            """plot the damage results
                            """
                            fig, ax = plt.subplots(1,1,
                                            subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                                            figsize=(12,8),
                                            dpi=500)
                            ax = get_axes(ax,extent=JAMAICA_EXTENT)
                            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
                            scale_bar_and_direction(ax)

                            if sector["sector_label"] in ["Roads","Potable water"]:
                                ax = line_map_plotting_colors_width(
                                                                    ax,edges_damages,damage_column,
                                                                    legend_label=lt,
                                                                    no_value_label=no_value_string,
                                                                    width_step=40,
                                                                    # interpolation="log",
                                                                    interpolation="fisher-jenks",
                                                                    plot_title=f"{sector['sector_label']} multi-hazard {titles}"
                                                                    )

                            else:
                                ax = line_map_plotting_colors_width(
                                                                    ax,edges_damages,damage_column,
                                                                    edge_classify_column=sector["edge_classify_column"],
                                                                    edge_categories=sector["edge_categories"],
                                                                    edge_colors=sector["edge_categories_colors"],
                                                                    edge_labels=sector["edge_categories_labels"],
                                                                    edge_zorder=sector["edge_categories_zorder"],
                                                                    legend_label=lt,
                                                                    no_value_label=no_value_string,
                                                                    line_steps=6,
                                                                    width_step=200,
                                                                    # interpolation="log",
                                                                    interpolation="fisher-jenks",
                                                                    plot_title=f"{sector['sector_label']} multi-hazard {titles}"
                                                                    )
                    
                            save_fig(
                                    os.path.join(
                                        figures_data_path, 
                                        f"{sector['sector_label'].lower().replace(' ','_')}_{sector['edge_layer']}_{dt}.png"
                                        )
                                    )
            del (edges_damages)

        
        if sector["sector_label"] == "Potable water":
            sector["sector_gpkg"] = "potable_facilities_NWC.gpkg"
            nodes = get_sector_layer(sector,asset_data_path,"node")
        elif sector["sector_label"] == "Ports and Airports":
            nodes = jamaica_port_and_airport_nodes()
        else:
            nodes = get_sector_layer(sector,asset_data_path,"node")
        if len(nodes) > 0:
            if sector["sector_label"] == "Ports and Airports":
                port_sector = [s for s in sector_attributes if s["sector_label"] == "Ports"][0]
                port_nodes = get_sector_layer(port_sector,asset_data_path,"area")
                port_damages = get_asset_total_damage_values(port_sector,damage_data_path,damage_string,
                                                    port_nodes,
                                                    damages_filter_columns,
                                                    damages_filter_values,
                                                    damage_groupby,damage_types,"area",damage_divisor=1.0e6)
                airport_sector = [s for s in sector_attributes if s["sector_label"] == "Airports"][0]
                airport_nodes = get_sector_layer(airport_sector,asset_data_path,"area")
                airport_damages = get_asset_total_damage_values(airport_sector,damage_data_path,damage_string,
                                                    airport_nodes,
                                                    damages_filter_columns,
                                                    damages_filter_values,
                                                    damage_groupby,damage_types,"area",damage_divisor=1.0e6)
                nodes_damages = pd.concat(
                                            [
                                            port_damages.drop("geometry",axis=1),
                                            airport_damages.drop("geometry",axis=1)
                                            ],
                                            axis=0,
                                            ignore_index=True
                                        )
                nodes_damages = pd.merge(
                                            nodes,
                                            nodes_damages,how='left',on=[sector["node_id_column"]]).fillna(0)
            else:
                nodes_damages = get_asset_total_damage_values(sector,damage_data_path,damage_string,
                                                    nodes,
                                                    damages_filter_columns,
                                                    damages_filter_values,
                                                    damage_groupby,damage_types,"node",damage_divisor=1.0e6)
            
            if len(nodes_damages) > 0:
                for idx, (dt,lt,titles) in enumerate(list(zip(damage_types,damage_legends,damage_titles))):
                    damage_column = [c for c in nodes_damages if dt in c]
                    if len(damage_column) > 0:
                        damage_column = damage_column[0]
                        if nodes_damages[damage_column].sum() > 0:
                            """plot the damage results
                            """
                            fig, ax = plt.subplots(1,1,
                                            subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                                            figsize=(12,8),
                                            dpi=500)
                            ax = get_axes(ax,extent=JAMAICA_EXTENT)
                            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
                            scale_bar_and_direction(ax)

                            if sector["sector_label"] in ["Roads","Ports and Airports","Railways","Potable water","Irrigation","Wastewater Treatment"]:
                                ax = point_map_plotting_colors_width(
                                                                    ax,nodes_damages,damage_column,
                                                                    legend_label=lt,
                                                                    no_value_label=no_value_string,
                                                                    width_step=20,
                                                                    # interpolation="log",
                                                                    interpolation="fisher-jenks",
                                                                    plot_title=f"{sector['sector_label']} multi-hazard {titles}"
                                                                    )

                            else:
                                ax = point_map_plotting_colors_width(
                                                                    ax,nodes_damages,damage_column,
                                                                    point_classify_column=sector["node_classify_column"],
                                                                    point_categories=sector["node_categories"],
                                                                    point_colors=sector["node_categories_colors"],
                                                                    point_labels=sector["node_categories_labels"],
                                                                    point_zorder=sector["node_categories_zorder"],
                                                                    legend_label=lt,
                                                                    no_value_label=no_value_string,
                                                                    width_step=20,
                                                                    # interpolation="log",
                                                                    interpolation="fisher-jenks",
                                                                    plot_title=f"{sector['sector_label']} multi-hazard {titles}",
                                                                    )
                
                            if len(edges) > 0:
                                ax = plot_line_assets(ax,JAMAICA_GRID_EPSG,edges,"#969696",0.5,4)
                            save_fig(
                                    os.path.join(
                                        figures_data_path, 
                                        f"{sector['sector_label'].lower().replace(' ','_')}_{sector['node_layer']}_{dt}.png"
                                        )
                                    )
            del (nodes_damages)



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
