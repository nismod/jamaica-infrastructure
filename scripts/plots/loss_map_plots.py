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
                                                jamaica_port_and_airport_nodes,
                                                jamaica_port_nodes
                                                )
from collections import namedtuple, OrderedDict

JAMAICA_EXTENT = (598251, 838079, 610353, 714779)
JAMAICA_GRID_EPSG = 3448

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
    damages.loc[damages["epoch"] == 2010,"rcp"] = "baseline"
    damages = damages.set_index(damages_filter_columns)
    # print (damages)
    # print (list(set(damages.index.values.tolist())))
    damages = damages[damages.index.isin(damages_filter_values)].reset_index()
    # print (damages)
    # damages = damages[damages.set_index(damages_filter_columns).index.isin(damages_filter_values)].reset_index()
    if asset_filter_column is not None:
        asset_ids = asset_dataframe[asset_dataframe[asset_filter_column].isin(asset_filter_list)][asset_id_column].values.tolist()
        damages = damages[damages[asset_id_column].isin(asset_ids)]
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

def get_asset_total_loss_values(sector,asset_dataframe,losses_data_path,layer_key,damage_sum_columns,damage_divisor=1.0):
    asset_id_column = sector[f"{layer_key}_id_column"]
    asset_filter_column = sector[f"{layer_key}_damage_filter_column"]
    asset_filter_list = sector[f"{layer_key}_damage_categories"]
    sector_string = sector[f"{layer_key}_loss_result"]
    damages = pd.read_csv(
                    os.path.join(
                        losses_data_path,
                        f"single_point_failure_{sector_string}_economic_losses.csv"
                        )
                    )[[asset_id_column,"economic_loss"]]
    if asset_filter_column is not None:
        asset_ids = asset_dataframe[asset_dataframe[asset_filter_column].isin(asset_filter_list)][asset_id_column].values.tolist()
        damages = damages[damages[asset_id_column].isin(asset_ids)]
    
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


    losses_data_path = os.path.join(output_data_path,
                                    "economic_losses",
                                    "single_failure_scenarios"
                                    )
    loss_columns = ["economic_loss"]
    sector_attributes = jamaica_sector_attributes()

    no_value_string = "No Estimated Loss"
    for sector in sector_attributes:
        if sector["sector_label"] not in ["Airports","Ports and Airports","Wastewater"]:
            asset_data_path = os.path.join(processed_data_path,
                            "networks",
                            sector["sector"])
            if sector["sector_label"] == "Railways":
                damage_divisor = 1.0
                legend_title = "Losses (J$/day)" 
            else:
                damage_divisor = 1.0e6
                legend_title = "Losses (J$ million/day)"    
            if sector["sector_label"] == "Potable water":
                sector["sector_gpkg"] = "pipelines_NWC.gpkg"
                edges = get_sector_layer(sector,asset_data_path,"edge")
            else:
                edges = get_sector_layer(sector,asset_data_path,"edge")

            if len(edges) > 0:
                edges_damages = get_asset_total_loss_values(sector,edges,
                                                    losses_data_path,"edge",loss_columns,damage_divisor=damage_divisor)
                edges_damages.to_file(os.path.join(losses_data_path, sector["sector_gpkg"]),layer="edges",driver="GPKG")

                """plot the damage results
                """
                fig, ax = plt.subplots(1,1,
                                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                                figsize=(12,8),
                                dpi=500)
                ax = get_axes(ax,extent=JAMAICA_EXTENT)
                plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
                scale_bar_and_direction(ax)

                if len(edges_damages) > 0:
                    if sector["sector_label"] in ["Roads","Potable water"]:
                        ax = line_map_plotting_colors_width(
                                                            ax,edges_damages,"economic_loss",
                                                            legend_label=legend_title,
                                                            no_value_label=no_value_string,
                                                            width_step=40,
                                                            interpolation="fisher-jenks",
                                                            line_steps=8,
                                                            # interpolation="log",
                                                            plot_title=f"{sector['sector_label']} criticality impacts"
                                                            )

                    else:
                        # print ("* Don't plot")
                        ax = line_map_plotting_colors_width(
                                                            ax,edges_damages,"economic_loss",
                                                            edge_classify_column=sector["edge_classify_column"],
                                                            edge_categories=sector["edge_categories"],
                                                            edge_colors=sector["edge_categories_colors"],
                                                            edge_labels=sector["edge_categories_labels"],
                                                            edge_zorder=sector["edge_categories_zorder"],
                                                            legend_label=legend_title,
                                                            no_value_label=no_value_string,
                                                            line_steps=6,
                                                            width_step=200,
                                                            # interpolation="log",
                                                            interpolation="fisher-jenks",
                                                            plot_title=f"{sector['sector_label']} criticality impacts",
                                                            )
                
                    save_fig(
                            os.path.join(
                                figures_data_path, 
                                f"{sector['sector_label'].lower().replace(' ','_')}_{sector['edge_layer']}_losses_JD.png"
                                )
                            )
                    del (edges_damages)
            if sector["sector_label"] == "Potable water":
                sector["sector_gpkg"] = "potable_facilities_NWC.gpkg"
                nodes = get_sector_layer(sector,asset_data_path,"node")
            elif sector["sector_label"] == "Ports and Airports":
                nodes = jamaica_port_and_airport_nodes()
            elif sector["sector_label"] == "Ports":
                nodes = jamaica_port_nodes()
            else:
                nodes = get_sector_layer(sector,asset_data_path,"node")
            if len(nodes) > 0:
                if sector["sector_label"] == "Ports and Airports":
                    port_sector = [s for s in sector_attributes if s["sector_label"] == "Ports"][0]
                    port_nodes = get_sector_layer(port_sector,asset_data_path,"area")
                    port_damages = get_asset_total_loss_values(port_sector,port_nodes,
                                                    losses_data_path,"area",loss_columns,damage_divisor=damage_divisor)
                    airport_sector = [s for s in sector_attributes if s["sector_label"] == "Airports"][0]
                    airport_nodes = get_sector_layer(airport_sector,asset_data_path,"area")
                    airport_damages = get_asset_total_loss_values(airport_sector,airport_nodes,
                                                    losses_data_path,"area",loss_columns,damage_divisor=damage_divisor)
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
                    nodes_damages = get_asset_total_loss_values(sector,nodes,
                                                    losses_data_path,"node",loss_columns,damage_divisor=damage_divisor)

                nodes_damages.to_file(os.path.join(losses_data_path, sector["sector_gpkg"]),layer="nodes",driver="GPKG")
                """plot the damage results
                """

                fig, ax = plt.subplots(1,1,
                                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                                figsize=(12,8),
                                dpi=500)
                ax = get_axes(ax,extent=JAMAICA_EXTENT)
                plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
                scale_bar_and_direction(ax)

                if len(nodes_damages) > 0:
                    if sector["sector_label"] in ["Roads","Ports","Ports and Airports","Railways","Potable water","Irrigation","Wastewater Treatment"]:
                        ax = point_map_plotting_colors_width(
                                                            ax,nodes_damages,"economic_loss",
                                                            legend_label=legend_title,
                                                            no_value_label=no_value_string,
                                                            width_step=20,
                                                            # interpolation="log",
                                                            # interpolation="equal bins",
                                                            interpolation="fisher-jenks",
                                                            plot_title=f"{sector['sector_label']} criticality impacts"
                                                            )

                    else:
                        ax = point_map_plotting_colors_width(
                                                            ax,nodes_damages,"economic_loss",
                                                            point_classify_column=sector["node_classify_column"],
                                                            point_categories=sector["node_categories"],
                                                            point_colors=sector["node_categories_colors"],
                                                            point_labels=sector["node_categories_labels"],
                                                            point_zorder=sector["node_categories_zorder"],
                                                            legend_label=legend_title,
                                                            no_value_label=no_value_string,
                                                            width_step=20,
                                                            # interpolation="log",
                                                            # interpolation="equal bins",
                                                            interpolation="fisher-jenks",
                                                            plot_title=f"{sector['sector_label']} criticality impacts",
                                                            )
                
                    if len(edges) > 0:
                        ax = plot_line_assets(ax,JAMAICA_GRID_EPSG,edges,"#969696",0.5,4)
                    save_fig(
                            os.path.join(
                                figures_data_path, 
                                f"{sector['sector_label'].lower().replace(' ','_')}_{sector['node_layer']}_losses_JD.png"
                                )
                            )
                    del (nodes_damages)



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
