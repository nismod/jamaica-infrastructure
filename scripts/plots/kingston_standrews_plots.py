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

def get_asset_total_damage_values(sector,damage_data_path,
                            damage_string,asset_dataframe,
                            damages_filter_columns,damages_filter_values,
                            damage_groupby,
                            damage_sum_columns,layer_key):
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
    return pd.merge(
                    asset_dataframe[[asset_id_column,sector[f"{layer_key}_classify_column"],"geometry"]],
                    damages,how="left",on=[asset_id_column]).fillna(0)
def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
    figures_data_path = config['paths']['figures']

    boundaries = gpd.read_file(os.path.join(processed_data_path,
                                                "boundaries",
                                                "admin_boundaries.gpkg"),
                                layer='admin1').to_crs(JAMAICA_GRID_EPSG)
    boundaries= boundaries[boundaries["PARISH"].isin(["Kingston","St. Andrew"])]
    bounds = boundaries.geometry.total_bounds # this gives your boundaries of the map as (xmin,ymin,xmax,ymax)
    kst_extent = (bounds[0]-10000,bounds[2]+10000,bounds[1],bounds[3])

    damage_data_path = os.path.join(output_data_path,
                        "direct_damages_summary")
    damage_string = "EAD_without_confidence_value" 
    damage_columns = ["EAD_undefended_mean"]
    damage_groupby = ["exposure_unit","damage_cost_unit","epoch"]
    damages_filter_columns = ["epoch"]
    damages_filter_values = [2010]
    sector_attributes = jamaica_sector_attributes()

    no_value_string = "No risk/exposure/operation"
    for sector in sector_attributes:
        asset_data_path = os.path.join(processed_data_path,
                        "networks",
                        sector["sector"])
        if sector["sector"] == "energy":
            legend_title = "Expected Annual Damages (US$/MW)"
        else:
            legend_title = "Expected Annual Damages (US$)"

        if sector["sector_label"] == "Potable water":
            sector["sector_gpkg"] = "pipelines_NWC.gpkg"
            edges = get_sector_layer(sector,asset_data_path,"edge")
        else:
            edges = get_sector_layer(sector,asset_data_path,"edge")
        if len(edges) > 0:
            edges_damages = get_asset_total_damage_values(sector,
                                                damage_data_path,damage_string,
                                                edges,
                                                damages_filter_columns,
                                                damages_filter_values,
                                                damage_groupby,damage_columns,"edge")
            edges_damages["losses"] = edges_damages.progress_apply(lambda x:convert_to_usd(x,"EAD_undefended_mean"),axis=1)

            
            """plot the damage results
            """
            fig, ax = plt.subplots(1,1,
                            subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                            figsize=(12,8),
                            dpi=500)
            ax = get_axes(ax,extent=kst_extent)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
            scale_bar_and_direction(ax,scalebar_distance=5)

            if len(edges_damages) > 0:
                if sector["sector_label"] in ["Roads","Railways","Potable water"]:
                    ax = line_map_plotting_colors_width(
                                                        ax,edges_damages,"losses",
                                                        legend_label=legend_title,
                                                        no_value_label=no_value_string,
                                                        width_step=10,
                                                        interpolation="log",
                                                        plot_title=f"{sector['sector_label']} multi-hazard Expected Annual Damages"
                                                        )

                else:
                    # print ("* Don't plot")
                    ax = line_map_plotting_colors_width(
                                                        ax,edges_damages,"losses",
                                                        edge_classify_column=sector["edge_classify_column"],
                                                        edge_categories=sector["edge_categories"],
                                                        edge_colors=sector["edge_categories_colors"],
                                                        edge_labels=sector["edge_categories_labels"],
                                                        edge_zorder=sector["edge_categories_zorder"],
                                                        legend_label=legend_title,
                                                        no_value_label=no_value_string,
                                                        line_steps=6,
                                                        width_step=20,
                                                        # interpolation="log",
                                                        plot_title=f"{sector['sector_label']} multi-hazard Expected Annual Damages",
                                                        )
            
                save_fig(
                        os.path.join(
                            figures_data_path, 
                            f"k_st_{sector['sector_label'].lower().replace(' ','_')}_{sector['edge_layer']}_EAD.png"
                            )
                        )
            

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
                port_damages = get_asset_total_damage_values(port_sector,
                                                    damage_data_path,damage_string,
                                                    port_nodes,
                                                    damages_filter_columns,
                                                    damages_filter_values,
                                                    damage_groupby,damage_columns,"area")
                airport_sector = [s for s in sector_attributes if s["sector_label"] == "Airports"][0]
                airport_nodes = get_sector_layer(airport_sector,asset_data_path,"area")
                airport_damages = get_asset_total_damage_values(airport_sector,
                                                    damage_data_path,damage_string,
                                                    airport_nodes,
                                                    damages_filter_columns,
                                                    damages_filter_values,
                                                    damage_groupby,damage_columns,"area")
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
                nodes_damages = get_asset_total_damage_values(sector,
                                                    damage_data_path,damage_string,
                                                    nodes,
                                                    damages_filter_columns,
                                                    damages_filter_values,
                                                    damage_groupby,damage_columns,"node")
            nodes_damages["losses"] = nodes_damages.progress_apply(lambda x:convert_to_usd(x,"EAD_undefended_mean"),axis=1)

            """plot the damage results
            """
            fig, ax = plt.subplots(1,1,
                            subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                            figsize=(12,8),
                            dpi=500)
            ax = get_axes(ax,extent=kst_extent)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
            scale_bar_and_direction(ax,scalebar_distance=5)

            if len(nodes_damages) > 0:
                if sector["sector_label"] in ["Roads","Railways","Potable water","Irrigation","Wastewater Treatment"]:
                    ax = point_map_plotting_colors_width(
                                                        ax,nodes_damages,"losses",
                                                        legend_label=legend_title,
                                                        no_value_label=no_value_string,
                                                        width_step=10,
                                                        interpolation="log",
                                                        plot_title=f"{sector['sector_label']} multi-hazard Expected Annual Damages"
                                                        )

                else:
                    # print ("* Don't plot")
                    ax = point_map_plotting_colors_width(
                                                        ax,nodes_damages,"losses",
                                                        point_classify_column=sector["node_classify_column"],
                                                        point_categories=sector["node_categories"],
                                                        point_colors=sector["node_categories_colors"],
                                                        point_labels=sector["node_categories_labels"],
                                                        point_zorder=sector["node_categories_zorder"],
                                                        legend_label=legend_title,
                                                        no_value_label=no_value_string,
                                                        width_step=10,
                                                        interpolation="log",
                                                        plot_title=f"{sector['sector_label']} multi-hazard Expected Annual Damages",
                                                        )
            
                if len(edges) > 0:
                    ax = plot_line_assets(ax,JAMAICA_GRID_EPSG,edges,"#969696",0.5,4)
                save_fig(
                        os.path.join(
                            figures_data_path, 
                            f"k_st_{sector['sector_label'].lower().replace(' ','_')}_{sector['node_layer']}_EAD.png"
                            )
                        )



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
