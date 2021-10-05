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
from jamaica_sector_plotting_attributes import jamaica_sector_attributes, jamaica_currency_conversion
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


    # asset_data_path = os.path.join(processed_data_path,
    #                     "networks",
    #                     "transport")
    damage_data_path = os.path.join(output_data_path,
                        "direct_damages_summary")
    damage_string = "EAD_without_confidence_value" 
    damage_columns = ["EAD_undefended_mean"]
    damage_groupby = ["exposure_unit","damage_cost_unit","epoch"]
    # damages_filter_columns = ["rcp","epoch"]
    # damages_filter_values = [("baseline",2010),("4.5",2010)]
    damages_filter_columns = ["epoch"]
    damages_filter_values = [2010]
    sector_attributes = jamaica_sector_attributes()

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
            # sector["sector_gpkg"] = "potable_facilities_NWC.gpkg"
            # nodes = get_sector_layer(sector,asset_data_path,"node")
            # legend_font = 7
        else:
            edges = get_sector_layer(sector,asset_data_path,"edge")
            # nodes = get_sector_layer(sector,asset_data_path,"node")
            # legend_font = 12
        # edges = get_sector_layer(sector,asset_data_path,"edge")
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
            ax = get_axes(ax,extent=JAMAICA_EXTENT)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)

            if len(edges_damages) > 0:
                ax = line_map_plotting_colors_width(ax,edges_damages,"losses",
                                                    sector,
                                                    1.0,
                                                    legend_title,
                                                    "No risk/exposure/operation",
                                                    width_step = 200.0,
                                                    line_steps = 6,
                                                    # interpolation = 'log',
                                                    plot_title= f"{sector['sector_label']} multi-hazard Expected Annual Damages"
                                                    )
            
                save_fig(
                        os.path.join(
                            figures_data_path, 
                            f"{sector['sector_label'].lower().replace(' ','_')}_{sector['edge_layer']}_EAD.png"
                            )
                        )
            

        # nodes = get_sector_layer(sector,asset_data_path,"node")
        
        # if len(nodes) > 0:
        #     nodes_damages = get_asset_total_damage_values(sector,
        #                                         damage_data_path,damage_string,
        #                                         nodes,
        #                                         damages_filter_columns,
        #                                         damages_filter_values,
        #                                         damage_groupby,damage_columns,"node")
        #     nodes_damages["losses"] = nodes_damages.progress_apply(lambda x:convert_to_usd(x,"EAD_undefended_mean"),axis=1)

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
