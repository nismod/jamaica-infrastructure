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

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
    figures_data_path = config['paths']['figures']


    no_value_string = "No traffic counts"
    asset_data_path = os.path.join(processed_data_path,
                    "networks",
                    "transport")

    legend_title = "Traffic count (vehicles/day)"    
    
    edges = gpd.read_file(os.path.join(asset_data_path,"roads.gpkg"),layer="edges")

    """plot the damage results
    """
    fig, ax = plt.subplots(1,1,
                    subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                    figsize=(12,8),
                    dpi=500)
    ax = get_axes(ax,extent=JAMAICA_EXTENT)
    plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=False)
    scale_bar_and_direction(ax)

    ax = line_map_plotting_colors_width(
                                        ax,edges,"traffic_count",
                                        legend_label=legend_title,
                                        no_value_label=no_value_string,
                                        edge_colors=['#8c96c6','#88419d','#fb6a4a','#d7301f','#7f0000'],
                                        width_step=60,
                                        interpolation="log",
                                        plot_title="Roads - Traffic counts"
                                        )
    
    save_fig(
            os.path.join(
                figures_data_path, 
                "roads_traffic_counts.png"
                )
            )



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
