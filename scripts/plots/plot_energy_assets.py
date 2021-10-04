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
from jamaica_sector_plotting_attributes import jamaica_sector_attributes

JAMAICA_EXTENT = (598251, 838079, 610353, 714779)
JAMAICA_GRID_EPSG = 3448

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    figures_data_path = config['paths']['figures']


    energy_data_path = os.path.join(processed_data_path,
                        "networks",
                        "energy")
    sector_attributes = jamaica_sector_attributes()
    for sector in sector_attributes:
        if sector["sector_label"] in ["Energy"]:
            edges = get_sector_layer(sector,energy_data_path,"edge")
            nodes = get_sector_layer(sector,energy_data_path,"node")
            # print (edges, nodes)

            legend_handles = []
            fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(12,8),
                dpi=500)
            ax = get_axes(ax,extent=JAMAICA_EXTENT)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)

            if len(edges) > 0:
                ax, legend_handles = plot_lines_and_points(ax,legend_handles,sector,sector_dataframe=edges,layer_key="edge")
            if len(nodes) > 0:
                ax, legend_handles = plot_lines_and_points(ax,legend_handles,sector,sector_dataframe=nodes,layer_key="node")
            ax.legend(handles=legend_handles,fontsize=12,loc='lower left') 
            save_fig(os.path.join(figures_data_path, f"{sector['sector_label'].replace(' ','_').lower()}_network.png"))
            


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
