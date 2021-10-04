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


    transport_data_path = os.path.join(processed_data_path,
                        "networks",
                        "transport")
    sector_attributes = jamaica_sector_attributes()
    for sector in sector_attributes:
        if sector["sector_label"] in ["Roads","Railways"]:
            edges = get_sector_layer(sector,transport_data_path,"edge")
            nodes = get_sector_layer(sector,transport_data_path,"node")
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
            ax.legend(handles=legend_handles,fontsize=10,loc='lower left') 
            save_fig(os.path.join(figures_data_path, f"{sector['sector_label'].lower().replace(' ','_')}_network.png"))
    

    """Plot ports and airports
    """
    ports = gpd.read_file(os.path.join(transport_data_path,"port_polygon.gpkg"),layer="areas").to_crs(epsg=JAMAICA_GRID_EPSG)
    ports["asset_type"] = "port"
    airports = gpd.read_file(os.path.join(transport_data_path,"airport_polygon.gpkg"),layer="areas").to_crs(epsg=JAMAICA_GRID_EPSG)
    airports["asset_type"] = "airport"
    nodes = pd.concat([ports[["name","asset_type","geometry"]],
                        airports[["name","asset_type","geometry"]]],
                        axis=0,ignore_index=True)
    
    nodes["centroid"] = nodes.geometry.centroid
    nodes.drop("geometry",axis=1,inplace=True)
    nodes.rename(columns={"centroid":"geometry"},inplace=True)
    nodes = gpd.GeoDataFrame(nodes,geometry="geometry",crs=f"EPSG:{JAMAICA_GRID_EPSG}")

    nodes_categories_colors = {"port":"#08306b", 
                                "airport":"#8c510a", 
                                }

    legend_handles = []
    fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(12,8),
                dpi=500)
    ax = get_axes(ax,extent=JAMAICA_EXTENT)
    plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=False)
    for node_type,node_color in nodes_categories_colors.items():
        ax = plot_point_assets(ax,JAMAICA_GRID_EPSG,nodes[nodes["asset_type"] == node_type],node_color,15.0,'o',11)
        legend_handles.append(plt.plot([],[],
                                marker='.', 
                                ms=15.0, 
                                ls="",
                                color=node_color,
                                label=f"{node_type.upper()} AREAS")[0])
    labels = nodes.drop_duplicates("name",keep = "first")
    used_names = []
    for l in labels.itertuples():
        if "kingston" in l.name.lower() or l.name in ["Petrojam","Wherry Wharf","Tinson Pen Aerodrome"]:
            name = "Kingston Port and Aerodrome"
        else:
            name = l.name
        if name not in used_names:
            if name == "Kingston Port and Aerodrome":
                location_x = l.geometry.x - 10000
                location_y = l.geometry.y + 2000
            elif name == "Falmouth":
                location_x = l.geometry.x - 1000
                location_y = l.geometry.y + 1500
            elif name == "Port Rhoades":
                location_x = l.geometry.x
                location_y = l.geometry.y - 2000
            else:
                location_x = l.geometry.x - 10000
                location_y = l.geometry.y + 1000

            ax.text(location_x,location_y,name,size=6,weight="bold")
            used_names.append(name)

    ax.legend(handles=legend_handles,fontsize=10,loc='lower left') 
    save_fig(os.path.join(figures_data_path, "ports_airports.png"))


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
