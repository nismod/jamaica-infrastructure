"""Generate hazard-damage curves
"""
import os
import sys
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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

    boundaries = gpd.read_file(os.path.join(processed_data_path,
                                                "boundaries",
                                                "admin_boundaries.gpkg"),
                                layer='admin1').to_crs(JAMAICA_GRID_EPSG)
    boundaries= boundaries[boundaries["PARISH"].isin(["Kingston","St. Andrew"])]
    bounds = boundaries.geometry.total_bounds # this gives your boundaries of the map as (xmin,ymin,xmax,ymax)
    kst_extent = (bounds[0]-10000,bounds[2]+10000,bounds[1],bounds[3])

    sector_attributes = jamaica_sector_attributes()


    """Plot of hazard maps
    """
    flood_file = os.path.join(processed_data_path,
                            "hazards",
                            "Global Flood Map",
                            "Jamaica",
                            "Fluvial",
                            "Raw Depths",
                            "JM_FLRF_UD_Q100_RD_02.tif")
    fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
    ax = get_axes(ax,extent=JAMAICA_EXTENT)
    plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=False)
    scale_bar_and_direction(ax,scalebar_distance=25)

    proj_lat_lon = ccrs.epsg(JAMAICA_GRID_EPSG)
    im = plot_raster(ax, flood_file,cmap="terrain",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=None)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
    # cbar.set_clim(vmin=0,vmax=max_val)

    cbar.outline.set_color("none")
    cbar.ax.yaxis.set_tick_params(color='black')
    cbar.ax.set_xlabel('Flood depths (m)',fontsize=7,color='black')

    plt.title("Baseline - 1 in 100 year river flooding", fontsize = 10)
    save_fig(os.path.join(figures_data_path,"jamaica_river_flood_map.png"))

    tc_file = os.path.join(processed_data_path,
                            "hazards",
                            "TC_data_fixed_return",
                            "STORM_FIXED_RETURN_PERIODS_jamaica_baseline_rp100_mean.tif"
                            )

    fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
    ax = get_axes(ax,extent=JAMAICA_EXTENT)
    plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=False)
    scale_bar_and_direction(ax,scalebar_distance=25)

    proj_lat_lon = ccrs.epsg(JAMAICA_GRID_EPSG)
    im = plot_raster(ax, tc_file,cmap="cool",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=None)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax,fraction=0.05, shrink=0.7,pad=0.01, drawedges=False, orientation='horizontal')
    # cbar.set_clim(vmin=0,vmax=max_val)

    cbar.outline.set_color("none")
    cbar.ax.yaxis.set_tick_params(color='black')
    cbar.ax.set_xlabel('Wind speeds (m/s)',fontsize=7,color='black')

    plt.title("Baseline - 1 in 100 year Tropical Cyclone", fontsize = 10)
    save_fig(os.path.join(figures_data_path,"jamaica_tc_map.png"))

    coastal_file = os.path.join(processed_data_path,
                                "hazards",
                                "Coastal_flood_data",
                                "Flood_maps_current_climate",
                                "JamaicaJAM001RCP452010_epsg_32618_RP_100.tif"
                            )
    
    fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
    ax = get_axes(ax,extent=JAMAICA_EXTENT)
    plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
    scale_bar_and_direction(ax,scalebar_distance=25)

    proj_lat_lon = ccrs.epsg(JAMAICA_GRID_EPSG)
    im = plot_raster(ax, coastal_file,cmap="coolwarm",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=None)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
    # cbar.set_clim(vmin=0,vmax=max_val)

    cbar.outline.set_color("none")
    cbar.ax.yaxis.set_tick_params(color='black')
    cbar.ax.set_xlabel('Flood depth (m)',fontsize=7,color='black')

    plt.title("Baseline - 1 in 100 year Coastal flooding", fontsize = 10)
    save_fig(os.path.join(figures_data_path,"jamaica_coastal_flood_map.png"))

    for sector in sector_attributes:
        asset_data_path = os.path.join(processed_data_path,
                        "networks",
                        sector["sector"])
        if sector["sector_label"] in ["Roads"]:
            edges = get_sector_layer(sector,asset_data_path,"edge")
            nodes = get_sector_layer(sector,asset_data_path,"node")
            # print (edges, nodes)

            legend_handles = []
            fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
            ax = get_axes(ax,extent=kst_extent)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
            scale_bar_and_direction(ax,scalebar_distance=5)
            color = "#252525"
            if len(edges) > 0:
                ax = plot_line_assets(ax,JAMAICA_GRID_EPSG,edges,color,0.6,10)
                legend_handles.append(mpatches.Patch(color=color,
                                                label="ROADS"))
            if len(nodes) > 0:
                ax = plot_point_assets(ax,JAMAICA_GRID_EPSG,nodes,color,10,"o",11)
                legend_handles.append(plt.plot([],[],
                                        marker="o", 
                                        ms=10, 
                                        ls="",
                                        color=color,
                                        label="BRIDGES")[0])

            ax.legend(handles=legend_handles,fontsize=7,loc='lower left',).set_zorder(20) 

            im = plot_raster(ax, flood_file,cmap="terrain",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=kst_extent)

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
            # cbar.set_clim(vmin=0,vmax=max_val)

            cbar.outline.set_color("none")
            cbar.ax.yaxis.set_tick_params(color='black')
            cbar.ax.set_xlabel('Flood depths (m)',fontsize=7,color='black')

            plt.title(f"{sector['sector_label']} 1 in 100 year River Flooding exposure", fontsize = 10)
            save_fig(os.path.join(figures_data_path, f"kst_{sector['sector_label'].replace(' ','_').lower()}_network_river_flood_exposure.png"))

            
            legend_handles = []
            fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
            ax = get_axes(ax,extent=kst_extent)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
            scale_bar_and_direction(ax,scalebar_distance=5)

            color = "#252525"
            if len(edges) > 0:
                ax = plot_line_assets(ax,JAMAICA_GRID_EPSG,edges,color,0.6,10)
                legend_handles.append(mpatches.Patch(color=color,
                                                label="ROADS"))
            if len(nodes) > 0:
                ax = plot_point_assets(ax,JAMAICA_GRID_EPSG,nodes,color,10,"o",11)
                legend_handles.append(plt.plot([],[],
                                        marker="o", 
                                        ms=10, 
                                        ls="",
                                        color=color,
                                        label="BRIDGES")[0])
            ax.legend(handles=legend_handles,fontsize=7,loc='lower left',).set_zorder(20) 
            im = plot_raster(ax, coastal_file,cmap="coolwarm",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=None)

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
            # cbar.set_clim(vmin=0,vmax=max_val)

            cbar.outline.set_color("none")
            cbar.ax.yaxis.set_tick_params(color='black')
            cbar.ax.set_xlabel('Flood depth (m)',fontsize=7,color='black')
            plt.title(f"{sector['sector_label']} 1 in 100 year Coastal Flooding exposure", fontsize = 10)
            save_fig(os.path.join(figures_data_path, 
                            f"kst_{sector['sector_label'].replace(' ','_').lower()}_network_coastal_flood_exposure.png"))

            # legend_handles = []
            fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
            ax = get_axes(ax,extent=kst_extent)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
            scale_bar_and_direction(ax,scalebar_distance=5)

            color = "#252525"
            if len(edges) > 0:
                ax = line_map_plotting_colors_width(
                                                        ax,edges,"traffic_count",
                                                        legend_label="Traffic count (vehicles/day)",
                                                        no_value_label="No data",
                                                        width_step=10,
                                                        interpolation="log",
                                                        plot_title=None
                                                        )
            im = plot_raster(ax, coastal_file,cmap="coolwarm",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=None)

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
            # cbar.set_clim(vmin=0,vmax=max_val)

            cbar.outline.set_color("none")
            cbar.ax.yaxis.set_tick_params(color='black')
            cbar.ax.set_xlabel('Flood depth (m)',fontsize=7,color='black')
            plt.title(f"{sector['sector_label']} traffic counts within 1 in 100 year Coastal Flooding exposure", fontsize = 10)
            save_fig(os.path.join(figures_data_path, 
                            f"kst_{sector['sector_label'].replace(' ','_').lower()}_network_coastal_flood_exposure_with_traffic_counts.png"))

            # legend_handles = []
            fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
            ax = get_axes(ax,extent=kst_extent)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
            scale_bar_and_direction(ax,scalebar_distance=5)

            color = "#252525"
            if len(edges) > 0:
                ax = line_map_plotting_colors_width(
                                                        ax,edges,"traffic_count",
                                                        legend_label="Traffic count (vehicles/day)",
                                                        no_value_label="No data",
                                                        width_step=10,
                                                        interpolation="log",
                                                        plot_title=None
                                                        )
            im = plot_raster(ax, flood_file,cmap="coolwarm",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=None)

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
            # cbar.set_clim(vmin=0,vmax=max_val)

            cbar.outline.set_color("none")
            cbar.ax.yaxis.set_tick_params(color='black')
            cbar.ax.set_xlabel('Flood depth (m)',fontsize=7,color='black')
            plt.title(f"{sector['sector_label']} traffic counts within 1 in 100 year River Flooding exposure", fontsize = 10)
            save_fig(os.path.join(figures_data_path, 
                            f"kst_{sector['sector_label'].replace(' ','_').lower()}_network_river_flood_exposure_with_traffic_counts.png"))

        if sector["sector_label"] in ["Energy"]:
            edges = get_sector_layer(sector,asset_data_path,"edge")
            nodes = get_sector_layer(sector,asset_data_path,"node")
            # print (edges, nodes)

            legend_handles = []
            fig, ax = plt.subplots(1,1,
                subplot_kw={'projection': ccrs.epsg(JAMAICA_GRID_EPSG)},
                figsize=(8,6),
                dpi=500)
            ax = get_axes(ax,extent=kst_extent)
            plot_basemap(ax, processed_data_path, plot_regions=True, region_labels=True)
            scale_bar_and_direction(ax,scalebar_distance=5)

            if len(edges) > 0:
                ax, legend_handles = plot_lines_and_points(ax,
                                        legend_handles,sector,
                                        sector_dataframe=edges,
                                        layer_key="edge",
                                        marker_size_factor=0.5)
            if len(nodes) > 0:
                ax, legend_handles = plot_lines_and_points(ax,
                                        legend_handles,sector,
                                        sector_dataframe=nodes,
                                        layer_key="node",
                                        marker_size_factor=0.5)
            ax.legend(handles=legend_handles,fontsize=7,loc='lower left',).set_zorder(20) 

            im = plot_raster(ax, tc_file,cmap="cool",
                reproject_transform=JAMAICA_GRID_EPSG,
                clip_extent=kst_extent)

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal')
            # cbar.set_clim(vmin=0,vmax=max_val)

            cbar.outline.set_color("none")
            cbar.ax.yaxis.set_tick_params(color='black')
            cbar.ax.set_xlabel('Wind Speeds (m/s)',fontsize=7,color='black')

            plt.title(f"{sector['sector_label']} 1 in 100 year Tropical Cyclone exposure", fontsize = 10)
            save_fig(os.path.join(figures_data_path, f"kst_{sector['sector_label'].replace(' ','_').lower()}_network_tc_exposure.png"))



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
