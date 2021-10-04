"""Functions for plotting
"""
import os
import json
import warnings

from collections import namedtuple

import cartopy.crs as ccrs
import geopandas
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scalebar import scale_bar


def _get_palette():
    colors = {
        'TRANSPARENT': '#00000000',
        'WHITE': '#ffffff',
        'GREY_1': '#e4e4e3',
        'GREY_2': '#dededc',
        # 'BACKGROUND': '#6baed6',
        'BACKGROUND': '#9ecae1',
    }
    Palette = namedtuple('Palette', colors.keys())
    return Palette(**colors)


Palette = _get_palette()
JAMAICA_GRID_EPSG = 3448

def within_extent(x, y, extent):
    """Test x, y coordinates against (xmin, xmax, ymin, ymax) extent
    """
    xmin, xmax, ymin, ymax = extent
    return (xmin < x) and (x < xmax) and (ymin < y) and (y < ymax)

# def get_axes(extent=None):
#     """Get map axes

#     Parameters
#     ----------
#     extent, optional: tuple (x0, x1, y0, y1)
#         to be provided in Jamaica grid coordinates
#     """
#     ax_proj = ccrs.epsg(JAMAICA_GRID_EPSG)

#     plt.figure(figsize=(12, 8), dpi=500)
#     ax = plt.axes([0.025, 0.025, 0.95, 0.95], projection=ax_proj)
#     if extent is not None:
#         ax.set_extent(extent, crs=ax_proj)
#     ax.patch.set_facecolor(Palette.BACKGROUND)

#     return ax

def get_axes(ax,extent=None):
    """Get map axes

    Parameters
    ----------
    extent, optional: tuple (x0, x1, y0, y1)
        to be provided in Jamaica grid coordinates
    """
    ax_proj = ccrs.epsg(JAMAICA_GRID_EPSG)

    # plt.figure(figsize=(12, 8), dpi=500)
    # ax = plt.axes([0.025, 0.025, 0.95, 0.95], projection=ax_proj)
    if extent is not None:
        ax.set_extent(extent, crs=ax_proj)
    ax.patch.set_facecolor(Palette.BACKGROUND)

    return ax

def scale_bar_and_direction(ax,ax_crs=None,length=100,location=(0.5, 0.05), linewidth=3):
    """Draw a scale bar and direction arrow

    Parameters
    ----------
    ax : axes
    length : int
        length of the scalebar in km.
    ax_crs: projection system of the axis
        to be provided in Jamaica grid coordinates
    location: tuple
        center of the scalebar in axis coordinates (ie. 0.5 is the middle of the plot)
    linewidth: float
        thickness of the scalebar.
    """
    # lat-lon limits
    scale_bar(ax, (0.71, 0.05), 25, color='k')

    # Make the direction arrow specific for Jamaica
    if ax_crs is not None:
        x0, x1, y0, y1 = ax.get_extent()
        tmc = ccrs.epsg(ax_crs)

    else:
        llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
        # Transverse mercator for length
        x = (llx1 + llx0) / 2
        y = lly0 + (lly1 - lly0) * location[1]
        tmc = ccrs.TransverseMercator(x, y)

        # Extent of the plotted area in coordinates in metres
        x0, x1, y0, y1 = ax.get_extent(tmc)

    # Scalebar location coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]
    bar_xs = [sbx - length * 10, sbx - length * 90]

    ax.text(x=sbx - length * 150, y=sby + 50*length, s='N', fontsize=14)
    ax.arrow(sbx - length * 135, sby - 5*length, 0, 50*length, length_includes_head=True,
          head_width=50*length, head_length=60*length, overhang=0.2, facecolor='k')

def plot_basemap_labels(ax,ax_crs=None,labels=None,label_column=None,label_offset=0,include_zorder=20):
    """Plot countries and regions background
    """
    if ax_crs is None:
        proj = ccrs.PlateCarree()
    else:
        proj = ccrs.epsg(ax_crs)
    extent = ax.get_extent()
    if labels is not None:
        for label in labels.itertuples():
            text = getattr(label,label_column)
            x = float(label.geometry.centroid.x)
            y = float(label.geometry.centroid.y)
            size = 8
            if within_extent(x, y, extent):
                ax.text(
                    x - 10*label_offset, y - 10*label_offset,
                    text,
                    alpha=0.7,
                    size=size,
                    horizontalalignment='center',
                    zorder = include_zorder,
                    transform=proj)

def plot_basemap(ax, data_path, ax_crs=JAMAICA_GRID_EPSG, plot_regions=False, region_labels=False):
    """Plot countries and regions background
    """
    boundaries = os.path.join(data_path, 'boundaries', 'admin_boundaries.gpkg')
    states = geopandas.read_file(boundaries, layer='admin0')#.to_crs(ax_crs)

    states.plot(ax=ax, edgecolor=Palette.WHITE, facecolor=Palette.GREY_1, zorder=1)

    if plot_regions:
        regions = geopandas.read_file(boundaries, layer='admin1').to_crs(ax_crs)
        regions.plot(ax=ax, edgecolor=Palette.TRANSPARENT, facecolor=Palette.GREY_2)
        regions.plot(ax=ax, edgecolor=Palette.WHITE, facecolor=Palette.TRANSPARENT, zorder=2)
        if region_labels is True:
            plot_basemap_labels(ax,ax_crs=ax_crs,
                                labels=regions,label_column='PARISH',label_offset=100)
    scale_bar_and_direction(ax,ax_crs,location=(0.75, 0.05))

def plot_point_assets(ax,ax_crs,nodes,colors,size,marker,zorder):
    proj_lat_lon = ccrs.epsg(ax_crs)
    ax.scatter(
        list(nodes.geometry.x),
        list(nodes.geometry.y),
        transform=proj_lat_lon,
        facecolor=colors,
        s=size,
        marker=marker,
        zorder=zorder
    )
    return ax

def plot_line_assets(ax,ax_crs,edges,colors,size,zorder):
    proj_lat_lon = ccrs.epsg(ax_crs)
    ax.add_geometries(
        list(edges.geometry),
        crs=proj_lat_lon,
        linewidth=size,
        edgecolor=colors,
        facecolor='none',
        zorder=zorder
    )
    return ax

def get_sector_layer(sector,sector_data_path,layer_key):
    layer_name_key = f"{layer_key}_layer"
    layer_filter_key = f"{layer_key}_filter_column"
    layer_filter_values = f"{layer_key}_categories"
    if sector[layer_name_key] is not None:
        layer = geopandas.read_file(os.path.join(sector_data_path,sector["sector_gpkg"]),layer=sector[layer_name_key])
        layer = layer[layer[sector[layer_filter_key]].isin(sector[layer_filter_values])]
        return layer
    else:
        return []

def plot_lines_and_points(ax,legend_handles,sector,sector_dataframe=None,layer_key=None):  
    layer_details = list(zip(sector[f"{layer_key}_categories"],
                                        sector[f"{layer_key}_categories_colors"],
                                        sector[f"{layer_key}_categories_labels"],
                                        sector[f"{layer_key}_categories_zorder"]))
    use_labels = []
    for i,(cat,color,label,zorder) in enumerate(layer_details):
        if layer_key == "edge":
            ax = plot_line_assets(ax,JAMAICA_GRID_EPSG,
                                sector_dataframe[sector_dataframe[sector["edge_filter_column"]] == cat],
                                color,
                                sector["edge_categories_linewidth"][i],
                                zorder)
            if label not in use_labels:
                legend_handles.append(mpatches.Patch(color=color,
                                                label=label))
                use_labels.append(label)
        elif layer_key == "node":
            ax = plot_point_assets(ax,JAMAICA_GRID_EPSG,
                            sector_dataframe[sector_dataframe[sector["node_filter_column"]] == cat],
                            color,
                            sector["node_categories_markersize"][i],
                            sector["node_categories_marker"][i],
                            zorder)
            if label not in use_labels:
                legend_handles.append(plt.plot([],[],
                                        marker=sector["node_categories_marker"][i], 
                                        ms=sector["node_categories_markersize"][i], 
                                        ls="",
                                        color=color,
                                        label=label)[0])
                use_labels.append(label)
        
    return ax, legend_handles

def test_plot(data_path, figures_path):
    plt.figure(figsize=(12, 8), dpi=500)
    ax = plt.axes([0.025, 0.025, 0.95, 0.95], projection=ccrs.epsg(JAMAICA_GRID_EPSG))
    ax = get_axes(ax,extent=(598251, 838079, 610353, 714779))
    plot_basemap(ax, data_path, plot_regions=True, region_labels=True)
    save_fig(os.path.join(figures_path, "admin_map.png"))


def load_config():
    """Read config.json"""
    config_path = os.path.join(os.path.dirname(__file__), "..", "..", "config.json")
    with open(config_path, "r") as config_fh:
        config = json.load(config_fh)
    return config


def geopandas_read_file_type(file_path, file_layer, file_database=None):
    if file_database is not None:
        return geopandas.read_file(os.path.join(file_path, file_database), layer=file_layer)
    else:
        return geopandas.read_file(os.path.join(file_path, file_layer))


def save_fig(output_filename):
    print(" * Save", os.path.basename(output_filename))
    plt.savefig(output_filename,bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    # Ignore reading-geopackage warnings
    warnings.filterwarnings('ignore', message='.*Sequential read of iterator was interrupted.*')
    # Load config
    CONFIG = load_config()
    test_plot(CONFIG['paths']['data'], CONFIG['paths']['figures'])
    # Show for ease of check/test
    plt.show()
