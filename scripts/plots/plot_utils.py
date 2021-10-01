"""Functions for plotting
"""
import os
import json
import warnings

from collections import namedtuple

import cartopy.crs as ccrs
import geopandas
import matplotlib.pyplot as plt


def _get_palette():
    colors = {
        'TRANSPARENT': '#00000000',
        'WHITE': '#ffffff',
        'GREY_1': '#e4e4e3',
        'GREY_2': '#dededc',
        'BACKGROUND': '#bfc0bf',
    }
    Palette = namedtuple('Palette', colors.keys())
    return Palette(**colors)


Palette = _get_palette()
JAMAICA_GRID_EPSG = 3448


def save_fig(output_filename):
    plt.savefig(output_filename)


def get_axes(extent=None):
    """Get map axes

    Parameters
    ----------
    extent, optional: tuple (x0, x1, y0, y1)
        to be provided in Jamaica grid coordinates
    """
    ax_proj = ccrs.epsg(JAMAICA_GRID_EPSG)

    plt.figure(figsize=(12, 8), dpi=150)
    ax = plt.axes([0.025, 0.025, 0.95, 0.95], projection=ax_proj)
    if extent is not None:
        ax.set_extent(extent, crs=ax_proj)
    ax.patch.set_facecolor(Palette.BACKGROUND)

    return ax


def plot_basemap(ax, data_path, ax_crs=JAMAICA_GRID_EPSG, plot_regions=False):
    """Plot countries and regions background
    """
    boundaries = os.path.join(data_path, 'boundaries', 'admin_boundaries.gpkg')
    states = geopandas.read_file(boundaries, layer='admin0')#.to_crs(ax_crs)

    states.plot(ax=ax, edgecolor=Palette.WHITE, facecolor=Palette.GREY_1, zorder=1)

    if plot_regions:
        regions = geopandas.read_file(boundaries, layer='admin1').to_crs(ax_crs)
        regions.plot(ax=ax, edgecolor=Palette.TRANSPARENT, facecolor=Palette.GREY_2)
        regions.plot(ax=ax, edgecolor=Palette.WHITE, facecolor=Palette.TRANSPARENT, zorder=2)


def test_plot(data_path, figures_path):
    ax = get_axes(extent=(598251, 838079, 610353, 714779))
    plot_basemap(ax, data_path, plot_regions=True)
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
    plt.savefig(output_filename)


if __name__ == '__main__':
    # Ignore reading-geopackage warnings
    warnings.filterwarnings('ignore', message='.*Sequential read of iterator was interrupted.*')
    # Load config
    CONFIG = load_config()
    test_plot(CONFIG['paths']['data'], CONFIG['paths']['figures'])
    # Show for ease of check/test
    plt.show()
