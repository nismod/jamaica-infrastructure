"""Functions for plotting
"""

import os
import geopandas
import json
import pandas
import numpy


def load_config():
    """Read config.json"""
    config_path = os.path.join(os.path.dirname(__file__), "..", "..", "config.json")
    with open(config_path, "r") as config_fh:
        config = json.load(config_fh)
    return config


def geopandas_read_file_type(file_path, file_layer, file_database=None):
    if file_database is not None:
        return geopandas.read_file(
            os.path.join(file_path, file_database), layer=file_layer
        )
    else:
        return geopandas.read_file(os.path.join(file_path, file_layer))


def get_sector_layer(sector, sector_data_path, layer_key):
    layer_name_key = f"{layer_key}_layer"
    layer_classify_key = f"{layer_key}_classify_column"
    layer_classify_values = f"{layer_key}_categories"
    if sector[layer_name_key] is not None:
        layer = geopandas.read_file(
            os.path.join(sector_data_path, sector["sector_gpkg"]),
            layer=sector[layer_name_key],
        )
        layer = layer[
            layer[sector[layer_classify_key]].isin(sector[layer_classify_values])
        ]
        return layer
    else:
        return []
