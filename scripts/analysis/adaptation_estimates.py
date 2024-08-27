"""Estimate adaptation options costs and benefits

"""

import sys
import os

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import geopandas as gpd
import numpy as np
from analysis_utils import *
from tqdm import tqdm

tqdm.pandas()


def add_exposure_dimensions(dataframe, dataframe_type="nodes", epsg=4326):
    geo_dataframe = gpd.GeoDataFrame(
        dataframe, geometry="geometry", crs={"init": f"epsg:{epsg}"}
    )
    if dataframe_type == "edges":
        geo_dataframe["exposure"] = geo_dataframe.geometry.length
        geo_dataframe["exposure_unit"] = "m"
    elif dataframe_type == "areas":
        geo_dataframe["exposure"] = geo_dataframe.geometry.area
        geo_dataframe["exposure_unit"] = "m2"
    else:
        geo_dataframe["exposure"] = 1
        geo_dataframe["exposure_unit"] = "unit"
    geo_dataframe.drop("geometry", axis=1, inplace=True)

    index_columns = [
        c for c in geo_dataframe.columns.values.tolist() if c != "exposure"
    ]
    return (
        geo_dataframe.groupby(index_columns, dropna=False)["exposure"]
        .sum()
        .reset_index()
    )


def estimate_direct_damage_costs_and_units(
    dataframe,
    damage_ratio_columns,
    cost_unit_column,
    damage_cost_column="damage_cost",
    dataframe_type="nodes",
):
    if dataframe_type == "nodes":
        dataframe[damage_ratio_columns] = dataframe[damage_ratio_columns].multiply(
            dataframe[damage_cost_column], axis="index"
        )
        dataframe["damage_cost_unit"] = dataframe[cost_unit_column]
    else:
        dataframe[damage_ratio_columns] = dataframe[damage_ratio_columns].multiply(
            dataframe[damage_cost_column] * dataframe["exposure"], axis="index"
        )
        cost_unit = dataframe[cost_unit_column].values.tolist()[0]
        dataframe["damage_cost_unit"] = "/".join(cost_unit.split("/")[:-1])

    return dataframe


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    baseline_year = 2019
    projection_end_year = 2100
    discounting_rate = 10

    asset_data_details = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks",
            "network_layers_hazard_intersections_details.csv",
        )
    )
    asset_data_details = asset_data_details[asset_data_details["sector"] != "buildings"]
    hazard_asset_intersection_path = os.path.join(
        output_data_path, "hazard_asset_intersection"
    )

    flood_hazards = ["coastal", "fluvial", "surface"]
    flood_threshold = 0.5
    hazard_data_details = pd.read_csv(
        os.path.join(processed_data_path, "hazards", "hazard_layers.csv"),
        encoding="latin1",
    )
    hazard_keys = hazard_data_details[
        hazard_data_details["hazard"].isin(flood_hazards)
    ].key.values.tolist()

    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column
        index_columns = [asset_id, "damage_cost_unit", "hazard"]

        hazard_intersection_file = os.path.join(
            hazard_asset_intersection_path,
            f"{asset_info.asset_gpkg}_splits__hazard_layers__{asset_info.asset_layer}.geoparquet",
        )
        if os.path.isfile(hazard_intersection_file) is True:
            hazard_df = gpd.read_parquet(hazard_intersection_file)
            hazard_df = hazard_df.to_crs(epsg=epsg_jamaica)
            hazard_df = add_exposure_dimensions(
                hazard_df, dataframe_type=asset_info.asset_layer, epsg=epsg_jamaica
            )
            hazard_df = hazard_df[[asset_id, "exposure", "exposure_unit"] + hazard_keys]
            hazard_df["max_flood_depth"] = hazard_df[hazard_keys].max(axis=1)
            hazard_df = hazard_df[hazard_df["max_flood_depth"] > flood_threshold]

            print(hazard_df.sort_values(by="max_flood_depth", ascending=False).head(20))
            print(f"* Done with {asset_info.asset_gpkg}")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
