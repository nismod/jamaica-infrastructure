"""Estimate direct damages to physical assets exposed to hazards

"""

import sys
import os
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import geopandas as gpd
from shapely import wkb
import numpy as np
from SALib.sample import morris
import SALib.analyze.morris

from analysis_utils import *
from tqdm import tqdm

tqdm.pandas()


def select_damage_curves(
    asset_df, asset_type_column, damage_curve_lookup_df, assigned_sector, hazard_type
):
    selected_assets = damage_curve_lookup_df[
        (damage_curve_lookup_df["sector"] == assigned_sector)
        & (damage_curve_lookup_df["hazard_type"] == hazard_type)
    ][["asset_name", "asset_sheet"]]
    data = pd.read_excel(
        os.path.join(damage_data_path, f"damage_curves_{sector}_{hazard_type}.xlsx"),
        sheet_name=data_key.asset_sheet,
    )
    asset_df = pd.merge(
        asset_df,
        selected_assets,
        how="left",
        left_on=asset_type_column,
        right_on="asset_name",
    )
    asset_df = asset_df[~asset_df["asset_sheet"].isna()]


def get_damage_data(x, damage_data_path, uplift_factor=0, uncertainty_parameter=0):
    data = pd.read_excel(
        os.path.join(
            damage_data_path, f"damage_curves_{x.sector}_{x.hazard_type}.xlsx"
        ),
        sheet_name=x.asset_sheet,
    )
    if x.hazard_type == "flooding":
        x_data = data.flood_depth
    else:
        x_data = data.wind_speed

    y_data = data.damage_ratio_min + uncertainty_parameter * (
        data.damage_ratio_max - data.damage_ratio_min
    )
    # print(y_data)
    y_data = np.minimum(y_data + uplift_factor, 1.0)

    return x_data.values, y_data.values


def convert_cost_units(x, cost_value, cost_unit, conversion_rate):
    if ("US$" in x[cost_unit]) or ("USD" in x[cost_unit]):
        return conversion_rate * x[cost_value]
    else:
        return x[cost_value]


def modify_cost_units(x, cost_dimension, damage_cost_column="damage_cost"):
    if "/km" in str(x[cost_dimension]):
        return 0.001 * x[damage_cost_column]
    else:
        return x[damage_cost_column]


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


def create_damage_curves(
    damage_data_path, damage_curve_lookup_df, uplift_factor=0, uncertainty_parameter=0
):
    damage_curve_lookup_df["x_y_data"] = damage_curve_lookup_df.progress_apply(
        lambda x: get_damage_data(
            x, damage_data_path, uplift_factor, uncertainty_parameter
        ),
        axis=1,
    )
    damage_curve_lookup_df[["damage_x_data", "damage_y_data"]] = damage_curve_lookup_df[
        "x_y_data"
    ].apply(pd.Series)
    damage_curve_lookup_df.drop("x_y_data", axis=1, inplace=True)

    return damage_curve_lookup_df


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
    set_count = 0
    cost_uncertainty_parameter = 0
    damage_uncertainty_parameter = 1
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    jd_to_usd = 0.0067  # Jamaican dollar to USD conversion
    direct_damages_results = os.path.join(output_data_path, "direct_damages")

    nomangrove = pd.read_parquet(
        os.path.join(
            direct_damages_results,
            "airport_polygon_areas",
            "airport_polygon_areas_direct_damages_parameter_set_0.parquet",
        )
    )
    print(nomangrove)
    coastal_columns = [c for c in nomangrove.columns.values.tolist() if "coastal" in c]
    nomangrove[["node_id", "exposure"] + coastal_columns].to_csv(
        os.path.join(
            direct_damages_results,
            "airport_polygon_areas",
            "airport_polygon_areas_direct_damages_parameter_set_0.csv",
        )
    )

    direct_damages_results = os.path.join(
        output_data_path, "direct_damages_with_mangroves"
    )
    mangrove = pd.read_parquet(
        os.path.join(
            direct_damages_results,
            "airport_polygon_areas",
            "airport_polygon_areas_direct_damages_parameter_set_0.parquet",
        )
    )
    print(mangrove)
    mangrove[["node_id", "exposure"] + coastal_columns].to_csv(
        os.path.join(
            direct_damages_results,
            "airport_polygon_areas",
            "airport_polygon_areas_direct_damages_parameter_set_0.csv",
        )
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
