"""Estimate direct damages to physical assets exposed to hazard events
"""
import logging
import os
import sys
# import warnings
# warnings.simplefilter(action="ignore", category=FutureWarning)

import pandas as pd

# pd.options.mode.chained_assignment = None  # default='warn'
# warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

import geopandas as gpd
import numpy as np

from analysis_utils import *
from tqdm import tqdm

tqdm.pandas()


def create_damage_curves(
    damage_data_path, damage_curve_lookup_df, uplift_factor=0, uncertainty_parameter=0
):
    damage_curve_lookup_df["x_y_data"] = damage_curve_lookup_df.apply(
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

    y_data = np.minimum(y_data + uplift_factor, 1.0)

    return x_data.values, y_data.values


def convert_cost_units(x, cost_value, cost_unit):
    jd_to_usd = 0.0067  # Jamaican dollar to USD conversion
    usd_to_jd = 1 / jd_to_usd

    if ("US$" in x[cost_unit]) or ("USD" in x[cost_unit]):
        return usd_to_jd * x[cost_value]
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


def get_damage_curves(processed_data_path, hazard_attributes, damage_uncertainty_parameter):
    damage_data_path = os.path.join(processed_data_path, "damage_curves")
    damage_curve_lookup = pd.read_csv(
        os.path.join(damage_data_path, "asset_damage_curve_mapping.csv"),
        usecols=["sector", "hazard_type", "asset_name", "asset_sheet"]
    )

    damage_curves = []
    for _, hazard in hazard_attributes.iterrows():
        damage_curve_df = damage_curve_lookup[
            damage_curve_lookup["hazard_type"] == hazard["hazard_type"]
        ].copy()
        damage_curve_df["hazard"] = hazard["hazard"]

        damage_curve_df = create_damage_curves(
            damage_data_path,
            damage_curve_df,
            uplift_factor=hazard["uplift_factor"],
            uncertainty_parameter=damage_uncertainty_parameter,
        )
        damage_curves.append(damage_curve_df)

    return pd.concat(damage_curves, axis=0, ignore_index=True)


def get_asset_data(processed_data_path, cost_uncertainty_parameter, asset):
    asset_min_cost = asset.asset_min_cost_column
    asset_max_cost = asset.asset_max_cost_column
    asset_cost_unit = asset.asset_cost_unit_column

    asset_df = gpd.read_file(
        os.path.join(processed_data_path, asset.path),
        layer=asset.asset_layer,
    )
    asset_df[asset_min_cost] = asset_df.apply(
        lambda x: convert_cost_units(
            x, asset_min_cost, asset_cost_unit
        ),
        axis=1,
    )
    asset_df[asset_max_cost] = asset_df.apply(
        lambda x: convert_cost_units(
            x, asset_max_cost, asset_cost_unit
        ),
        axis=1,
    )
    asset_df[asset_cost_unit] = asset_df[asset_cost_unit].replace(
        ["USD", "US$"], "J$", regex=True
    )

    # We need to modify and correct energy asset costs of nodes to J$ and costs of of edges to $J/m
    if asset.sector == "energy" and asset.asset_layer == "edges":
        asset_df[asset_min_cost] = asset_df[asset_min_cost] / asset_df["length"]
        asset_df[asset_max_cost] = asset_df[asset_max_cost] / asset_df["length"]
        asset_df[asset_cost_unit] = "J$/m"

    asset_df["damage_cost"] = (
        asset_df[asset_min_cost] +
        cost_uncertainty_parameter * (asset_df[asset_max_cost] - asset_df[asset_min_cost])
    )
    asset_df["damage_cost"] = asset_df.progress_apply(
        lambda x: modify_cost_units(x, asset_cost_unit), axis=1
    )

    return asset_df


def main(config, set_count, cost_uncertainty_parameter, damage_uncertainty_parameter):
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    event_damages_results_dir = os.path.join(output_data_path, "direct_damages_events")
    if os.path.exists(event_damages_results_dir) == False:
        os.mkdir(event_damages_results_dir)

    hazard_attributes = pd.DataFrame([
        {
            "hazard": "fluvial",
            "hazard_type": "flooding",
            "hazard_threshold": 0.5,
            "uplift_factor": 0.0,
        },
        {
            "hazard": "surface",
            "hazard_type": "flooding",
            "hazard_threshold": 0.5,
            "uplift_factor": 0.0,
        },
    ])

    # DataFrame with columns
    # 'sector', 'hazard_type', 'asset_name', 'asset_sheet', 'hazard', 'damage_x_data', 'damage_y_data'
    # 'damage_{x,y}_data' have lists in each cell value, defining the damage curves
    damage_curves = get_damage_curves(processed_data_path, hazard_attributes, damage_uncertainty_parameter)

    splits_dir = os.path.join(
        output_data_path, "hazard_asset_intersection"
    )
    # DataFrame with columns
    # 'sector', 'asset_description', 'asset_gpkg', 'asset_layer', 'asset_id_column',
    # 'asset_min_cost_column', 'asset_max_cost_column', 'asset_mean_cost_column',
    # 'asset_reopen_cost_column', 'asset_cost_unit_column', 'asset_reopen_cost_unit_column',
    # 'flooding_asset_damage_lookup_column', 'TC_asset_damage_lookup_column', 'coastal_asset_damage_lookup_column',
    # 'fluvial_asset_damage_lookup_column', 'surface_asset_damage_lookup_column', 'cyclone_asset_damage_lookup_column',
    # 'path', 'single_failure_scenarios'
    network_hazard_intersections = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks",
            "network_layers_hazard_intersections_details.csv",
        )
    )

    # TODO pick from command-line parameters
    hazard = hazard_attributes[0]
    asset = network_hazard_intersections[0]
    breakpoint()

    logging.info("Processing %s")
    assets = get_asset_data(processed_data_path, cost_uncertainty_parameter, asset)

    splits_fname = os.path.join(
        splits_dir,
        f"{asset.asset_gpkg}_splits__hazard_layers__{asset.asset_layer}.geoparquet",
    )
    if not os.path.isfile(splits_fname):
        logging.info("No splits file at %s", splits_fname)
        return

    splits = gpd.read_parquet(splits_fname)
    splits = splits.to_crs(epsg=epsg_jamaica)
    splits = add_exposure_dimensions(
        splits, dataframe_type=asset.asset_layer, epsg=epsg_jamaica
    )
    id_col = asset.asset_id_column
    cost_unit_col = asset.asset_cost_unit_column

    asset_hazard_col = getattr(asset, f"{hazard.hazard}_asset_damage_lookup_column")

    if (asset_hazard_col == "none"):
        logging.info("%s %s not affected by %s", asset.asset_gpkg, asset.asset_layer, hazard.hazard)
        return

    # Get damage curves for this sector, hazard
    asset_hazard_damage_curves = damage_curves[
        (damage_curves["sector"] == asset.sector)
        & (damage_curves["hazard"] == hazard.hazard)
    ]
    # Find the set of asset types with damage curves
    asset_types = set(assets[asset_hazard_col].unique()) & set(asset_hazard_damage_curves.asset_name.unique())

    asset_cols = [id_col, asset_hazard_col, cost_unit_col, "damage_cost"]
    hazard_cols = None
    # TODO hazard_cols = list of events
    exposure = None
    # TODO exposure = splits join event depths
    # TODO exposure.exposure_factor = 1 or length or area (to multiply by damage_cost for damage)
    # TODO Filter splits to only relevant assets (by ID, by type?)
    # TODO Subtract `hazard.hazard_threshold` from event depth
    # TODO Calculate damages at the level of splits

    for curve in asset_hazard_damage_curves.itertuples():
        hazard_intensity, damage_fraction = curve.damage_x_data, curve.damage_y_data
        bounds = tuple(f(damage_fraction) for f in (min, max))
        interpolated_damage_curve = interp1d(
            hazard_intensity,
            damage_fraction,
            kind='linear',
            fill_value=bounds,
            bounds_error=False
        )

        # pick out rows of asset type and columns of hazard intensity
        # TODO check col to use for .asset_type
        asset_type_mask: gpd.GeoDataFrame = exposure.asset_type == curve.asset_name
        asset_exposure: pd.DataFrame = pd.DataFrame(exposure.loc[asset_type_mask, hazard_cols])

        # apply damage_curve function to exposure table
        # the return value of interpolated_damage_curve is a numpy array
        logging.info("Calculating damage fractions")
        damage_fraction_for_asset_type = pd.DataFrame(
            interpolated_damage_curve(asset_exposure),
            index=asset_exposure.index,
            columns=asset_exposure.columns
        )

        # TODO put asset_cols back on damage_fraction

        direct_damages_only = damage_fraction_for_asset_type[hazard_cols] \
            .multiply(damage_fraction_for_asset_type.damage_cost, axis="index") \
            .multiply(damage_fraction_for_asset_type.exposure_factor, axis="index")

        # TODO save to file - pick a good name
        direct_damages_only.to_parquet(out_fname)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    CONFIG = load_config()
    try:
        SET_COUNT = str(sys.argv[1])
        COST_UNCERTAINTY_PARAMETER = float(sys.argv[2])
        DAMAGE_UNCERTAINTY_PARAMETER = float(sys.argv[3])
        # TODO asset_gpkg, asset_layer, hazard
    except IndexError:
        logging.warning("""Expected to be called with
    ['direct_damages_calculations_events.py', set_count, cost_uncertainty, damage_uncertainty]
Got arguments
    %s""", sys.argv)
        exit()

    main(CONFIG, SET_COUNT, COST_UNCERTAINTY_PARAMETER, DAMAGE_UNCERTAINTY_PARAMETER)
