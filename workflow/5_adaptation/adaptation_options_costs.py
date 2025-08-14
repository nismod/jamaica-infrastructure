import os
import logging
import warnings

import click
import pandas as pd
import geopandas as gpd
import numpy as np
from tqdm import tqdm

from jamaica_infrastructure.adaptation import (
    calculate_discounting_rate_factor,
)

warnings.simplefilter(action="ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

tqdm.pandas()


def assign_maintenance_cost_over_time(
    df,
    maintenance_intervals_years,
    maintenance_cost,
    start_year=2019,
    end_year=2100
):
    maintain_intervals = list(
        set(df[maintenance_intervals_years].values.tolist())
    )

    if len(maintain_intervals) > 1:
        df_maintain = []
        for interval in routine_maintain_intervals:
            if interval > 0:
                maintain_years = np.arange(start_year, end_year + 1, interval)
                df_mod = df[df[maintenance_intervals_years] == interval]
                df_mod[maintain_years[1:]] = df_mod[maintain_years[1:]].add(
                    df_mod[maintenance_cost], axis=0
                )
                df_maintain.append(df_mod)
        if len(df_maintain) > 0:
            df_maintain = pd.concat(
                df_maintain, axis=0, ignore_index=True
            ).fillna(0)
        else:
            df_maintain = df.copy()
    else:
        if maintain_intervals[0] > 0:
            df_maintain = df.copy()
            maintain_years = np.arange(
                start_year, end_year + 1, maintain_intervals[0]
            )
            df_maintain[maintain_years[1:]] = df_maintain[
                maintain_years[1:]
            ].add(df_maintain[maintenance_cost], axis=0)
        else:
            df_maintain = df.copy()

    return df_maintain


def assign_costs_over_time(
    df, asset_id, start_year=2019, end_year=2100, discounting_rate=10
):
    timeseries = np.arange(start_year, end_year + 1, 1)
    df[timeseries] = 0
    df[timeseries[0]] += df["initial_investment_cost"]

    df = assign_maintenance_cost_over_time(
        df,
        "routine_maintenance_intervals_years",
        "routine_maintenance_cost",
        start_year=start_year,
        end_year=end_year,
    )
    df = assign_maintenance_cost_over_time(
        df,
        "periodic_maintenance_intervals_years",
        "periodic_maintenance_cost",
        start_year=start_year,
        end_year=end_year,
    )
    return df[
        [
            asset_id,
            "adaptation_option",
            "asset_adaptation_cost",
        ] + list(timeseries)
    ]


def get_dimension_factor(x):
    dimension_type = x["asset_dimensions"]
    if dimension_type in ["length", "perimeter"]:
        dimension = x.geometry.length
    elif dimension_type == "area":
        dimension = x.geometry.area
    else:
        dimension = 1

    change_type = x["change_parameter"]
    if change_type == "flood depth":
        cost_unit = "J$/m"
    else:
        cost_unit = "J$"

    return dimension, cost_unit


def get_coastal_dimension_factor(x, flood_params, protect_feature, protect_dict, asset_id):
    rcp, rp, epoch = flood_params
    flood_id_col = f"flood_id_rcp_{int(rcp*10)}{epoch}_rp_{rp}"
    
    # Look up flood_id for this asset
    asset_match = protect_dict[protect_dict[asset_id] == x[asset_id]]
    if asset_match.empty:
        dimension = 1
    else:
        value = asset_match[flood_id_col].iloc[0]
        if pd.isna(value):  # Handle NaN values explicitly
            dimension = 1
        else:
            flood_id = int(value)
            # Look up coastline length for this flood_id
            feature_match = protect_feature[protect_feature["polygon_id"] == flood_id]
            if feature_match.empty:
                dimension = 1
            else:
                dimension = feature_match["coastline_length"].iloc[0]
            # print (f"{x[asset_id]}----{flood_id}------{dimension}")
    
    # Determine cost unit
    change_type = x["change_parameter"]
    if change_type == "flood depth":
        cost_unit = "J$/m"
    else:
        cost_unit = "J$"
    
    return dimension, cost_unit


def get_adaptation_options_costs(asset_df, asset_id, hazard_label, flood_params, protection_feature_breakdown, protection_asset_dict):
    if hazard_label != 'coastal':
        asset_df["dimension_cost_factor"] = asset_df.progress_apply(
            lambda x: get_dimension_factor(x), axis=1
        )
    else:
        protect_feature = pd.read_csv(protection_feature_breakdown)
        protect_dict = pd.read_parquet(protection_asset_dict)
        asset_df["dimension_cost_factor"] = asset_df.progress_apply(
            lambda x: get_coastal_dimension_factor(x, flood_params, protect_feature, protect_dict, asset_id), axis=1
        )
    asset_df[["dimension_factor", "asset_adaptation_cost"]] = asset_df[
        "dimension_cost_factor"
    ].apply(pd.Series)
    asset_df["cost_multiplier"] = (
        asset_df["dimension_factor"]
        * asset_df["currency_conversion"]
        * asset_df["asset_dimension_coversion"]
    )
    asset_df[
        [
            "initial_investment_cost",
            "periodic_maintenance_cost",
            "routine_maintenance_cost",
        ]
    ] = asset_df[
        [
            "initial_investment_cost_per_unit",
            "periodic_maintenance_cost_per_unit",
            "routine_maintenance_cost_per_unit",
        ]
    ].multiply(
        asset_df["cost_multiplier"], axis="index"
    )
    
    return asset_df[
        [
            asset_id,
            "adaptation_option",
            "option_unit_cost",
            "initial_investment_cost_per_unit",
            "periodic_maintenance_cost_per_unit",
            "routine_maintenance_cost_per_unit",
            "periodic_maintenance_intervals_years",
            "routine_maintenance_intervals_years",
            "asset_adaptation_cost",
            "initial_investment_cost",
            "periodic_maintenance_cost",
            "routine_maintenance_cost",
        ]
    ]


def get_adaptation_options_costs_roads(asset_df, adapt_costs, asset_id, hazard_label, flood_params, protection_feature_breakdown, protection_asset_dict):
    road_costs = adapt_costs[adapt_costs["asset_description"] == "roads"]
    roads_df = []
    for rc in road_costs.itertuples():
        if rc.asset_details == "Roads-2L":
            df = asset_df[asset_df["lanes"] == 2]
            df["lane_factor"] = 1
        elif rc.asset_details == "Roads-4L":
            df = asset_df[asset_df["lanes"] != 2]
            df["lane_factor"] = df["lanes"] / 4.0
        else:
            df = asset_df.copy()
            df["lane_factor"] = 1

        for column in [
            "asset_dimensions",
            "change_parameter",
            "currency_conversion",
            "asset_dimension_coversion",
            "adaptation_option",
            "option_unit_cost",
            "initial_investment_cost_per_unit",
            "periodic_maintenance_cost_per_unit",
            "routine_maintenance_cost_per_unit",
            "periodic_maintenance_intervals_years",
            "routine_maintenance_intervals_years",
        ]:
            if column == "asset_dimension_coversion":
                df[column] = getattr(rc, column) * df["lane_factor"]
            else:
                df[column] = getattr(rc, column)
        df = get_adaptation_options_costs(df, asset_id, hazard_label, flood_params, protection_feature_breakdown, protection_asset_dict)
        roads_df.append(df)

    roads_df = pd.concat(roads_df, axis=0, ignore_index=True)
    return roads_df


def write_empty_output(asset_unit_costs_csv: str, asset_timeseries_csv: str, baseline_year: int, projection_end_year: int) -> None:
    pd.DataFrame([]).to_csv(asset_unit_costs_csv, index=False)
    # Schema necessary for subsequent script to open timeseries file,
    # so write the header but nothing else
    pd.DataFrame(
        [],
        columns=[
            "undefined_asset_id",
            "adaptation_option",
            "asset_adaptation_cost",
        ] + list(map(str, range(baseline_year, projection_end_year + 1)))
    ).to_csv(asset_timeseries_csv, index=False)
    return


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        readable=True
    ),
    help="Path to the asset layer data",
)
@click.option(
    "--asset-file",
    "-b",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        readable=True
    ),
    help="Path to asset gpkg file",
)
@click.option(
    "--cost-file",
    "-c",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        readable=True
    ),
    help="Path to adaptation cost data",
)
@click.option(
    "--protection-asset-dict",
    "-pa",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        readable=True
    ),
    help="Path to directory with network assets to flood protection area map",
)
@click.option(
    "--protection-feature-breakdown",
    "-pfb",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        readable=True
    ),
    help="Path to directory with breakdown of network assets for each flood protection feature",
)
@click.option(
    "--hazard-label",
    "-h",
    required=True,
    help="hazard label",
)
@click.option(
    "--asset-gpkg",
    "-g",
    required=True,
    help="asset_gpkg value in the network CSV",
)
@click.option(
    "--asset-layer",
    "-l",
    required=True,
    help="asset_layer value in the network CSV",
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        readable=True
    ),
    help="Path to the output directory",
)
@click.option(
    "--baseline-year",
    "-y",
    default=2019,
    required=False,
    type=int,
    help="Baseline year",
)
@click.option(
    "--projection-end-year",
    "-p",
    default=2100,
    required=False,
    type=int,
    help="Projection end year",
)
@click.option(
    "--rcp",
    "-rcp",
    required=True,
    type=float,
    help="RPS value",
)
@click.option(
    "--rp",
    "-rp",
    required=True,
    type=int,
    help="RP value",
)
@click.option(
    "--discounting-rate",
    "-d",
    default=10,
    required=False,
    type=float,
    help="Discounting rate",
)
@click.option(
    "--epsg",
    "-e",
    default=3448,
    required=False,
    type=int,
    help="EPSG",
)
def adaptation_options_costs(
    network_csv,
    asset_file,
    cost_file,
    protection_asset_dict,
    protection_feature_breakdown,
    hazard_label,
    asset_gpkg,
    asset_layer,
    output_dir,
    baseline_year,
    projection_end_year,
    rcp,
    rp,
    discounting_rate,
    epsg,
):

    # make output folders
    adaptation_results = os.path.join(output_dir, "adaptation_costs")
    os.makedirs(adaptation_results, exist_ok=True)
    hazard_outputs = os.path.join(adaptation_results, f"{hazard_label}_costs")
    os.makedirs(hazard_outputs, exist_ok=True)
    # output filepaths
    asset_unit_costs_csv = os.path.join(hazard_outputs, f"{asset_gpkg}_{asset_layer}_adaptation_unit_costs.csv")
    asset_timeseries_csv = os.path.join(hazard_outputs, f"{asset_gpkg}_{asset_layer}_adaptation_timeseries_and_npvs.csv")

    logging.info("Read network metadata")
    network_metadata = pd.read_csv(network_csv)
    network_csv_null_value = "none"
    network_layer = network_metadata[
        (network_metadata["asset_gpkg"] == asset_gpkg)
        & (network_metadata["asset_layer"] == asset_layer)
    ].squeeze()

    logging.info("Read adaptation option costs")
    all_adaptation_costs = pd.read_excel(cost_file, sheet_name="Sheet1").fillna(0)
    hazard_adaptation_costs = all_adaptation_costs[all_adaptation_costs["hazard"] == hazard_label]

    asset_id_col = network_layer.asset_id_column
    asset_type_lookup = f"{hazard_label}_asset_damage_lookup_column"
    asset_type_col = network_layer[asset_type_lookup]

    logging.info("Read network layer (per-asset attributes)")
    assets = gpd.read_file(asset_file, layer=asset_layer).to_crs(epsg=epsg)

    logging.info("Lookup costs for network assets")
    if network_layer.asset_description != "roads":

        if asset_type_col == network_csv_null_value:
            logging.info(f"{network_csv}::{asset_type_lookup} for {asset_gpkg}::{asset_layer} is {network_csv_null_value}, skipping...")
            write_empty_output(asset_unit_costs_csv, asset_timeseries_csv, baseline_year, projection_end_year)
            return

        assets = assets[assets[asset_type_col].isin(hazard_adaptation_costs.asset_name)]
        if assets.empty:
            logging.info("No adaptation options listed for assets, skipping...")
            write_empty_output(asset_unit_costs_csv, asset_timeseries_csv, baseline_year, projection_end_year)
            return

        # Join adaptation option costs 'asset_name' on assets `asset_type_col`
        assets = pd.merge(
            assets,
            hazard_adaptation_costs,
            how="left",
            left_on=asset_type_col,
            right_on="asset_name",
        )
        assets = get_adaptation_options_costs(
            assets,
            asset_id_col,
            hazard_label,
            [rcp, rp, projection_end_year],
            protection_feature_breakdown,
            protection_asset_dict
        )
    else:
        assets = get_adaptation_options_costs_roads(
            assets,
            hazard_adaptation_costs,
            asset_id_col,
            hazard_label,
            [rcp, rp, projection_end_year],
            protection_feature_breakdown,
            protection_asset_dict
        )

    logging.info("Write out per asset costs")
    protect_dict = pd.read_parquet(protection_asset_dict)
    assets = assets[assets[asset_id_col].isin(protect_dict[asset_id_col])]
    logging.info(f"\n{assets}")
    assets.to_csv(asset_unit_costs_csv, index=False)

    logging.info("Calculate costs over time")
    assets = assign_costs_over_time(
        assets,
        asset_id_col,
        start_year=baseline_year,
        end_year=projection_end_year,
        discounting_rate=discounting_rate,
    )
    discount_rate = calculate_discounting_rate_factor(
        discount_rate=discounting_rate,
        start_year=baseline_year,
        end_year=projection_end_year,
        maintain_period=1,
    )
    discounted_assets = assets.copy()
    years: np.ndarray = np.arange(baseline_year, projection_end_year + 1, 1)
    discounted_assets[years] = np.multiply(discounted_assets[years], discount_rate)
    assets["adapt_cost_npv"] = discounted_assets[years].sum(axis=1)

    logging.info("Write out per asset per year costs")
    assets = assets[assets[asset_id_col].isin(protect_dict[asset_id_col])]
    logging.info(f"\n{assets}")
    assets.to_csv(asset_timeseries_csv, index=False)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )
    adaptation_options_costs()
