"""Estimate adaptation options costs and benefits

"""

import click
import os
import logging
import warnings

import pandas as pd

import geopandas as gpd
import numpy as np
from tqdm import tqdm

from jamaica_infrastructure.analysis.utils import (
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


def get_adaptation_options_costs(asset_df, asset_id):
    asset_df["dimension_cost_factor"] = asset_df.progress_apply(
        lambda x: get_dimension_factor(x), axis=1
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


def get_adaptation_options_costs_roads(asset_df, adapt_costs, asset_id):
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
        df = get_adaptation_options_costs(df, asset_id)
        roads_df.append(df)

    roads_df = pd.concat(roads_df, axis=0, ignore_index=True)
    return roads_df


@click.command()
@click.version_option("1.0")
@click.option(
    "--asset-data",
    "-a",
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
def adaptation_options_costs(
    asset_data,
    asset_file,
    cost_file,
    hazard_label,
    asset_gpkg,
    asset_layer,
    output_dir,
):
    epsg_jamaica = 3448
    baseline_year = 2019
    projection_end_year = 2100
    discounting_rate = 10

    cost_df = pd.read_excel(
        cost_file,
        sheet_name="Sheet1",
    ).fillna(0)

    asset_df = pd.read_csv(asset_data)
    asset_data_details = asset_df[
        (asset_df["asset_gpkg"] == asset_gpkg)
        & (asset_df["asset_layer"] == asset_layer)
    ]
    asset_df = gpd.read_file(
        asset_file,
        layer=asset_layer,
    )

    dsc_rate = calculate_discounting_rate_factor(
        discount_rate=discounting_rate,
        start_year=baseline_year,
        end_year=projection_end_year,
        maintain_period=1,
    )
    cost_timeseries = np.arange(baseline_year, projection_end_year + 1, 1)

    adaptation_results = os.path.join(output_dir, "adaptation_costs")
    if os.path.exists(adaptation_results) is False:
        os.mkdir(adaptation_results)

    hazard_outputs = os.path.join(adaptation_results, f"{hazard_label}_costs")
    if os.path.exists(hazard_outputs) is False:
        os.mkdir(hazard_outputs)

    adapt_costs = cost_df[cost_df["hazard"] == hazard_label]
    cost_description = list(
        set(adapt_costs["asset_description"].values.tolist())
    )
    costed_assets = list(set(adapt_costs["asset_name"].values.tolist()))

    adapt_assets = asset_data_details[
        asset_data_details["asset_description"].isin(cost_description)
    ]

    asset_info = adapt_assets.squeeze()
    asset_id = asset_info.asset_id_column
    asset_hazard = getattr(
        asset_info, f"{hazard_label}_asset_damage_lookup_column"
    )
    asset_df = asset_df.to_crs(epsg=epsg_jamaica)
    if asset_info.asset_description != "roads":
        asset_df = asset_df[asset_df[asset_hazard].isin(costed_assets)]
        asset_df = pd.merge(
            asset_df,
            adapt_costs,
            how="left",
            left_on=asset_hazard,
            right_on="asset_name",
        )
        asset_df = get_adaptation_options_costs(asset_df, asset_id)
    else:
        asset_df = get_adaptation_options_costs_roads(
            asset_df, adapt_costs, asset_id
        )

    asset_unit_costs_csv = os.path.join(
        hazard_outputs,
        f"{asset_gpkg}_{asset_layer}_adaptation_unit_costs.csv",
    )
    asset_df.to_csv(
        asset_unit_costs_csv,
        index=False,
    )
    logging.info(asset_unit_costs_csv)
    asset_df = assign_costs_over_time(
        asset_df,
        asset_id,
        start_year=baseline_year,
        end_year=projection_end_year,
        discounting_rate=discounting_rate,
    )

    df = asset_df.copy()
    df[cost_timeseries] = np.multiply(df[cost_timeseries], dsc_rate)
    asset_df["adapt_cost_npv"] = df[cost_timeseries].sum(axis=1)
    asset_timeseries_csv = os.path.join(
        hazard_outputs,
        f"{asset_gpkg}_{asset_layer}_adaptation_timeseries_and_npvs.csv",
    )
    asset_df.to_csv(
        asset_timeseries_csv,
        index=False,
    )
    logging.info(asset_timeseries_csv)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )
    adaptation_options_costs()
