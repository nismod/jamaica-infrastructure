"""Combine results of bridge and other edge failures to estimate economic losses."""

import glob
import logging
import os

import pandas as pd
import geopandas as gpd
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import load_config

tqdm.pandas()
epsg_jamaica = 3448


def get_failure_estimates(failure_df, id_column, hourly_wage):
    failure_df["mean_labour_rerouting_loss"] = (
        hourly_wage * failure_df["mean_trip_time_loss"] * failure_df["working_trips"]
    )
    failure_df["labour_gdp_loss"] = failure_df["no_access"] * failure_df["GDP_to_trips"]
    failure_df["rerouting_loss"] = (1 - failure_df["no_access"]) * (
        failure_df["mean_labour_rerouting_loss"] + failure_df["trade_rerouting_loss"]
    )
    failure_df["isolation_loss"] = failure_df["no_access"] * (
        failure_df["labour_gdp_loss"] + failure_df["trade_loss"]
    )
    failure_df["economic_loss"] = (
        failure_df["rerouting_loss"] + failure_df["isolation_loss"]
    )

    failure_df = (
        failure_df.groupby([id_column])[
            ["rerouting_loss", "isolation_loss", "economic_loss"]
        ]
        .sum()
        .reset_index()
    )
    failure_df["loss_unit"] = "JD/day"

    return failure_df


def main(config):
    results_path = config["paths"]["output"]

    hourly_wage = (
        0.4 * (1 + 0.454) * 235.25
    )  # Between 200 - 500 JMD for 2012 stats, 45.4% inflation in currency
    results_dir = os.path.join(results_path, "transport_failures")

    logging.info("Reading road failure results")
    all_failures = []
    for file_path in glob.glob(os.path.join(results_dir, "scenario_results", "single_link_failures_scenarios*.csv")):
        df = pd.read_csv(file_path)
        all_failures.append(
            df.loc[
                :,
                [
                    "edge_id",
                    "no_access",
                    "time_loss",
                    "labour_rerouting_loss",
                    "trade_rerouting_loss",
                    "labour_gdp_loss",
                    "trade_loss",
                    "min_trip_time_loss",
                    "max_trip_time_loss",
                    "mean_trip_time_loss",
                ]
            ]
        )
    all_failures = pd.concat(all_failures, axis=0, ignore_index=True).fillna(0)

    logging.info("Reading labour flows")
    labour_flows = pd.read_parquet(os.path.join(results_path, "flow_mapping", "labour_trips_and_activity.pq"))
    breakpoint()

    all_failures = pd.merge(
        all_failures,
        labour_flows[["edge_id", "working_trips", "GDP_to_trips"]],
        how="left",
        on=["edge_id"],
    ).fillna(0)

    all_failures = get_failure_estimates(all_failures, "edge_id", hourly_wage)

    results_dir = os.path.join(
        results_path, "economic_losses", "single_failure_scenarios"
    )
    all_failures.to_csv(
        os.path.join(
            results_dir, "single_point_failure_road_rail_edges_economic_losses.csv"
        ),
        index=False,
    )


if __name__ == "__main__":
    CONFIG = load_config()
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    main(CONFIG)
