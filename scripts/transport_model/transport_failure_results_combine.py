"""Do a transport failure analysis with rerouting
"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
import igraph as ig
import ast
from utils import *
from tqdm import tqdm

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
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    hourly_wage = (
        0.4 * (1 + 0.454) * 235.25
    )  # Between 200 - 500 JMD for 2012 stats, 45.4% inflation in currency
    results_dir = os.path.join(results_path, "transport_failures")
    all_failures = []
    bridge_failures = []
    for root, dirs, files in os.walk(os.path.join(results_dir, "scenario_results")):
        for file in files:
            if file.startswith("single_bridge_failures_scenarios"):
                df = pd.read_csv(os.path.join(root, file))
                bridge_failures.append(df)
            elif file.endswith(".csv"):
                df = pd.read_csv(os.path.join(root, file))
                df.rename(columns={"no_access_x": "no_access"}, inplace=True)
                all_failures.append(
                    df[
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
    bridge_failures = pd.concat(bridge_failures, axis=0, ignore_index=True).fillna(0)
    labour_flows = gpd.read_file(
        os.path.join(results_path, "flow_mapping", "labour_trips_and_activity.gpkg"),
        layer="edges",
    )
    bridges = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="nodes",
    )
    bridges = bridges[bridges["asset_type"] == "bridge"]["node_id"].values.tolist()
    edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    bridge_roads = []
    for bridge in bridges:
        bridge_edges = edges[
            (edges["from_node"] == bridge) | (edges["to_node"] == bridge)
        ]
        if len(bridge_edges.index) > 0:
            bridge_edges["node_id"] = bridge
            bridge_roads.append(bridge_edges[["node_id", "edge_id"]])
    bridge_roads = pd.concat(bridge_roads, axis=0, ignore_index=True)
    bridge_labour_trips = pd.read_csv(
        os.path.join(
            results_path,
            "flow_mapping",
            "origins_destinations_labour_economic_activity.csv",
        )
    )
    bridge_failures = pd.merge(
        bridge_failures, bridge_labour_trips, how="left", on=["node_id"]
    ).fillna(0)
    # bridge_trade_trips = pd.read_csv(os.path.join(results_path,
    #                                 'flow_mapping',
    #                                 'origins_destinations_trade_economic_activity.csv'))
    # bridge_failures = pd.merge(bridge_failures,
    #                         bridge_trade_trips,
    #                         how="left",on=["node_id"])
    all_failures = pd.merge(
        all_failures,
        labour_flows[["edge_id", "working_trips", "GDP_to_trips"]],
        how="left",
        on=["edge_id"],
    ).fillna(0)
    bridge_roads = pd.merge(bridge_roads, all_failures, how="left", on=["edge_id"])
    bridge_roads["node_degree"] = bridge_roads.groupby("node_id")["node_id"].transform(
        "count"
    )
    bridge_roads_access = (
        bridge_roads.groupby(["node_id", "node_degree"])["no_access"]
        .sum()
        .reset_index()
    )
    bridge_roads_access["no_access"] = bridge_roads_access.apply(
        lambda x: 1 if x.no_access > 0 else 0, axis=1
    )
    bridge_roads_trips = (
        bridge_roads.groupby(["node_id", "node_degree"])["working_trips"]
        .sum()
        .reset_index()
    )
    bridge_roads_gdp = (
        bridge_roads.groupby(["node_id", "node_degree"])["GDP_to_trips"]
        .max()
        .reset_index()
    )
    bridge_roads_access = pd.merge(
        bridge_roads_access,
        bridge_roads_trips,
        how="left",
        on=["node_id", "node_degree"],
    )
    bridge_roads_access = pd.merge(
        bridge_roads_access, bridge_roads_gdp, how="left", on=["node_id", "node_degree"]
    )
    del bridge_roads, bridge_roads_trips, bridge_roads_gdp
    bridge_roads_access["working_trips_thru"] = (
        1 - bridge_roads_access["no_access"]
    ) * bridge_roads_access["working_trips"]
    bridge_roads_access["labour_GDP_thru"] = (
        bridge_roads_access["no_access"] * bridge_roads_access["GDP_to_trips"]
    )

    bridge_failures = pd.merge(
        bridge_failures,
        bridge_roads_access[
            ["node_id", "node_degree", "working_trips_thru", "labour_GDP_thru"]
        ],
        how="left",
        on=["node_id"],
    ).fillna(0)

    bridge_failures["working_trips"] = (
        bridge_failures["working_trips_thru"]
        + bridge_failures["d_trip"]
        - bridge_failures["o_trip"]
    ) / bridge_failures["node_degree"]
    bridge_failures["GDP_to_trips"] = bridge_failures.apply(
        lambda x: max(x["GDP_to_trips"], x["labour_GDP_thru"]), axis=1
    )

    bridge_failures = get_failure_estimates(bridge_failures, "node_id", hourly_wage)
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
    bridge_failures.to_csv(
        os.path.join(
            results_dir, "single_point_failure_road_bridges_economic_losses.csv"
        ),
        index=False,
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
