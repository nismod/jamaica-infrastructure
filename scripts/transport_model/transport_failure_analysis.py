"""Transport failure analysis with rerouting."""

import json
import logging
import os
import sys

import pandas as pd
import geopandas as gpd
import igraph as ig
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import (
    load_config,
    map_nearest_locations_and_create_lines,
    network_od_paths_assembly,
)
import jamaica_infrastructure.transport.flow as tf

tqdm.pandas()
epsg_jamaica = 3448


def route_areas_to_nearest_ports(
    areas,
    areas_id,
    areas_gdp,
    nodes,
    edges,
    ports,
    port_weight,
    connection_type="areas",
    trade_type="import",
    include_rail=False,
):
    network_columns = [
        "from_node",
        "to_node",
        "edge_id",
        "from_mode",
        "to_mode",
        "length_m",
        "speed",
        "time",
        "geometry",
    ]
    nearest_roads = map_nearest_locations_and_create_lines(
        areas.copy(),
        nodes[nodes["mode"] == "road"].copy(),
        areas_id,
        "node_id",
        connection_type,
        "road",
    )
    nearest_roads["edge_id"] = nearest_roads.progress_apply(
        lambda x: f"{connection_type}roade_{x.name}", axis=1
    )
    nearest_roads["speed"] = 10.0
    nearest_roads["time"] = 0.001 * nearest_roads["length_m"] / nearest_roads["speed"]

    if include_rail is True:
        nearest_stations = map_nearest_locations_and_create_lines(
            areas.copy(),
            nodes[nodes["mode"] == "rail"].copy(),
            areas_id,
            "node_id",
            connection_type,
            "rail",
        )
        nearest_stations["edge_id"] = nearest_stations.progress_apply(
            lambda x: f"{connection_type}raile_{x.name}", axis=1
        )
        nearest_stations["speed"] = 10.0
        nearest_stations["time"] = (
            0.001 * nearest_stations["length_m"] / nearest_stations["speed"]
        )

        multi_edges = edges[edges["from_mode"] != edges["to_mode"]]
        multi_edges = multi_edges[multi_edges["length_m"] >= 5000][
            "edge_id"
        ].values.tolist()
        edges = edges[~edges["edge_id"].isin(multi_edges)]
        network = pd.concat(
            [
                edges[network_columns],
                nearest_stations[network_columns],
                nearest_roads[network_columns],
            ],
            axis=0,
            ignore_index=True,
        )[network_columns]
        area_edges = pd.concat(
            [nearest_roads[network_columns], nearest_stations[network_columns]],
            axis=0,
            ignore_index=True,
        )[network_columns]
    else:
        edges = edges[(edges["from_mode"] != "rail") & (edges["to_mode"] != "rail")]
        network = pd.concat(
            [edges[network_columns], nearest_roads[network_columns]],
            axis=0,
            ignore_index=True,
        )[network_columns]
        area_edges = nearest_roads[network_columns]
        # edges = edges[(edges["from_mode"] != "rail") & (edges["to_mode"] != "rail")]

    G = ig.Graph.TupleList(
        network.itertuples(index=False), edge_attrs=list(network.columns)[2:]
    )

    all_ports = ports["node_id"].values.tolist()

    # od_pairs = [list(zip([b]*len(all_ports),all_ports)) for b in areas[areas_id].values.tolist()]
    od_pairs = [
        list(zip(all_ports, [b] * len(all_ports)))
        for b in areas[areas_id].values.tolist()
    ]
    od_pairs = [item for sublist in od_pairs for item in sublist]
    od_pairs = pd.DataFrame(od_pairs, columns=["origin_id", "destination_id"])
    od_pairs = pd.merge(
        od_pairs,
        areas[[areas_id, areas_gdp]],
        how="left",
        left_on=["destination_id"],
        right_on=[areas_id],
    )
    flow_paths = network_od_paths_assembly(
        od_pairs[["origin_id", "destination_id", areas_gdp]], G, "time", areas_gdp
    )
    flow_paths = flow_paths.sort_values(by="gcost")
    flow_paths = flow_paths.drop_duplicates(subset=["destination_id"], keep="first")
    flow_paths = flow_paths[
        ["origin_id", "destination_id", "edge_path", "gcost", areas_gdp]
    ]

    flow_paths_areas = flow_paths.groupby(["origin_id"])[areas_gdp].sum().reset_index()
    flow_paths_areas.rename(columns={areas_gdp: "tot_GDP"}, inplace=True)
    flow_paths = pd.merge(flow_paths, flow_paths_areas, how="left", on=["origin_id"])
    del flow_paths_areas
    flow_paths = pd.merge(
        flow_paths,
        ports[["node_id", port_weight]],
        how="left",
        left_on=["origin_id"],
        right_on=["node_id"],
    )
    flow_paths["trade_wt"] = (
        flow_paths[port_weight] * flow_paths[areas_gdp] / flow_paths["tot_GDP"]
    )
    flow_paths = flow_paths[
        ["origin_id", "destination_id", "edge_path", "gcost", "trade_wt"]
    ]

    if trade_type == "export":
        flow_paths.columns = [
            "destination_id",
            "origin_id",
            "edge_path",
            "gcost",
            "trade_wt",
        ]

    return flow_paths, area_edges


def port_import_exports(ports, tons_column, trade_type):
    ports = ports.drop_duplicates(subset=["name"], keep="first")
    ports = ports[["node_id", "name", tons_column, "geometry"]]
    ports[f"{trade_type}_wt"] = ports[tons_column] / ports[tons_column].sum()
    ports["geometry"] = ports.progress_apply(lambda x: x.geometry.centroid, axis=1)
    ports = ports.to_crs(epsg=epsg_jamaica)

    return ports


def filter_sector_from_buildings(buildings_dataframe, sector_code, subsector_code):
    get_sector = buildings_dataframe[buildings_dataframe[f"{sector_code}_GDP"] > 0]
    get_sector["find_subsector"] = get_sector.progress_apply(
        lambda x: 1 if subsector_code in str(x.subsector_code) else 0, axis=1
    )
    return get_sector[get_sector["find_subsector"] == 1]


def read_flow_data(data_path: str):
    """
    Read combined flow data from disk ready for transport failure disruption.
    """

    network_data: dict = {}
    for flow_type in ("labour", "trade"):
        output_flow_dir = os.path.join(data_path, flow_type)
        logging.info(f"Reading {flow_type} network")
        network_df = gpd.read_parquet(os.path.join(output_flow_dir, "network.gpq"))
        network: ig.Graph = ig.Graph.TupleList(
            network_df.itertuples(index=False),
            edge_attrs=['edge_id', 'from_mode', 'to_mode', 'length_m', 'speed', 'time', 'geometry']
        )
        logging.info(f"Reading {flow_type} flows")
        flows = pd.read_parquet(os.path.join(output_flow_dir, "flows.pq"))
        logging.info(f"Reading {flow_type} edge indices")
        edge_indexes = pd.read_parquet(os.path.join(output_flow_dir, "edge_indexes.pq"))
        network_data[flow_type] = {
            "network": network,
            "flows": flows,
            "edge_indexes": edge_indexes.to_dict()["edge_indexes"],
        }

    logging.info("Reading combined flows")
    all_flows = pd.read_parquet(os.path.join(data_path, "all_flows.pq"))

    logging.info("Reading trade sectors")
    with open(os.path.join(data_path, "trade", "trade_sectors.json"), "r") as fp:
        trade_sectors = json.load(fp)

    return network_data, all_flows, trade_sectors


def main(config, min_edge_number, max_edge_number):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    logging.info(
        f"Running transport failure analysis for edge positions {min_edge_number} -> {max_edge_number}"
    )

    logging.info("Reading network data")

    network_dictionary, all_flows, trade_sectors = read_flow_data(
        os.path.join(results_path, "transport_failures", "nominal")
    )

    edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    rail_edges = edges[(edges["from_mode"] == "rail") & (edges["to_mode"] == "rail")][
        "edge_id"
    ].values.tolist()
    road_edges = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")][
        "edge_id"
    ].values.tolist()
    edge_fail = rail_edges + road_edges

    logging.info("Failing edges and reallocating flows")

    if max_edge_number > len(edge_fail):
        max_edge_number = len(edge_fail)

    edge_fail_results = []
    # TODO: extremely slow, can we parallelise on this loop? are we already?
    for edge_number in range(min_edge_number, max_edge_number):
        edge = edge_fail[edge_number]
        logging.info(f"Failing {edge}")
        for networks in network_dictionary.values():
            edge_fail_results += tf.igraph_scenario_edge_failures_premade_network(
                # we will remove edges from the graph, only operate on a copy
                networks["network"].copy(),
                [edge],
                networks["flows"],
                networks["edge_indexes"],
                "edge_path",
                "time",
            )

    logging.info("Done failing edges")

    logging.info("Calculating resulting economic losses")

    # 0.6 - 2.1% of the value per day
    # so we need to make an assumption on the average wage per working person,
    # say 200 USD per day. Then if a road is disrupted which has 1000 daily trips and
    # they have to be rerouted with an hour, the cost would be: 0.4 * 200 * 1/24 * 100 = 333 USD
    # So corrected for inflation in 2019 values, this would be 1.2-2.9 USD per hour of value of time for business related trips
    hourly_wage = (
        0.4 * (1 + 0.454) * 235.25
    )  # Between 200 - 500 JMD for 2012 stats, 45.4% inflation in currency
    trade_effect = 0.02  # 2% of the value of trade will be affected by rerouting

    edge_fail_results = pd.DataFrame(edge_fail_results)
    edge_fail_results = pd.merge(
        edge_fail_results, all_flows, how="left", on=["origin_id", "destination_id"]
    ).fillna(0)
    edge_fail_results["total_trade"] = edge_fail_results[
        [f"{t}_trade" for t in trade_sectors]
    ].sum(axis=1)
    # print (edge_fail_results)

    edge_fail_results["time_loss"] = (1 - edge_fail_results["no_access"]) * (
        edge_fail_results["new_cost"] - edge_fail_results["gcost"]
    )
    edge_fail_results["labour_rerouting_loss"] = (
        hourly_wage
        * edge_fail_results["time_loss"]
        * edge_fail_results["working_trips"]
    )
    edge_fail_results["trade_rerouting_loss"] = (
        trade_effect * edge_fail_results["time_loss"] * edge_fail_results["total_trade"]
    )
    edge_fail_results["labour_gdp_loss"] = (
        edge_fail_results["no_access"] * edge_fail_results["GDP_to_trips"]
    )
    edge_fail_results["trade_loss"] = (
        edge_fail_results["no_access"] * edge_fail_results["total_trade"]
    )
    # print (edge_fail_results[["edge_id","no_access","time_loss",
    #                         "labour_rerouting_loss","trade_rerouting_loss",
    #                         "labour_gdp_loss","trade_loss"]])

    index_cols = ["edge_id", "no_access"]
    losses = (
        edge_fail_results.loc[
            :,
            index_cols
            + [
                "time_loss",
                "labour_rerouting_loss",
                "trade_rerouting_loss",
                "labour_gdp_loss",
                "trade_loss",
            ],
        ]
        .groupby(index_cols)
        .sum()
        .reset_index()
    )
    rerouting_times_min = (
        edge_fail_results.loc[:, index_cols + ["time_loss"]]
        .groupby(index_cols)
        .min()
        .reset_index()
    )
    rerouting_times_min = rerouting_times_min.rename(
        columns={"time_loss": "min_trip_time_loss"}
    )
    rerouting_times_max = (
        edge_fail_results.loc[:, index_cols + ["time_loss"]]
        .groupby(index_cols)
        .max()
        .reset_index()
    )
    rerouting_times_max = rerouting_times_max.rename(
        columns={"time_loss": "max_trip_time_loss"}
    )
    rerouting_times_mean = (
        edge_fail_results.loc[:, index_cols + ["time_loss"]]
        .groupby(index_cols)
        .mean()
        .reset_index()
    )
    rerouting_times_mean = rerouting_times_mean.rename(
        columns={"time_loss": "mean_trip_time_loss"}
    )

    logging.info("Writing results to disk")

    losses = pd.merge(losses, rerouting_times_min.drop(columns=["no_access"]), how="left", on=["edge_id"])
    losses = pd.merge(losses, rerouting_times_max.drop(columns=["no_access"]), how="left", on=["edge_id"])
    losses = pd.merge(losses, rerouting_times_mean.drop(columns=["no_access"]), how="left", on=["edge_id"])

    output_dir = os.path.join(results_path, "transport_failures", "scenario_results")
    losses.to_csv(
        os.path.join(output_dir, f"single_link_failures_scenarios_{min_edge_number}_{max_edge_number}.csv"),
        index=False,
    )


if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )

    CONFIG = load_config()
    try:
        min_edge_number = int(sys.argv[2])
        max_edge_number = int(sys.argv[3])
    except IndexError:
        logging.info("Got arguments", sys.argv)
        exit()

    main(CONFIG, min_edge_number, max_edge_number)
