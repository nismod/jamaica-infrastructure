"""Map labour to work trips as a proxy to map GDP onto roads
"""

import os

import pandas as pd
import geopandas as gpd
import numpy as np
import igraph as ig
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import (
    ckdnearest,
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


def main(config):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    """Get the ports and multimodal network
    """
    nodes = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="nodes",
    )
    nodes = nodes[nodes["mode"] == "road"]
    nodes = nodes.to_crs(epsg=epsg_jamaica)
    columns = [
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
    edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    network = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")][
        columns
    ]
    # 0.6 - 2.1% of the value per day
    # so we need to make an assumption on the average wage per working person,
    # say 200 USD per day. Then if a road is disrupted which has 1000 daily trips and
    # they have to be rerouted with an hour, the cost would be: 0.4 * 200 * 1/24 * 100 = 333 USD
    # So corrected for inflation in 2019 values, this would be 1.2-2.9 USD per hour of value of time for business related trips

    G = ig.Graph.TupleList(
        network.itertuples(index=False), edge_attrs=list(network.columns)[2:]
    )
    """Get the buildings with populations and commercial activites assigned to them
    """
    buildings = gpd.read_file(
        os.path.join(
            processed_data_path,
            "buildings",
            "buildings_assigned_economic_activity.gpkg",
        ),
        layer="areas",
    )
    buildings["osm_id"] = buildings.progress_apply(
        lambda x: f"building_{x.osm_id}", axis=1
    )
    buildings["geometry"] = buildings.progress_apply(
        lambda x: x.geometry.centroid, axis=1
    )
    buildings = buildings.to_crs(epsg=epsg_jamaica)

    print(buildings)

    population_year = 2019
    population = gpd.read_file(
        os.path.join(processed_data_path, "population", "population_projections.gpkg"),
        layer="mean",
    )
    population.columns = population.columns.map(str)
    population["working_frac"] = (
        population[f"working_{population_year}"] / population[f"{population_year}"]
    )
    buildings = pd.merge(
        buildings,
        population[["ED_ID", "ED", "working_frac"]],
        how="left",
        on=["ED_ID", "ED"],
    )
    buildings["working_population"] = (
        buildings["residential_population"] * buildings["working_frac"]
    )
    del population
    """Assign the closest roads to buildings
    """
    buildings_to_roads = ckdnearest(buildings, nodes[["node_id", "geometry"]])
    print(buildings_to_roads)

    population_threshold = 50
    gdp_threshold = 50000
    nodes_population = (
        buildings_to_roads.groupby(["node_id"])["working_population"]
        .sum()
        .reset_index()
    )
    nodes_population = pd.merge(
        nodes_population, nodes[["node_id", "geometry"]], how="left", on=["node_id"]
    )
    nodes_population = gpd.GeoDataFrame(
        nodes_population, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
    )
    nodes_population.to_file(
        os.path.join(
            results_path,
            "flow_mapping",
            "road_nodes_labour_economic_activity_aggregations.gpkg",
        ),
        layer="working_population",
        driver="GPKG",
    )

    nodes_economic_activity = (
        buildings_to_roads.groupby(["node_id"])["total_GDP"].sum().reset_index()
    )
    nodes_economic_activity = gpd.GeoDataFrame(
        pd.merge(
            nodes_economic_activity,
            nodes[["node_id", "geometry"]],
            how="left",
            on=["node_id"],
        ),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    nodes_economic_activity.to_file(
        os.path.join(
            results_path,
            "flow_mapping",
            "road_nodes_labour_economic_activity_aggregations.gpkg",
        ),
        layer="economic_activity",
        driver="GPKG",
    )

    """Build the radiation model
    """
    nodes_population = nodes_population[
        nodes_population["working_population"] >= population_threshold
    ]
    nodes_economic_activity = nodes_economic_activity[
        nodes_economic_activity["total_GDP"] >= gdp_threshold
    ]
    print(nodes_population)
    print(nodes_economic_activity)
    buffer_distance = 10000  # 10 km distance buffer
    # nodes_population = pd.merge(nodes_population,nodes[["node_id","geometry"]],how="left",on=["node_id"])
    nodes_population["geometry"] = nodes_population.apply(
        lambda x: x.geometry.buffer(buffer_distance), axis=1
    )
    nodes_population.rename(columns={"node_id": "origin_id"}, inplace=True)
    # nodes_population = gpd.GeoDataFrame(nodes_population,
    #                                     geometry="geometry",
    #                                     crs=f"EPSG:{epsg_jamaica}")

    # nodes_economic_activity = gpd.GeoDataFrame(pd.merge(nodes_economic_activity,nodes[["node_id","geometry"]],
    #                                     how="left",on=["node_id"]),
    #                                     geometry="geometry",
    #                                     crs=f"EPSG:{epsg_jamaica}")
    nodes_economic_activity.rename(columns={"node_id": "destination_id"}, inplace=True)
    od_pairs = gpd.sjoin(
        nodes_economic_activity, nodes_population, how="inner", predicate="within"
    ).reset_index()

    flow_paths = network_od_paths_assembly(
        od_pairs[["origin_id", "destination_id", "total_GDP"]], G, "time", "total_GDP"
    )
    print(flow_paths)
    flow_paths = flow_paths[flow_paths["gcost"] <= 1.0]
    print(flow_paths)
    flow_paths = pd.merge(
        flow_paths,
        nodes_population[["origin_id", "working_population"]],
        how="left",
        on=["origin_id"],
    )
    flow_paths.to_csv(
        os.path.join(results_path, "flow_mapping", "labour_to_sectors_flow_paths.csv"),
        index=False,
    )
    flow_radius = flow_paths.groupby(["origin_id"])["total_GDP"].sum().reset_index()
    flow_radius.rename(columns={"total_GDP": "radius_GDP"}, inplace=True)

    flow_paths = pd.merge(flow_paths, flow_radius, how="left", on=["origin_id"])

    flow_paths["t_ij_ext"] = flow_paths.progress_apply(
        lambda x: x["total_GDP"] / (x["radius_GDP"] - x["total_GDP"]), axis=1
    )
    flow_paths_sums = flow_paths.groupby(["origin_id"])["t_ij_ext"].sum().reset_index()
    flow_paths_sums.rename(columns={"t_ij_ext": "t_ij_ext_sums"}, inplace=True)
    flow_paths = pd.merge(
        flow_paths, flow_paths_sums, how="left", on=["origin_id"]
    ).fillna(0)
    del flow_paths_sums
    flow_paths["working_trips"] = flow_paths.progress_apply(
        lambda x: x["working_population"] * (x["t_ij_ext"] / x["t_ij_ext_sums"]), axis=1
    )
    # flow_paths.drop("t_ij_ext_sums",axis=1,inplace=True)
    flow_paths_sums = (
        flow_paths.groupby(["destination_id"])["working_trips"].sum().reset_index()
    )
    flow_paths_sums.rename(
        columns={"working_trips": "working_trips_sums"}, inplace=True
    )
    flow_paths = pd.merge(
        flow_paths, flow_paths_sums, how="left", on=["destination_id"]
    ).fillna(0)
    del flow_paths_sums
    flow_paths["GDP_to_trips"] = flow_paths.progress_apply(
        lambda x: x["total_GDP"] * (x["working_trips"] / x["working_trips_sums"]),
        axis=1,
    )
    flow_paths.drop(["t_ij_ext_sums", "working_trips_sums"], axis=1, inplace=True)
    flow_paths.to_csv(
        os.path.join(
            results_path, "flow_mapping", "labour_to_sectors_trips_and_activity.csv"
        ),
        index=False,
    )

    flow_paths = pd.read_csv(
        os.path.join(
            results_path, "flow_mapping", "labour_to_sectors_trips_and_activity.csv"
        )
    )
    # flow_paths["edge_path"] = flow_paths.progress_apply(lambda x:ast.literal_eval(x.edge_path),axis=1)
    # flow_path_indexes = tf.get_flow_paths_indexes_of_edges(flow_paths,"edge_path")
    # flow_path_counts = tf.get_flow_paths_numbers_of_edges(flow_path_indexes)
    # pd.DataFrame([(k,v) for k,v in flow_path_counts.items()],
    #             columns=["edge_id","no_of_paths"]).to_csv(os.path.join(results_path,
    #                                     'flow_mapping',
    #                                     'roads_flow_paths_count.csv'),index=False)
    common_nodes = flow_paths[flow_paths["origin_id"] == flow_paths["destination_id"]]
    common_nodes = (
        common_nodes.groupby(["origin_id"])[["working_trips", "GDP_to_trips"]]
        .sum()
        .reset_index()
    )

    uncommon_nodes = flow_paths[flow_paths["origin_id"] != flow_paths["destination_id"]]
    origin_trips = (
        uncommon_nodes.groupby(["origin_id"])[["working_trips", "GDP_to_trips"]]
        .sum()
        .reset_index()
    )
    destination_trips = (
        uncommon_nodes.groupby(["destination_id"])[["working_trips", "GDP_to_trips"]]
        .sum()
        .reset_index()
    )
    od_diff = pd.DataFrame(
        list(
            set(
                origin_trips["origin_id"].values.tolist()
                + destination_trips["destination_id"].values.tolist()
            )
        ),
        columns=["node_id"],
    )
    od_diff = pd.merge(
        od_diff,
        origin_trips[["origin_id", "working_trips"]],
        how="left",
        left_on=["node_id"],
        right_on=["origin_id"],
    ).fillna(0)
    od_diff.rename(columns={"working_trips": "o_trip"}, inplace=True)
    od_diff = pd.merge(
        od_diff,
        destination_trips[["destination_id", "working_trips"]],
        how="left",
        left_on=["node_id"],
        right_on=["destination_id"],
    ).fillna(0)
    od_diff.rename(columns={"working_trips": "d_trip"}, inplace=True)

    destination_trips.rename(columns={"destination_id": "origin_id"}, inplace=True)
    node_activity = pd.concat(
        [common_nodes, origin_trips, destination_trips], axis=0, ignore_index=True
    )
    node_activity.rename(columns={"origin_id": "node_id"}, inplace=True)

    node_activity = (
        node_activity.groupby(["node_id"])[["working_trips", "GDP_to_trips"]]
        .sum()
        .reset_index()
    )
    node_activity = pd.merge(
        node_activity,
        od_diff[["node_id", "o_trip", "d_trip"]],
        how="left",
        on=["node_id"],
    ).fillna(0)
    node_activity.to_csv(
        os.path.join(
            results_path,
            "flow_mapping",
            "origins_destinations_labour_economic_activity.csv",
        ),
        index=False,
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
