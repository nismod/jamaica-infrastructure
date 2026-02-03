import logging
import os

import click
import pandas as pd
import geopandas as gpd
import igraph as ig
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import (
    ckdnearest,
    get_flow_on_edges,
    network_od_paths_assembly,
)
from jamaica_infrastructure.geo import LOCAL_PROJ_CRS_EPSG


tqdm.pandas()


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--network-file",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Jamaican Multimodal Network GeoPackage file",
)
@click.option(
    "--buildings-file",
    "-b",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the building economic activity GeoPackage file.",
)
@click.option(
    "--population-file",
    "-p",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the population projections GeoPackage file.",
)
@click.option(
    "--out-dir",
    "-o",
    required=True,
    help="Directory to save the results.",
)
def commuter_flow_mapping(network_file, buildings_file, population_file, out_dir):
    """Map labour to work trips as a proxy to map GDP onto roads"""

    logging.info("Read input data")
    nodes = gpd.read_file(network_file, layer="nodes")
    edges = gpd.read_file(network_file, layer="edges")
    buildings = gpd.read_file(buildings_file, layer="areas")
    population = gpd.read_file(population_file, layer="mean")

    # Set up variables
    nodes = nodes[nodes["mode"] == "road"]
    nodes = nodes.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)

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
    network = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")][columns]

    logging.info("Build network graph")
    graph = ig.Graph.TupleList(network.itertuples(index=False), edge_attrs=list(network.columns)[2:])

    buildings["osm_id"] = buildings.progress_apply(lambda x: f"building_{x.osm_id}", axis=1)
    buildings["geometry"] = buildings.progress_apply(lambda x: x.geometry.centroid, axis=1)
    buildings = buildings.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)

    population_year = 2019
    population.columns = population.columns.map(str)
    population["working_frac"] = population[f"working_{population_year}"] / population[f"{population_year}"]

    logging.info("Assign buildings to road nodes")
    buildings = pd.merge(
        buildings,
        population[["ED_ID", "ED", "working_frac"]],
        how="left",
        on=["ED_ID", "ED"],
    )
    buildings["working_population"] = buildings["residential_population"] * buildings["working_frac"]

    # Assign the closest roads to buildings
    buildings_to_roads = ckdnearest(buildings, nodes[["node_id", "geometry"]])

    population_threshold = 50
    gdp_threshold = 50000
    nodes_population = buildings_to_roads.groupby(["node_id"])["working_population"].sum().reset_index()
    nodes_population = pd.merge(nodes_population, nodes[["node_id", "geometry"]], how="left", on=["node_id"])
    nodes_population = gpd.GeoDataFrame(nodes_population, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}")

    logging.info("Write out working population at road nodes")
    nodes_population.to_file(
        os.path.join(out_dir, "road_nodes_labour_economic_activity_aggregations.gpkg"),
        layer="working_population",
        driver="GPKG",
    )

    nodes_economic_activity = buildings_to_roads.groupby(["node_id"])["total_GDP"].sum().reset_index()
    nodes_economic_activity = gpd.GeoDataFrame(
        pd.merge(
            nodes_economic_activity,
            nodes[["node_id", "geometry"]],
            how="left",
            on=["node_id"],
        ),
        geometry="geometry",
        crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}",
    )
    logging.info("Write out economic activity at road nodes")
    nodes_economic_activity.to_file(
        os.path.join(out_dir, "road_nodes_labour_economic_activity_aggregations.gpkg"),
        layer="economic_activity",
        driver="GPKG",
    )

    logging.info("Build radiation model")
    nodes_population = nodes_population[nodes_population["working_population"] >= population_threshold]
    nodes_economic_activity = nodes_economic_activity[nodes_economic_activity["total_GDP"] >= gdp_threshold]
    buffer_distance = 10000  # 10 km distance buffer
    logging.info("Buffer nodes")
    nodes_population["geometry"] = nodes_population.apply(lambda x: x.geometry.buffer(buffer_distance), axis=1)
    nodes_population.rename(columns={"node_id": "origin_id"}, inplace=True)
    nodes_economic_activity.rename(columns={"node_id": "destination_id"}, inplace=True)
    logging.info("Spatially join economic activity proximate population (create OD)")
    od_pairs = gpd.sjoin(nodes_economic_activity, nodes_population, how="inner", predicate="within").reset_index()

    logging.info("Find shortest paths for OD over network")
    flow_paths = network_od_paths_assembly(od_pairs[["origin_id", "destination_id", "total_GDP"]], graph, "time", "total_GDP")
    flow_paths = flow_paths[flow_paths["gcost"] <= 1.0]
    flow_paths = pd.merge(
        flow_paths,
        nodes_population[["origin_id", "working_population"]],
        how="left",
        on=["origin_id"],
    )
    logging.info("Write out OD and shortest paths to disk")
    flow_paths.to_parquet(
        os.path.join(out_dir, "labour_to_sectors_flow_paths.pq"),
        index=False,
    )

    logging.info("Sum up flows")
    flow_radius = flow_paths.groupby(["origin_id"])["total_GDP"].sum().reset_index()
    flow_radius.rename(columns={"total_GDP": "radius_GDP"}, inplace=True)

    flow_paths = pd.merge(flow_paths, flow_radius, how="left", on=["origin_id"])

    flow_paths["t_ij_ext"] = flow_paths.progress_apply(lambda x: x["total_GDP"] / (x["radius_GDP"] - x["total_GDP"]), axis=1)
    flow_paths_sums = flow_paths.groupby(["origin_id"])["t_ij_ext"].sum().reset_index()
    flow_paths_sums.rename(columns={"t_ij_ext": "t_ij_ext_sums"}, inplace=True)
    flow_paths = pd.merge(flow_paths, flow_paths_sums, how="left", on=["origin_id"]).fillna(0)

    flow_paths["working_trips"] = flow_paths.progress_apply(
        lambda x: x["working_population"] * (x["t_ij_ext"] / x["t_ij_ext_sums"]), axis=1
    )
    # flow_paths.drop("t_ij_ext_sums",axis=1,inplace=True)
    flow_paths_sums_wt = flow_paths.groupby(["destination_id"])["working_trips"].sum().reset_index()
    flow_paths_sums_wt.rename(columns={"working_trips": "working_trips_sums"}, inplace=True)
    flow_paths = pd.merge(flow_paths, flow_paths_sums_wt, how="left", on=["destination_id"]).fillna(0)

    flow_paths["GDP_to_trips"] = flow_paths.progress_apply(
        lambda x: x["total_GDP"] * (x["working_trips"] / x["working_trips_sums"]),
        axis=1,
    )
    flow_paths.drop(["t_ij_ext_sums", "working_trips_sums"], axis=1, inplace=True)

    logging.info("Write out flows to disk")
    flow_paths.to_parquet(
        os.path.join(out_dir, "labour_to_sectors_trips_and_activity.pq"),
        index=False,
    )

    common_nodes = flow_paths[flow_paths["origin_id"] == flow_paths["destination_id"]]
    common_nodes = common_nodes.groupby(["origin_id"])[["working_trips", "GDP_to_trips"]].sum().reset_index()

    uncommon_nodes = flow_paths[flow_paths["origin_id"] != flow_paths["destination_id"]]
    origin_trips = uncommon_nodes.groupby(["origin_id"])[["working_trips", "GDP_to_trips"]].sum().reset_index()
    destination_trips = uncommon_nodes.groupby(["destination_id"])[["working_trips", "GDP_to_trips"]].sum().reset_index()
    od_diff = pd.DataFrame(
        list(set(origin_trips["origin_id"].values.tolist() + destination_trips["destination_id"].values.tolist())),
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
    node_activity = pd.concat([common_nodes, origin_trips, destination_trips], axis=0, ignore_index=True)
    node_activity.rename(columns={"origin_id": "node_id"}, inplace=True)

    logging.info("Write out nodal economic activity")
    node_activity = node_activity.groupby(["node_id"])[["working_trips", "GDP_to_trips"]].sum().reset_index()
    node_activity = pd.merge(
        node_activity,
        od_diff[["node_id", "o_trip", "d_trip"]],
        how="left",
        on=["node_id"],
    ).fillna(0)
    node_activity.to_csv(
        os.path.join(out_dir, "origins_destinations_labour_economic_activity.csv"),
        index=False,
    )

    logging.info("Accumulate flows to edges")
    flow_paths = flow_paths[flow_paths["working_trips"] >= 1]
    edge_flows_trips = get_flow_on_edges(flow_paths, "edge_id", "edge_path", "working_trips")
    edge_flows_gdp = get_flow_on_edges(flow_paths, "edge_id", "edge_path", "GDP_to_trips")
    network = network.merge(edge_flows_trips, how="left", on=["edge_id"]) \
        .merge(edge_flows_gdp, how="left", on=["edge_id"]).fillna(0)

    logging.info("Write out accumulated flows on edges")
    network.to_parquet(os.path.join(out_dir, "labour_trips_and_activity.gpq"))


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    commuter_flow_mapping()
