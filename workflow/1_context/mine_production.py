"""Assign mining GDP values to quarry and bauxite mine areas"""

import logging

import click
import pandas as pd
import geopandas as gpd
import igraph as ig

from jamaica_infrastructure.transport.utils import (
    map_nearest_locations_and_create_lines,
    network_od_paths_assembly,
)
from jamaica_infrastructure.geo import LOCAL_PROJ_CRS_EPSG


def route_mines_to_nearest_ports(areas, nodes, edges, ports, connection_type="mines"):
    nearest_stations = map_nearest_locations_and_create_lines(
        areas.copy(),
        nodes[nodes["mode"] == "rail"].copy(),
        "mining_id",
        "node_id",
        connection_type,
        "rail",
    )
    nearest_roads = map_nearest_locations_and_create_lines(
        areas.copy(),
        nodes[nodes["mode"] == "road"].copy(),
        "mining_id",
        "node_id",
        connection_type,
        "road",
    )
    nearest_stations["edge_id"] = nearest_stations.apply(
        lambda x: f"{connection_type}raile_{x.name}", axis=1
    )
    nearest_stations["speed"] = 10.0
    nearest_stations["time"] = (
        0.001 * nearest_stations["length_m"] / nearest_stations["speed"]
    )
    nearest_roads["edge_id"] = nearest_roads.apply(
        lambda x: f"{connection_type}roade_{x.name}", axis=1
    )
    nearest_roads["speed"] = 10.0
    nearest_roads["time"] = 0.001 * nearest_roads["length_m"] / nearest_roads["speed"]

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
    network = pd.concat(
        [edges[columns], nearest_roads[columns], nearest_stations[columns]],
        axis=0,
        ignore_index=True,
    )[columns].copy()
    logging.debug(network)
    # network = network[network["to_mode"] != "road"]
    G = ig.Graph.TupleList(
        network.itertuples(index=False), edge_attrs=list(network.columns)[2:]
    )

    all_ports = ports["node_id"].values.tolist()

    od_pairs = [
        list(zip([b] * len(all_ports), all_ports))
        for b in areas.mining_id.values.tolist()
    ]
    od_pairs = [item for sublist in od_pairs for item in sublist]
    od_pairs = pd.DataFrame(od_pairs, columns=["origin_id", "destination_id"])
    od_pairs = pd.merge(
        od_pairs,
        areas[["mining_id", "area_m2"]],
        how="left",
        left_on=["origin_id"],
        right_on=["mining_id"],
    )
    flow_paths = network_od_paths_assembly(
        od_pairs[["origin_id", "destination_id", "area_m2"]], G, "time", "area_m2"
    )
    flow_paths = flow_paths.sort_values(by="gcost")
    flow_paths = flow_paths.drop_duplicates(subset=["origin_id"], keep="first")

    return (
        flow_paths,
        pd.concat(
            [nearest_roads[columns], nearest_stations[columns]],
            axis=0,
            ignore_index=True,
        )[columns],
    )


@click.command()
@click.version_option("1.0")
@click.option(
    "--land_use_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to land use GeoPackage.",
)
@click.option(
    "--ports_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to ports GeoPackage.",
)
@click.option(
    "--network_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to multimodal network GeoPackage.",
)
@click.option(
    "--economic_output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to macroeconomic CSV file.",
)
@click.option(
    "--output_path",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Path to output mining GDP GeoPackage.",
)
def main(land_use_path, ports_path, network_path, economic_output_path, output_path):
    """Assign GDP to mining areas, separately for quarry and bauxite, in
    proportion to export volume of nearest ports, then in proportion to mining
    area footprint.

    Example usage:
        python workflow/1_context/mine_production.py \
            --land_use_path processed_data/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg \
            --ports_path processed_data/networks/transport/port_polygon.gpkg \
            --network_path processed_data/networks/transport/multi_modal_network.gpkg \
            --economic_output_path processed_data/macroeconomic_data/NIP_2023.csv \
            --output_path processed_data/mining_data/mining_gdp.gpkg
    """
    land_use_types = gpd.read_file(land_use_path, layer="areas")
    land_use_types = land_use_types.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)

    if "index_left" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_left", axis=1, inplace=True)
    elif "index_right" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_right", axis=1, inplace=True)

    logging.debug(land_use_types.columns)

    mining_areas = land_use_types[
        (land_use_types["sector_code_forest"] == "C")
        | (land_use_types["sector_code_tnc"] == "C")
        | (land_use_types["global_LU_type"] == "Bauxite Extraction")
    ].copy()
    mining_areas["mining_id"] = mining_areas.index.values.tolist()
    mining_areas["mining_id"] = mining_areas["mining_id"].astype(str)
    logging.debug(mining_areas)

    quarry_areas = mining_areas[
        (mining_areas["subsector_code_forest"].isin([141, "141"]))
        | (mining_areas["subsector_code_tnc"].isin([141, "141"]))
    ].copy()
    logging.debug(quarry_areas)
    bauxite_areas = mining_areas[
        ~mining_areas["mining_id"].isin(quarry_areas["mining_id"].values.tolist())
    ].copy()
    logging.debug(bauxite_areas)

    quarry_areas["geometry"] = quarry_areas.apply(lambda x: x.geometry.centroid, axis=1)
    quarry_areas = quarry_areas.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    logging.debug(quarry_areas)

    bauxite_areas["geometry"] = bauxite_areas.apply(
        lambda x: x.geometry.centroid, axis=1
    )
    bauxite_areas = bauxite_areas.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    logging.debug(bauxite_areas)

    del land_use_types

    """Get the ports and multimodal network
    """
    ports = gpd.read_file(ports_path, layer="areas")
    ports = ports[ports["commodity"] == "alumina"]
    ports = ports.drop_duplicates(subset=["name"], keep="first")
    ports = ports[["node_id", "name", "export_tonnes", "geometry"]].copy()
    ports["export_wt"] = ports["export_tonnes"] / ports["export_tonnes"].sum()
    ports["geometry"] = ports.apply(lambda x: x.geometry.centroid, axis=1)
    ports = ports.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    logging.debug(ports)

    nodes = gpd.read_file(network_path, layer="nodes")
    nodes = nodes.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    edges = gpd.read_file(network_path, layer="edges")
    multi_edges = edges[edges["from_mode"] != edges["to_mode"]].copy()
    multi_edges = multi_edges[multi_edges["length_m"] >= 5000][
        "edge_id"
    ].values.tolist()
    edges = edges[~edges["edge_id"].isin(multi_edges)].copy()

    """Find routes to quarries and mines to their closest to ports to get total areas proportioned to ports
    """

    quarry_to_ports, quarry_network = route_mines_to_nearest_ports(
        quarry_areas, nodes, edges, ports, connection_type="quarry"
    )
    bauxite_to_ports, bauxite_network = route_mines_to_nearest_ports(
        bauxite_areas, nodes, edges, ports
    )

    quarry_to_ports_areas = (
        quarry_to_ports.groupby(["destination_id"])["area_m2"].sum().reset_index()
    )
    quarry_to_ports_areas.rename(columns={"area_m2": "tot_area"}, inplace=True)
    bauxite_to_ports_areas = (
        bauxite_to_ports.groupby(["destination_id"])["area_m2"].sum().reset_index()
    )
    bauxite_to_ports_areas.rename(columns={"area_m2": "tot_area"}, inplace=True)
    quarry_to_ports = pd.merge(
        quarry_to_ports, quarry_to_ports_areas, how="left", on=["destination_id"]
    )
    bauxite_to_ports = pd.merge(
        bauxite_to_ports, bauxite_to_ports_areas, how="left", on=["destination_id"]
    )
    del quarry_to_ports_areas, bauxite_to_ports_areas

    #
    # Allocate the GDP values to the mining areas
    #
    economic_output = pd.read_csv(economic_output_path, comment="#")

    quarry_output = economic_output[economic_output["Subsector"] == "Quarrying"][
        "GVA JMD Millions (2023)"
    ].sum()
    bauxite_output = economic_output[
        economic_output["Subsector"] == "Bauxite & Alumina"
    ]["GVA JMD Millions (2023)"].sum()

    logging.info("Given GDP %f", quarry_output + bauxite_output)
    quarry_to_ports = pd.merge(
        quarry_to_ports,
        ports[["node_id", "export_wt"]],
        how="left",
        left_on=["destination_id"],
        right_on=["node_id"],
    )
    quarry_to_ports["C_GDP"] = (
        quarry_output
        * quarry_to_ports["export_wt"]
        * quarry_to_ports["area_m2"]
        / quarry_to_ports["tot_area"]
    )
    bauxite_to_ports = pd.merge(
        bauxite_to_ports,
        ports[["node_id", "export_wt"]],
        how="left",
        left_on=["destination_id"],
        right_on=["node_id"],
    )
    bauxite_to_ports["C_GDP"] = (
        bauxite_output
        * bauxite_to_ports["export_wt"]
        * bauxite_to_ports["area_m2"]
        / bauxite_to_ports["tot_area"]
    )
    quarry_to_ports["GDP_persqm"] = (
        quarry_to_ports["C_GDP"] / quarry_to_ports["area_m2"]
    )
    bauxite_to_ports["GDP_persqm"] = (
        bauxite_to_ports["C_GDP"] / bauxite_to_ports["area_m2"]
    )
    quarry_to_ports["mining_class"] = "quarry"
    bauxite_to_ports["mining_class"] = "bauxite"

    cols = ["origin_id", "C_GDP", "GDP_persqm", "mining_class"]
    mining_gdp = pd.concat(
        [
            quarry_to_ports[cols],
            bauxite_to_ports[cols],
        ],
        axis=0,
        ignore_index=True,
    )
    mining_gdp.rename(columns={"origin_id": "mining_id"}, inplace=True)
    mining_gdp["GDP_unit"] = "JD/day"
    logging.debug(mining_gdp)

    cols = ["mining_id", "area_m2", "geometry"]
    mining_areas = pd.merge(
        mining_areas[cols], mining_gdp, how="left", on=["mining_id"]
    )
    mining_areas = gpd.GeoDataFrame(
        mining_areas, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    )

    tot_gpd = mining_areas["C_GDP"].sum()
    tot_area = mining_areas["area_m2"].sum()
    logging.info("Estimated GDP %f", tot_gpd)
    logging.info("Estimated Areas %f", tot_area)

    mining_areas.to_file(output_path, layer="areas", engine="pyogrio")


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    logging.info("Start %s", __file__)
    main()
