"""Map mining flows onto the rail network of Jamaica
"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
import igraph as ig
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import (
    get_flow_on_edges,
    load_config,
    map_nearest_locations_and_create_lines,
    network_od_paths_assembly,
)

tqdm.pandas()
epsg_jamaica = 3448


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
    nearest_stations["edge_id"] = nearest_stations.progress_apply(
        lambda x: f"{connection_type}raile_{x.name}", axis=1
    )
    nearest_stations["speed"] = 10.0
    nearest_stations["time"] = (
        0.001 * nearest_stations["length_m"] / nearest_stations["speed"]
    )
    nearest_roads["edge_id"] = nearest_roads.progress_apply(
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
    )[columns]
    print(network)
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


def main(config):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]
    """Get the mining areas
    """
    land_use_types = gpd.read_file(
        os.path.join(
            processed_data_path,
            "land_type_and_use",
            "jamaica_land_use_combined_with_sectors.gpkg",
        ),
        layer="areas",
    )
    land_use_types = land_use_types.to_crs(epsg=epsg_jamaica)

    if "index_left" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_left", axis=1, inplace=True)
    elif "index_right" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_right", axis=1, inplace=True)

    print(land_use_types.columns)

    mining_areas = land_use_types[
        (land_use_types["sector_code_forest"] == "C")
        | (land_use_types["sector_code_tnc"] == "C")
        | (land_use_types["global_LU_type"] == "Bauxite Extraction")
    ]
    mining_areas["mining_id"] = mining_areas.index.values.tolist()
    mining_areas["mining_id"] = mining_areas["mining_id"].astype(str)
    print(mining_areas)

    quarry_areas = mining_areas[
        (mining_areas["subsector_code_forest"].isin([141, "141"]))
        | (mining_areas["subsector_code_tnc"].isin([141, "141"]))
    ]
    print(quarry_areas)
    bauxite_areas = mining_areas[
        ~mining_areas["mining_id"].isin(quarry_areas["mining_id"].values.tolist())
    ]
    print(bauxite_areas)

    quarry_areas["geometry"] = quarry_areas.progress_apply(
        lambda x: x.geometry.centroid, axis=1
    )
    quarry_areas = quarry_areas.to_crs(epsg=epsg_jamaica)
    print(quarry_areas)

    bauxite_areas["geometry"] = bauxite_areas.progress_apply(
        lambda x: x.geometry.centroid, axis=1
    )
    bauxite_areas = bauxite_areas.to_crs(epsg=epsg_jamaica)
    print(bauxite_areas)

    del land_use_types

    """Get the ports and multimodal network
    """
    ports = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "port_polygon.gpkg"),
        layer="areas",
    )
    ports = ports[ports["commodity"] == "alumina"]
    ports = ports.drop_duplicates(subset=["name"], keep="first")
    ports = ports[["node_id", "name", "export_tonnes", "geometry"]]
    ports["export_wt"] = ports["export_tonnes"] / ports["export_tonnes"].sum()
    ports["geometry"] = ports.progress_apply(lambda x: x.geometry.centroid, axis=1)
    ports = ports.to_crs(epsg=epsg_jamaica)
    print(ports)

    nodes = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="nodes",
    )
    nodes = nodes.to_crs(epsg=epsg_jamaica)
    edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    multi_edges = edges[edges["from_mode"] != edges["to_mode"]]
    multi_edges = multi_edges[multi_edges["length_m"] >= 5000][
        "edge_id"
    ].values.tolist()
    edges = edges[~edges["edge_id"].isin(multi_edges)]

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

    """Allocate the GDP values to the mining areas
    """
    financial_year = 2019
    economic_output = pd.read_excel(
        os.path.join(
            processed_data_path,
            "macroeconomic_data",
            "detailed_sector_GVA_GDP_current_prices.xlsx",
        ),
        sheet_name="2019",
    )
    economic_output.columns = [
        str(c).strip() for c in economic_output.columns.values.tolist()
    ]
    economic_output["subsector_code"] = economic_output["subsector_code"].apply(str)
    totat_gva = economic_output[economic_output["sector_code"] == "GVA"][
        f"{financial_year}"
    ].sum()
    total_tax = economic_output[economic_output["sector_code"] == "TAX"][
        f"{financial_year}"
    ].sum()
    tax_rate = 1.0 * total_tax / totat_gva

    quarry_output = (
        (1 + tax_rate)
        * (1.0e6 / 365.0)
        * economic_output[
            (economic_output["sector_code"] == "C")
            & (economic_output["subsector_code"] == "141")
        ][f"{financial_year}"].sum()
    )

    bauxite_output = (
        (1 + tax_rate)
        * (1.0e6 / 365.0)
        * economic_output[
            (economic_output["sector_code"] == "C")
            & (economic_output["subsector_code"] == "132")
        ][f"{financial_year}"].sum()
    )
    print("* Given GDP", quarry_output + bauxite_output)
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

    mining_gdp = pd.concat(
        [
            quarry_to_ports[["origin_id", "C_GDP", "GDP_persqm"]],
            bauxite_to_ports[["origin_id", "C_GDP", "GDP_persqm"]],
        ],
        axis=0,
        ignore_index=True,
    )
    mining_gdp.rename(columns={"origin_id": "mining_id"}, inplace=True)
    mining_gdp["GDP_unit"] = "JD/day"
    print(mining_gdp)

    mining_areas = pd.merge(mining_areas, mining_gdp, how="left", on=["mining_id"])
    mining_areas = gpd.GeoDataFrame(
        mining_areas, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
    )

    tot_gpd = mining_areas["C_GDP"].sum()
    tot_area = mining_areas["area_m2"].sum()
    print("* Estimated GDP", tot_gpd)
    print("* Estimated Areas", tot_area)

    mining_areas.to_file(
        os.path.join(processed_data_path, "mining_data", "mining_gdp.gpkg"),
        layer="areas",
        driver="GPKG",
    )
    flow_paths = pd.concat(
        [quarry_to_ports, bauxite_to_ports], axis=0, ignore_index=True
    )
    flow_paths[["origin_id", "destination_id", "edge_path", "gcost", "C_GDP"]].to_csv(
        os.path.join(results_path, "flow_mapping", "mines_flow_paths.csv"), index=False
    )
    edge_flows = get_flow_on_edges(flow_paths, "edge_id", "edge_path", "C_GDP")
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
        [edges[columns], quarry_network[columns], bauxite_network[columns]],
        axis=0,
        ignore_index=True,
    )[columns]
    edge_flows = gpd.GeoDataFrame(
        pd.merge(edge_flows, network, how="left", on=["edge_id"]),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    edge_flows.to_file(
        os.path.join(results_path, "flow_mapping", "mines_flows.gpkg"),
        layer="edges",
        driver="GPKG",
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
