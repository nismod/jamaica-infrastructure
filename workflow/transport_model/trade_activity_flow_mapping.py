"""Map import and export trade activites to Ports in Jamaica
"""

import os

import pandas as pd
import geopandas as gpd
import numpy as np
import igraph as ig
from jamaica_infrastructure.utils import (
    get_flow_on_edges,
    load_config,
    map_nearest_locations_and_create_lines,
    network_od_paths_assembly,
)
from tqdm import tqdm

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
        nearest_stations = nearest_stations[nearest_stations["length_m"] < 5000]

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
    config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    results_path = processed_data_path

    """Get the ports and multimodal network
    """
    jam_ports = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "port_polygon.gpkg"),
        layer="areas",
    )
    jam_ports["export_tonnes"] = jam_ports["export_tonnes"].fillna(0)
    jam_ports["import_tonnes"] = jam_ports["import_tonnes"].fillna(0)

    mining_ports = port_import_exports(
        jam_ports[jam_ports["commodity"] == "alumina"], "export_tonnes", "mining"
    )
    non_mining_export_ports = port_import_exports(
        jam_ports[
            (jam_ports["commodity"] != "alumina") & (jam_ports["export_tonnes"] > 0)
        ],
        "export_tonnes",
        "export",
    )
    non_mining_import_ports = port_import_exports(
        jam_ports[
            (jam_ports["commodity"] != "alumina") & (jam_ports["import_tonnes"] > 0)
        ],
        "import_tonnes",
        "import",
    )
    petrojam_port = port_import_exports(
        jam_ports[jam_ports["name"] == "Petrojam"], "import_tonnes", "petrojam"
    )

    all_ports = pd.concat(
        [mining_ports, non_mining_export_ports, non_mining_import_ports, petrojam_port],
        axis=0,
        ignore_index=True,
    ).fillna(0)
    del mining_ports, non_mining_export_ports, non_mining_import_ports, petrojam_port
    print(all_ports)

    nodes = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="nodes",
    )
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
    multi_edges = edges[edges["from_mode"] != edges["to_mode"]]
    multi_edges = multi_edges[multi_edges["length_m"] >= 5000][
        "edge_id"
    ].values.tolist()
    edges = edges[~edges["edge_id"].isin(multi_edges)]
    del multi_edges
    print(edges)

    """Read the export datasets and create the specific cases of exports
        The exporting sectors are Agriculture, Fisheries, Mining & Quarrying, Manufacturing, Fuels
        For agriculture we consider all the areas of crops
        For Mining & Quarrying we need to consider the separation between Bauxite/Alumina and Rest of Quarrying
        For Manufacturing we need to consider every sector, except the sector D-240 which is for Fuels
        Fuels are pretty much exported from the Petrojam facility    
    """
    financial_year = 2019
    exports = pd.read_excel(
        os.path.join(
            processed_data_path, "macroeconomic_data", "domestic_export_by_sector.xlsx"
        ),
        sheet_name=str(financial_year),
    )
    exports.columns = [str(c).strip() for c in exports.columns.values.tolist()]
    exports = exports[~exports["subsector_code"].isin(["ALL", "SUB"])]
    exports = (
        exports.groupby(["sector_code", "sector_type"])[str(financial_year)]
        .sum()
        .reset_index()
    )
    exports = exports[
        exports["sector_code"].isin(["A", "B", "C", "D"])
    ]  # These are mostly the sectors exporting commodities
    exports["trade_value"] = 1.0e3 / 365.0 * exports[str(financial_year)]
    exports["trade_type"] = "export"
    """Read the import datasets and create the specific cases of imports
        The importing sectors are Retail and Trade, Construction, Fuels
        For Retail and Trade we consider all the buildings assigned to th sector/subsector - G/500-1
        Fuels are being imported 
            Into the Petrojam facility, 
            Into the different ports that then connects to the mining locations - C/132
            Towards the different locations belonging to the automobile trade sector - G/500-2/3  
    """
    imports = pd.read_excel(
        os.path.join(
            processed_data_path, "macroeconomic_data", "import_by_industry.xlsx"
        ),
        sheet_name=str(financial_year),
    )
    fuel_shares = pd.read_excel(
        os.path.join(
            processed_data_path, "macroeconomic_data", "import_by_industry.xlsx"
        ),
        sheet_name="fuel_shares",
    )
    imports.columns = [str(c).strip() for c in imports.columns.values.tolist()]
    imports = imports[~imports["subsector_code"].isin(["ALL", "SUB"])]
    imports_groups = (
        imports.groupby(["sector_code", "sector_type"])[f"{financial_year}"]
        .sum()
        .reset_index()
    )
    imports_fuels = imports_groups[imports_groups["sector_code"] == "C,D,G"][
        f"{financial_year}"
    ].sum()
    fuel_shares[str(financial_year)] = imports_fuels * fuel_shares["weight"]
    imports_groups = imports_groups[imports_groups["sector_code"] != "C,D,G"]
    imports = pd.concat(
        [
            imports_groups,
            fuel_shares[["sector_code", "sector_type", str(financial_year)]],
        ],
        axis=0,
        ignore_index=True,
    )
    imports = (
        imports.groupby(["sector_code", "sector_type"])[f"{financial_year}"]
        .sum()
        .reset_index()
    )
    imports["trade_value"] = 1.0e3 / 365.0 * imports[str(financial_year)]
    imports["trade_type"] = "import"

    trade = pd.concat([imports, exports], axis=0, ignore_index=True)
    # trade.columns = [str(c).strip() for c in trade.columns.values.tolist()]
    print(trade)

    """Find routes of economic activites to their closest to ports to get total areas proportioned to ports
    """
    buildings = gpd.read_file(
        os.path.join(
            processed_data_path,
            "buildings",
            "buildings_assigned_economic_activity.gpkg",
        ),
        layer="areas",
    )
    id_column = "osm_id"
    buildings["osm_id"] = buildings.progress_apply(
        lambda x: f"building_{x.osm_id}", axis=1
    )

    sector_network = [edges[columns]]
    flow_network = []
    sector_flows = []
    for trade_details in trade.itertuples():
        # print (trade_details)
        include_rail = False
        id_column = "osm_id"
        port_wt = f"{trade_details.trade_type}_wt"
        gdp_value = f"{trade_details.sector_code}_GDP"
        if trade_details.sector_type == "Agriculture":
            gdp_areas = gpd.read_file(
                os.path.join(
                    processed_data_path, "agriculture_data", "agriculture_gdp.gpkg"
                ),
                layer="areas",
            )
            gdp_areas = gdp_areas[gdp_areas["crop_tons"] > 0]
            id_column = "land_id"
            port_wt = f"{trade_details.trade_type}_wt"
            # ports = all_ports[all_ports[f"{trade_details.trade_type}_wt"] > 0]
        elif trade_details.sector_type in ("Mining", "Quarrying", "Fuel for mines"):
            mining_areas = gpd.read_file(
                os.path.join(processed_data_path, "mining_data", "mining_gdp.gpkg"),
                layer="areas",
            )
            id_column = "mining_id"
            # mining_areas[id_column] = mining_areas[id_column].astype(str)
            mining_areas["mining_id"] = mining_areas.progress_apply(
                lambda x: f"mines_{x.mining_id}", axis=1
            )
            quarry_areas = mining_areas[
                (mining_areas["subsector_code_forest"].isin([141, "141"]))
                | (mining_areas["subsector_code_tnc"].isin([141, "141"]))
            ]
            if trade_details.sector_type == "Quarrying":
                gdp_areas = quarry_areas.copy()
            else:
                gdp_areas = mining_areas[
                    ~mining_areas["mining_id"].isin(
                        quarry_areas["mining_id"].values.tolist()
                    )
                ].copy()

            # ports = all_ports[all_ports["mining_wt"] > 0]
            port_wt = "mining_wt"
            include_rail = True

        elif trade_details.sector_type in ("Equipment", "Fuel for parts"):
            gdp_areas = filter_sector_from_buildings(buildings, "G", "500-2/3")
            # ports = all_ports[all_ports[f"{trade_details.trade_type}_wt"] > 0]
        elif trade_details.sector_type == "Retail":
            gdp_areas = filter_sector_from_buildings(buildings, "G", "500-1")
            # ports = all_ports[all_ports[f"{trade_details.trade_type}_wt"] > 0]
        elif trade_details.sector_type == "Petrojam":
            gdp_areas = jam_ports[jam_ports["name"] == "Petrojam"]
            id_column = "node_id"
            port_wt = "petrojam_wt"
            gdp_value = "import_tonnes"
            # ports = all_ports[all_ports["petrojam_port"] > 0]
        else:
            gdp_areas = buildings[buildings[f"{trade_details.sector_code}_GDP"] > 0]

        gdp_areas["geometry"] = gdp_areas.progress_apply(
            lambda x: x.geometry.centroid, axis=1
        )
        gdp_areas = gdp_areas.to_crs(epsg=epsg_jamaica)

        print(gdp_areas)
        print(all_ports[all_ports[port_wt] > 0])

        sector_to_ports, sector_edges = route_areas_to_nearest_ports(
            gdp_areas,
            id_column,
            gdp_value,
            nodes.copy(),
            edges.copy(),
            all_ports[all_ports[port_wt] > 0],
            port_wt,
            connection_type=trade_details.sector_code,
            trade_type=trade_details.trade_type,
            include_rail=include_rail,
        )
        sector_to_ports[f"{trade_details.sector_code}_trade"] = (
            trade_details.trade_value * sector_to_ports["trade_wt"]
        )
        sector_to_ports["trade_type"] = trade_details.trade_type
        sector_to_ports["sector_type"] = trade_details.sector_type
        sector_flows.append(sector_to_ports)
        # sector_to_ports.to_csv(os.path.join(results_path,
        #                                 'flow_mapping',
        #                                 f'sector_{sector}_{trade_type}_to_ports_flow_paths.csv'),
        #                         index=False)
        edge_flows = get_flow_on_edges(
            sector_to_ports,
            "edge_id",
            "edge_path",
            f"{trade_details.sector_code}_trade",
        )
        sector_network.append(sector_edges[columns])
        flow_network.append(edge_flows)

        print(f"* Done with {trade_details.sector_type} {trade_details.trade_type}")
    pd.concat(sector_flows, axis=0, ignore_index=True).to_csv(
        os.path.join(results_path, "flow_mapping", f"sector_to_ports_flow_paths.csv"),
        index=False,
    )
    sector_network = pd.concat(sector_network, axis=0, ignore_index=True)[columns]
    flow_network = pd.concat(flow_network, axis=0, ignore_index=True)
    flow_gdp_columns = [
        c for c in flow_network.columns.values.tolist() if "_trade" in c
    ]
    for f in flow_gdp_columns:
        flow_network[f] = flow_network[f].fillna(0)
    flow_index_columns = [
        c for c in flow_network.columns.values.tolist() if c not in flow_gdp_columns
    ]
    flow_network = (
        flow_network.groupby(flow_index_columns)[flow_gdp_columns].sum().reset_index()
    )
    flow_network["total_trade"] = flow_network[flow_gdp_columns].sum(axis=1)

    edge_flows = gpd.GeoDataFrame(
        pd.merge(sector_network, flow_network, how="left", on=["edge_id"]),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    ).fillna(0)
    # edge_flows["total_trade"] = edge_flows[[f"{s}_trade" for s in list(set([sc["sector"] for sc in sector_details]))]].sum(axis=1)
    edge_flows.to_file(
        os.path.join(
            results_path, "flow_mapping", "sector_imports_exports_to_ports_flows.gpkg"
        ),
        layer="edges",
        driver="GPKG",
    )

    flow_paths = pd.read_csv(
        os.path.join(results_path, "flow_mapping", "sector_to_ports_flow_paths.csv")
    ).fillna(0)
    trade_columns = [c for c in flow_paths.columns.values.tolist() if "_trade" in c]
    common_nodes = flow_paths[flow_paths["origin_id"] == flow_paths["destination_id"]]
    common_nodes = (
        common_nodes.groupby(["origin_id"])[trade_columns].sum().reset_index()
    )

    uncommon_nodes = flow_paths[flow_paths["origin_id"] != flow_paths["destination_id"]]
    origin_trips = (
        uncommon_nodes.groupby(["origin_id"])[trade_columns].sum().reset_index()
    )
    destination_trips = (
        uncommon_nodes.groupby(["destination_id"])[trade_columns].sum().reset_index()
    )
    destination_trips.rename(columns={"destination_id": "origin_id"}, inplace=True)
    node_activity = pd.concat(
        [common_nodes, origin_trips, destination_trips], axis=0, ignore_index=True
    )
    node_activity.rename(columns={"origin_id": "node_id"}, inplace=True)

    node_activity = (
        node_activity.groupby(["node_id"])[trade_columns].sum().reset_index()
    )
    node_activity["total_trade"] = node_activity[trade_columns].sum(axis=1)
    node_activity.to_csv(
        os.path.join(
            results_path,
            "flow_mapping",
            "origins_destinations_trade_economic_activity.csv",
        ),
        index=False,
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
