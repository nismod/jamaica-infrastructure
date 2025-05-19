"""
Rules for creating and analysing multi-modal transport networks.
"""


rule preprocess_road_network:
    """
    Preprocess road network data.
    Test with:
    snakemake -c1 processed_data/networks/transport/roads.gpkg
    """
    input:
        raw_roads = f"{RAW}/networks/transport/roads.gpkg",
    output:
        processed_roads = f"{DATA}/networks/transport/roads.gpkg",
    run:
        """
        Take an open-gira (or otherwise produced) road network and process it
        prior to inclusion in the multi-modal network.

        - Project to EPSG 3448
        - Explode multilinestrings
        - Calculate edge lengths
        - Assign asset_type for damage curve lookup
        - Estimate NWA road classification from OSM tag_highway
        - Gap-fill speed limit data
        - Add topology and component labels
        """

        import logging

        import geopandas as gpd
        import numpy as np
        import pandas as pd
        import snkit

        logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

        logging.info("Read raw road network data")
        edges = gpd.read_file(input.raw_roads, layer="edges")
        edges = edges.to_crs(epsg=3448)
        for col_name in ("from_node", "to_node", "edge_id"):
            if col_name in edges.columns:
                edges = edges.drop(columns=[col_name])
        logging.info(f"{len(edges)} edges")
        # cast multi-linestrings to linestrings
        exploded_geometry = edges.geometry.explode(index_parts=False)
        if len(edges == len(exploded_geometry)):
            edges.geometry = exploded_geometry
        else:
            raise ValueError("Multi-part geometries in raw road edges")

        logging.info("Remove edges with missing geometry")
        edges = edges[~edges.geometry.isna()]
        logging.info(f"{len(edges)} edges remain")

        logging.info("Calculate edge lengths in meters")
        edges["length_m"] = edges.geometry.length

        logging.info("Gap-filling speed limits with road classification modal value")
        class_to_speed_kph = edges.loc[:, ["tag_highway", "tag_maxspeed"]].groupby("tag_highway").agg(pd.Series.mode)
        edges["speed_kph"] = edges["tag_maxspeed"].copy()
        for road_class, modal_limit_kph in class_to_speed_kph.itertuples():
            if not 30 < modal_limit_kph < 130:
                raise ValueError(f"Suspicious speed limit: {modal_limit_kph} km/h for {road_class}")
            edges.loc[(edges.tag_highway == road_class) & edges.tag_maxspeed.isna(), "speed_kph"] = modal_limit_kph

        logging.info("Gap-filling # lanes")
        edges.lanes = edges.lanes.astype(int)
        edges[edges.lanes == 0] = 1

        logging.info("Set asset_type (for damage curve lookup)")
        edges["asset_type"] = edges.tag_highway.map(config["road_classification"]["OSM_to_damage_curve"])

        logging.info("Guess at NWA road classification from supplied mapping")
        edges["road_class"] = edges.tag_highway.map(config["road_classification"]["OSM_to_NWA"]) \
            .fillna(config["road_classification"]["default_NWA"])

        logging.info("Assign rebuild costs from config")
        USD_per_JMD = config["USD_per_JMD"]
        bridge_cost_JMD_per_meter = config["damages"]["rehabilitation_costs"]["bridge_cost_USD_per_meter"] / USD_per_JMD
        road_cost_JMD_per_lane_per_meter = (config["damages"]["rehabilitation_costs"]["road_cost_USD_per_lane_per_km"] / 1.0e3) / USD_per_JMD
        edges["mean_damage_cost"] = np.where(
            edges["asset_type"] == "road_bridge",
            bridge_cost_JMD_per_meter,
            road_cost_JMD_per_lane_per_meter * edges["lanes"],
        ) 
        edges["min_damage_cost"] = 0.8 * edges["mean_damage_cost"]
        edges["max_damage_cost"] = 1.2 * edges["mean_damage_cost"]
        edges["cost_unit"] = "J$/m"

        logging.info("Reading nodes")
        nodes = gpd.read_file(input.raw_roads, layer="nodes")
        nodes = nodes.to_crs(epsg=3448)
        if "node_id" in nodes.columns:
            nodes = nodes.drop(columns=["node_id"])
        logging.info(f"{len(nodes)} nodes")

        logging.info("Label road network with topology")
        network = snkit.Network(nodes=nodes, edges=edges)
        network = snkit.network.add_ids(network, edge_prefix="roade", node_prefix="roadn")
        network = snkit.network.add_topology(network, id_col="id")

        logging.info("Label bridge nodes (from_id of bridge edge)")
        # previously, bridge data from the NWA (2016) was used
        # that data has some of features present in OSM and some that are not
        # the NWA dataset contains some errors (false positives, imprecise locations)
        # and has fewer features in total
        # we decide in this rewrite to use the new OSM data exclusively

        # the downstream workflow is written assuming that bridges are nodes
        # we take the from_node of the OSM bridge edge to be the bridge node
        # this is not central to the span, but to one side of it
        bridge_edges = network.edges[network.edges.bridge]
        bridge_edges_to_merge = bridge_edges.loc[
            :,
            ["from_id", "length_m", "min_damage_cost", "mean_damage_cost", "max_damage_cost", "cost_unit"]
        ].rename(columns={"from_id": "id"})
        network.nodes = network.nodes.merge(bridge_edges_to_merge, how="outer", on="id")
        bridge_mask = ~network.nodes.mean_damage_cost.isna()
        network.nodes.loc[bridge_mask, "asset_type"] = "bridge"
        for cost_column in [f"{agg}_damage_cost" for agg in ["min", "mean", "max"]]:
            network.nodes[cost_column] *= network.nodes["length_m"]
        network.nodes.loc[bridge_mask, "cost_unit"] = "J$"

        logging.info("Label road network with components")
        network = snkit.network.add_component_ids(network)
        network.edges = network.edges.rename(
            columns={"from_id": "from_node", "to_id": "to_node", "id": "edge_id"},
        )
        network.nodes = network.nodes.rename(columns={"id": "node_id"})

        logging.info(f"Edges:\n{network.edges}")
        logging.info(f"Nodes:\n{network.nodes}")

        logging.info("Write processed road network to disk")
        network.edges.to_file(output.processed_roads, layer="edges", driver="GPKG")
        network.nodes.to_file(output.processed_roads, layer="nodes", driver="GPKG")


rule create_multi_modal_network:
    """
    Create a multi-modal transport network.

    Test with:
    snakemake -c1 processed_data/networks/transport/multi_modal_network.gpkg
    """
    input:
        script = "workflow/2a_transport/multi_modal_network_creation.py",
        airports = f"{DATA}/networks/transport/airport_polygon.gpkg",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        rail = f"{DATA}/networks/transport/rail.gpkg",
        road = f"{DATA}/networks/transport/roads.gpkg",
    output:
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
    shell:
        f"""
        python {{input.script}} \
            --road-network-path {{input.road}} \
            --rail-network-path {{input.rail}} \
            --port-path {{input.ports}} \
            --airport-path {{input.airports}} \
            --output-path {{output.multi_modal_network}}
        """


rule trade_activity_flow_mapping:
    """
    Create a mapping of trade activity to flows.

    Test with:
    snakemake -c1 results/flow_mapping/sector_imports_exports_to_ports_flows.gpkg
    """
    input:
        script = "workflow/2a_transport/trade_activity_flow_mapping.py",
        jam_ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        imports = f"{DATA}/macroeconomic_data/import_by_industry.xlsx",
        exports = f"{DATA}/macroeconomic_data/domestic_export_by_sector.xlsx",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        agriculture = f"{DATA}/agriculture_data/agriculture_gdp.gpkg",
        mining = f"{DATA}/mining_data/mining_gdp.gpkg",
    output:
        [
            "{output_path}/flow_mapping/sector_to_ports_flow_paths.csv",
            "{output_path}/flow_mapping/sector_to_ports_flow_paths.pq",
            "{output_path}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
            "{output_path}/flow_mapping/origins_destinations_trade_economic_activity.csv"
        ]
    shell:
        """
        python {input.script} \
            --ports {input.jam_ports} \
            --network {input.network} \
            --imports-xlsx {input.imports} \
            --exports-xlsx {input.exports} \
            --buildings-file {input.buildings} \
            --agriculture-file {input.agriculture} \
            --mining-file {input.mining} \
            --out-dir {wildcards.output_path}/flow_mapping
        """


rule labour_to_work_flow_mapping:
    """
    Create a mapping of labour to work flows.

    Test with:
    snakemake -c1 results/flow_mapping/labour_to_sectors_trips_and_activity.pq
    """
    input:
        script = "workflow/2a_transport/labour_to_work_flow_mapping.py",
        network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        population = f"{DATA}/population/population_projections.gpkg",
    output:
        [
            "{output_path}/flow_mapping/origins_destinations_labour_economic_activity.csv",
            "{output_path}/flow_mapping/road_nodes_labour_economic_activity_aggregations.gpkg",
            "{output_path}/flow_mapping/labour_to_sectors_flow_paths.pq",
            "{output_path}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
            "{output_path}/flow_mapping/labour_trips_and_activity.gpq",
        ]
    shell:
        """
        python {input.script} \
            --network-file {input.network} \
            --buildings-file {input.buildings} \
            --population-file {input.population} \
            --out-dir {wildcards.output_path}/flow_mapping
        """
