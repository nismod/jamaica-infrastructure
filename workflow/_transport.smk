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
        - Gap-fill speed limit data
        - Add topology and component labels

        TODO: Add government road categorisation data (could spatial join
        against NWA as Raghav did in scripts/preprocess/updated_road_network.py)
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

        logging.info("Assign rebuild costs")
        USD_per_JMD = 0.0064  # TODO: factor out cost assumptions
        # From NWA - Cost of two lane road reconstruction is US$ 1.5 million/km
        road_cost_USD_per_lane_per_km = 0.75E6
        bridge_cost_USD_per_meter = 1.5E6

        road_cost_JMD_per_lane_per_meter = (road_cost_USD_per_lane_per_km / 1.0e3) / USD_per_JMD
        bridge_cost_JMD_per_meter = bridge_cost_USD_per_meter / USD_per_JMD
        edges["cost_unit"] = "J$/m"
        edges["mean_damage_cost"] = np.where(
            edges["asset_type"] == "road_bridge",
            bridge_cost_JMD_per_meter * edges["length_m"],
            road_cost_JMD_per_lane_per_meter * edges["length_m"] * edges["lanes"],
        )
        edges["min_damage_cost"] = 0.8 * edges["mean_damage_cost"]
        edges["max_damage_cost"] = 1.2 * edges["mean_damage_cost"]

        nodes = gpd.read_file(input.raw_roads, layer="nodes")
        nodes = nodes.to_crs(epsg=3448)
        if "node_id" in nodes.columns:
            nodes = nodes.drop(columns=["node_id"])
        logging.info(f"{len(nodes)} nodes")

        logging.info("Label road network with topology")
        network = snkit.Network(nodes=nodes, edges=edges)
        network = snkit.network.add_ids(network, edge_prefix="roade", node_prefix="roadn")
        network = snkit.network.add_topology(network, id_col="id")
        logging.info("Label road network with components")
        network = snkit.network.add_component_ids(network)
        network.edges = network.edges.rename(
            columns={"from_id": "from_node", "to_id": "to_node", "id": "edge_id"},
        )
        network.nodes = network.nodes.rename(columns={"id": "node_id"})

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
        airports = f"{DATA}/networks/transport/airport_polygon.gpkg",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        rail = f"{DATA}/networks/transport/rail.gpkg",
        road = f"{DATA}/networks/transport/roads.gpkg",
    output:
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
    script:
        "../scripts/transport_model/multi_modal_network_creation.py"