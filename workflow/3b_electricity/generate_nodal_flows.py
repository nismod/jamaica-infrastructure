#!/usr/bin/env python
# coding: utf-8
"""
Generate nodal flow data for electricity network.

This script reads the electricity network nodes and creates a file with
supply (from sources) and demand (at sinks).

Usage:
    python generate_nodal_flows.py \\
        --network-file <path_to_network.gpkg> \\
        --output-file <path_to_output.csv>
"""

import argparse
import logging
from pathlib import Path

import geopandas as gpd


logging.basicConfig(
    format="%(asctime)s %(process)d %(levelname)s %(message)s", level=logging.INFO
)


def compute_supply_demand(nodes: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Calculate supply and demand flows in consistent units (kW).

    Args:
        nodes: GeoDataFrame containing network nodes

    Returns:
        GeoDataFrame with 'flow' column added
    """
    nodes["flow"] = 0

    nodes.loc[nodes.asset_type == "source", "flow"] = (
        nodes.loc[nodes.asset_type == "source", "capacity"] * 10**3
    )

    nodes.loc[nodes.asset_type == "sink", "flow"] = (
        nodes.loc[nodes.asset_type == "sink", "population"]
        * nodes.loc[nodes.asset_type == "sink", "ei"]
        * 10**3
    )

    return nodes


def main():
    parser = argparse.ArgumentParser(
        description="Generate nodal flow data for electricity network"
    )
    parser.add_argument(
        "--network-file",
        type=str,
        required=True,
        help="Path to network GeoPackage file",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path to output CSV file",
    )

    args = parser.parse_args()

    logging.info(f"Reading network from {args.network_file}")
    nodes = gpd.read_file(args.network_file, layer="nodes")

    logging.info("Filtering to source and sink nodes")
    flow_nodes = nodes[nodes.asset_type.isin(["source", "sink"])].copy()
    flow_nodes = flow_nodes.reset_index(drop=True)

    logging.info("Computing supply and demand flows")
    flow_nodes = compute_supply_demand(flow_nodes)

    logging.info("Reshaping data to wide format")
    flow_nodes = flow_nodes[["id", "flow"]]
    list_of_nodes = flow_nodes.id.to_list()
    flow_nodes = flow_nodes.pivot_table(columns="id").reset_index(drop=True)
    flow_nodes["timestep"] = 1
    flow_nodes = flow_nodes[["timestep"] + list_of_nodes]
    flow_nodes = flow_nodes.round(1)

    logging.info(f"Saving flow data to {args.output_file}")
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    flow_nodes.to_csv(args.output_file, index=False)

    logging.info(f"Generated flow data for {len(list_of_nodes)} nodes")
    logging.info("Done")


if __name__ == "__main__":
    main()
