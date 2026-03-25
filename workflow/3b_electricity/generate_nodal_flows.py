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
    format="%(asctime)s %(process)d %(message)s", level=logging.INFO
)


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

    logging.info("Reshaping data to wide format")
    flow_nodes = flow_nodes[["id", "capacity"]].copy().rename(columns={"capacity": "flow"})
    list_of_nodes = flow_nodes.id.to_list()
    flow_nodes = flow_nodes.pivot_table(columns="id").reset_index(drop=True)
    flow_nodes["timestep"] = 1
    flow_nodes = flow_nodes[["timestep"] + list_of_nodes]
    flow_nodes = flow_nodes.round(2)

    logging.info(f"Saving flow data to {args.output_file}")
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    flow_nodes.to_csv(args.output_file, index=False)

    logging.info(f"Generated flow data for {len(list_of_nodes)} nodes")
    logging.info("Done")


if __name__ == "__main__":
    main()
