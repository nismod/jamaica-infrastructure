#!/usr/bin/env python
# coding: utf-8
"""
Script to run nominal electricity network flows and save results.

Runs the JEM model with no asset failures and saves the resulting power flows
and network topology to a geoparquet file for visualization.

Usage:
    python electricity_nominal_flow_diagnostic.py \\
        --network-path processed_data/networks/energy/electricity_network_v3.2.gpkg \\
        --flows-path processed_data/networks/energy/generated_nodal_flows.csv \\
        --output-path results/electricity_failures/diagnostics/nominal_flows.geoparquet
"""

import logging

import click
import geopandas as gpd
import numpy as np
from pathlib import Path

from jem.model import jem
from jem.statistics import statistics

logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)


def log_nonzero(x):
    """Apply log transform to non-zero values."""
    with np.errstate(divide="ignore"):
        return np.log10(x, out=np.zeros_like(x), where=0 < x)


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--network-path",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the electricity network GPKG file.",
)
@click.option(
    "--flows-path",
    "-f",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the nodal flows CSV file.",
)
@click.option(
    "--output-edges-path",
    "-oe",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Path for output edges geoparquet file.",
)
@click.option(
    "--output-nodes-path",
    "-on",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Path for output nodes geoparquet file.",
)
def main(network_path, flows_path, output_edges_path, output_nodes_path):
    """Run JEM model and save nominal network flows to geoparquet."""

    logging.info("Reading network topology")
    nodes = gpd.read_file(network_path, layer="nodes", engine="pyogrio")
    edges = gpd.read_file(network_path, layer="edges", engine="pyogrio")

    logging.info(f"Network: {len(nodes)} nodes, {len(edges)} edges")

    logging.info("Building and optimizing JEM model (nominal, no failures)")
    run = jem(
        str(network_path),
        str(network_path),
        str(flows_path),
        print_to_console=False,
        nodes_to_attack=[],
        edges_to_attack=[],
        super_source=True,
        super_sink=False,
    )
    run.build()
    run.optimise(print_to_console=False)

    results = statistics(model_run=run)
    edge_flows = results.edge_flows
    nodes = nodes.merge(results.nodes_with_shortfall().loc[:, ["node", "shortfall"]].rename(columns={"node": "id"}), on="id", how="outer")

    edges_idx = edges.set_index(["from_id", "to_id"])
    flows_idx = edge_flows.set_index(["from_id", "to_id"])
    edges_with_flows = edges_idx.join(flows_idx, how="left")

    edges_with_flows["flow"] = edges_with_flows["flow"].fillna(0)
    edges_with_flows = edges_with_flows.reset_index()
    edges_with_flows["log10_flow"] = log_nonzero(edges_with_flows["flow"])
    max_flow = edges_with_flows["flow"].max()

    edges_with_flow = (edges_with_flows["flow"] > 0).sum()
    shortfall_W = nodes.shortfall.sum() * 1E3
    total_demand_W = nodes[nodes.asset_type == "sink"].capacity.sum() * 1E3
    total_available_supply_W = nodes[nodes.asset_type == "source"].capacity.sum() * 1E3

    logging.info("Network statistics:")
    logging.info(f"  Total edges: {len(edges_with_flows)}")
    logging.info(f"  Edges with flow: {edges_with_flow}")
    logging.info(f"  Max edge flow: {max_flow * 1E3:,.0f} W")
    logging.info(f"  Total available supply: {total_available_supply_W:,.0f} W")
    logging.info(f"  Total demand: {total_demand_W:,.0f} W")
    logging.info(f"  Supply shortfall: {shortfall_W:,.0f} W ({shortfall_W / total_demand_W * 100:,.2f}% of total demand)")

    logging.info(f"Saving edges with flows to {output_edges_path}")
    Path(output_edges_path).parent.mkdir(parents=True, exist_ok=True)
    edges_with_flows.to_parquet(output_edges_path, compression="gzip")
    edges_with_flows.to_file(
        str(output_edges_path).replace(".geoparquet", ".gpkg"), layer="edges"
    )
    nodes.to_parquet(output_nodes_path, compression="gzip")

    logging.info("Complete!")
    logging.info(f"  Edge flows saved to: {output_edges_path}")
    logging.info(f"  Node shortfalls saved to: {output_nodes_path}")


if __name__ == "__main__":
    main()
