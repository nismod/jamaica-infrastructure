#!/usr/bin/env python
# coding: utf-8
"""
Electricity network single-point edge failure analysis.

This script removes edges one at a time from the electricity network,
calculates the impact on service delivery, and outputs the results.

Usage:
    python electricity_edge_failure.py \\
        --nodes-file <path_to_nodes> \\
        --edges-file <path_to_edges> \\
        --flows-file <path_to_flows> \\
        --output-path <output_csv> \\
        --chunk-id <chunk_number> \\
        --chunk-count <total_chunks>
"""

import argparse
import logging
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from jem.model import jem
from jem.statistics import statistics

logging.basicConfig(
    format="%(asctime)s %(process)d %(levelname)s %(message)s", level=logging.WARN
)


def get_empty_edge_result(
    chunk_id: int,
    iteration_number: int,
    attacked_edge_id: str,
    attacked_edge_type: str,
    iteration_time_seconds: float,
) -> pd.DataFrame:
    """Create empty result dataframe when no nodes are affected."""
    return pd.DataFrame(
        {
            "chunk_id": [chunk_id],
            "iteration_number": [iteration_number],
            "attacked_edge_id": [attacked_edge_id],
            "affected_node_id": [np.nan],
            "attacked_edge_type": [attacked_edge_type],
            "affected_node_type": [np.nan],
            "total_nodes_affected": [np.nan],
            "population_affected": [np.nan],
            "demand_affected": [np.nan],
            "iteration_time_seconds": [iteration_time_seconds],
        }
    )


def analyse_edge_failure(
    edge_id: str,
    iteration_number: int,
    chunk_id: int,
    edges: gpd.GeoDataFrame,
    path_to_nodes: str,
    path_to_edges: str,
    path_to_flows: str,
) -> pd.DataFrame:
    """
    Run single-point failure analysis for a given edge.

    Args:
        edge_id: ID of the edge to attack
        iteration_number: Current iteration number
        edges: GeoDataFrame of all edges
        path_to_nodes: Path to nodes file
        path_to_edges: Path to edges file
        path_to_flows: Path to flows file

    Returns:
        DataFrame with failure analysis results
    """
    start_time = time.time()

    # Get attacked edge type
    attacked_edge_type = edges.loc[edges.id == edge_id, "asset_type"].iloc[0]

    # Run JEM model with edge failure
    try:
        run = jem(
            path_to_nodes,
            path_to_edges,
            path_to_flows,
            print_to_console=False,
            edges_to_attack=[edge_id],
            super_source=True,
            super_sink=False,
        )
        run.build()
        run.optimise(print_to_console=False)

        # Analyse results
        results = statistics(model_run=run)
        nodes_with_shortfall = results.nodes_with_shortfall().node.to_list()

        # If no nodes affected, return empty result
        if len(nodes_with_shortfall) == 0:
            time_taken = time.time() - start_time
            return get_empty_edge_result(
                chunk_id, iteration_number, edge_id, attacked_edge_type, time_taken
            )

        # Get population and demand impacts
        population_df = results.get_population_at_nodes(
            nodes_with_shortfall, col_id="affected_node_id"
        )
        demand_df = results.get_demand_at_nodes(
            nodes_with_shortfall, col_id="affected_node_id"
        )

        # Get affected node types
        node_types = results.nodes.loc[
            results.nodes.id.isin(nodes_with_shortfall), "asset_type"
        ].to_list()

        # Merge results
        df = population_df.merge(demand_df, on="affected_node_id")
        df["affected_node_type"] = node_types
        df["total_nodes_affected"] = len(nodes_with_shortfall)
        df["attacked_edge_id"] = edge_id
        df["attacked_edge_type"] = attacked_edge_type
        df["chunk_id"] = chunk_id
        df["iteration_number"] = iteration_number
        df["population_affected"] = df["population"]
        df["demand_affected"] = df["demand"]
        
        time_taken = time.time() - start_time
        df["iteration_time_seconds"] = time_taken
        
        # Reorder columns
        df = df[
            [
                "chunk_id",
                "iteration_number",
                "attacked_edge_id",
                "affected_node_id",
                "attacked_edge_type",
                "affected_node_type",
                "total_nodes_affected",
                "population_affected",
                "demand_affected",
                "iteration_time_seconds",
            ]
        ]

        return df

    except Exception as e:
        logging.error(f"Error analysing edge {edge_id}: {e}")
        time_taken = time.time() - start_time
        return get_empty_edge_result(
            chunk_id, iteration_number, edge_id, attacked_edge_type, time_taken
        )


def chunk_edges(
    edges: gpd.GeoDataFrame, chunk_id: int, chunk_count: int
) -> gpd.GeoDataFrame:
    """
    Split edges into chunks for parallel processing.

    Args:
        edges: GeoDataFrame of all edges
        chunk_id: Current chunk ID (0-indexed)
        chunk_count: Total number of chunks

    Returns:
        GeoDataFrame of edges in this chunk
    """
    if chunk_count == 1:
        return edges

    indices = np.round(np.linspace(0, len(edges), chunk_count + 1)).astype(int)
    start_idx = indices[chunk_id]
    end_idx = indices[chunk_id + 1]

    return edges.iloc[start_idx:end_idx].reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description="Electricity edge failure analysis")
    parser.add_argument(
        "--nodes-file",
        type=str,
        required=True,
        help="Path to electricity network nodes file (GeoPackage)",
    )
    parser.add_argument(
        "--edges-file",
        type=str,
        required=True,
        help="Path to electricity network edges file (GeoPackage)",
    )
    parser.add_argument(
        "--flows-file",
        type=str,
        required=True,
        help="Path to nodal flows CSV file",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        required=True,
        help="Path to output CSV file",
    )
    parser.add_argument(
        "--chunk-id",
        type=int,
        default=0,
        help="Chunk ID for parallel processing (0-indexed)",
    )
    parser.add_argument(
        "--chunk-count",
        type=int,
        default=1,
        help="Total number of chunks for parallel processing",
    )

    args = parser.parse_args()

    logging.info("Starting electricity edge failure analysis")
    logging.info(f"Chunk {args.chunk_id + 1} of {args.chunk_count}")

    # Read edges
    logging.info(f"Reading edges from {args.edges_file}")
    edges = gpd.read_file(args.edges_file, layer="edges", engine="pyogrio")

    # Chunk edges if needed
    if args.chunk_count > 1:
        edges = chunk_edges(edges, args.chunk_id, args.chunk_count)
        logging.info(f"Processing {len(edges)} edges in this chunk")

    # Get list of edges to attack
    edges_to_attack = edges.id.to_list()
    logging.info(f"Analysing {len(edges_to_attack)} edges")

    # Run analysis for each edge
    results_list = []
    for iteration, edge_id in enumerate(edges_to_attack, start=1):
        result_df = analyse_edge_failure(
            edge_id,
            iteration,
            args.chunk_id,
            edges,
            args.nodes_file,
            args.edges_file,
            args.flows_file,
        )
        results_list.append(result_df)

        if iteration % 100 == 0:
            logging.info(f"Completed {iteration}/{len(edges_to_attack)} edges")

    # Combine results
    logging.info("Combining results")
    results_df = pd.concat(results_list, ignore_index=True)

    # Save to CSV
    logging.info(f"Saving results to {args.output_path}")
    output_path = Path(args.output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, index=False)

    logging.info("Complete")


if __name__ == "__main__":
    main()
