"""
Disrupt nodes and edges within a hotspot grid cell.

This script runs multi-point failure analysis for a given grid cell,
removing all power infrastructure within that cell and calculating
the impact on population and demand.
"""

import logging
import time
from pathlib import Path

import click
import geopandas as gpd
import numpy as np
import pandas as pd
from jem.statistics import statistics
from jem.model import jem


def get_empty_results(cell_x, cell_y):
    """Create empty result when no nodes are affected."""
    return pd.DataFrame(
        {
            "cell_index_x": [cell_x],
            "cell_index_y": [cell_y],
            "affected_node_id": [np.nan],
            "population_affected": [np.nan],
            "demand_affected": [np.nan],
        }
    )


def process_grid_cell(
    cell_x,
    cell_y,
    network_path,
    flows_path,
    nodes_with_grid_path,
    edges_with_grid_path,
):
    """Run multi-point failure analysis for a grid cell."""
    logging.info(f"Processing grid cell ({cell_x}, {cell_y})")

    # Read the split network (with linked grid ids) to select failures
    logging.info("Reading network splits with grid IDs")
    nodes_split = gpd.read_parquet(nodes_with_grid_path)
    # Note: in the current model nodes of type 'sink' cannot be failed
    nodes_split = nodes_split[nodes_split["asset_type"] != "sink"].copy()
    edges_split = gpd.read_parquet(edges_with_grid_path)

    # Get nodes and edges to attack in this grid cell
    nodes_to_attack = nodes_split[
        (nodes_split.cell_index_0_x == cell_x) & (nodes_split.cell_index_0_y == cell_y)
    ].id.to_list()
    edges_to_attack = (
        edges_split[
            (edges_split.cell_index_0_x == cell_x)
            & (edges_split.cell_index_0_y == cell_y)
        ]
        .id.unique()
        .tolist()
    )

    if not (bool(nodes_to_attack) or bool(edges_to_attack)):
        logging.info(f"Grid cell ({cell_x}, {cell_y}): No edges or nodes present")
        return get_empty_results(cell_x, cell_y)

    logging.info(
        f"Grid cell ({cell_x}, {cell_y}): attacking {len(nodes_to_attack)} nodes "
        f"and {len(edges_to_attack)} edges"
    )
    # Run model
    logging.info("Building and optimizing JEM model")
    run = jem(
        str(network_path),
        str(network_path),
        str(flows_path),
        print_to_console=False,
        nodes_to_attack=nodes_to_attack,
        edges_to_attack=edges_to_attack,
        super_source=True,
        super_sink=False,
    )
    run.build()
    run.optimise(print_to_console=False)

    logging.info("Analyzing results")
    results = statistics(model_run=run)
    nodes_with_shortfall = results.nodes_with_shortfall().node.to_list()

    if results.nodes_with_shortfall().empty:
        logging.info(f"Grid cell ({cell_x}, {cell_y}): No nodes affected")
        return get_empty_results(cell_x, cell_y)

    logging.info(
        f"Grid cell ({cell_x}, {cell_y}): {len(nodes_with_shortfall)} nodes affected"
    )
    population = results.get_population_at_nodes(
        nodes_with_shortfall, col_id="affected_node_id"
    )
    demand = results.get_demand_at_nodes(
        nodes_with_shortfall, col_id="affected_node_id"
    )
    df = population.merge(demand, on="affected_node_id")
    df["cell_index_x"] = cell_x
    df["cell_index_y"] = cell_y
    df["population_affected"] = df["population"]
    df["demand_affected"] = df["demand"]
    df = df[
        [
            "cell_index_x",
            "cell_index_y",
            "affected_node_id",
            "population_affected",
            "demand_affected",
        ]
    ]
    return df


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--cell-x",
    "-x",
    required=True,
    type=int,
    help="Grid cell X coordinate to analyze.",
)
@click.option(
    "--cell-y",
    "-y",
    required=True,
    type=int,
    help="Grid cell Y coordinate to analyze.",
)
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
    "--nodes-with-grid-path",
    "-p",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to nodes with grid IDs.",
)
@click.option(
    "--edges-with-grid-path",
    "-e",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to edges with grid IDs.",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Path for output CSV file.",
)
def main(
    cell_x,
    cell_y,
    network_path,
    flows_path,
    nodes_with_grid_path,
    edges_with_grid_path,
    output_path,
):
    """Run multi-point failure analysis for a grid cell."""
    start_time = time.time()

    df = process_grid_cell(
        cell_x,
        cell_y,
        network_path,
        flows_path,
        nodes_with_grid_path,
        edges_with_grid_path,
    )
    logging.info(f"Writing results to {output_path}")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    elapsed = time.time() - start_time
    logging.info(f"Grid cell ({cell_x}, {cell_y}) complete in {elapsed:.1f} seconds")


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    main()
