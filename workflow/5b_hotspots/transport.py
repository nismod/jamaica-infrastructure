"""
Remove all edges within a hotspot grid cell and calculate the resulting economic
losses due to trade and labour disruption.
"""

import logging

import click
import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from jamaica_infrastructure.transport.econ import economic_losses_from_network_damage
from jamaica_infrastructure.transport.flow import (
    igraph_scenario_edge_failures_premade_network,
    read_flow_data,
)


tqdm.pandas()


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--labour-cost-JMD-per-hour", "-l", "labour_cost_JMD_per_hour", required=True, type=float,
    help="Average cost of labour in JMD per hour."
)
@click.option(
    "--trade-rerouting-fraction", "-t", required=True, type=float,
    help="Fraction of trade value that will be affected by rerouting."
)
@click.option(
    "--grid-y-index", "-y", required=True, type=int,
    help="Index of hotspots grid row (constant latitude) to process."
)
@click.option(
    "--road-splits-path", "-r", required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to road edges split on hotspots grid",
)
@click.option(
    "--rail-splits-path", "-a", required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to rail edges split on hotspots grid",
)
@click.option(
    "--flow-data-dir", "-f", required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Path to the directory containing nominal flow data (ODs, etc.).",
)
@click.option(
    "--output-path", "-o", required=True,
    type=click.Path(exists=False, file_okay=False, dir_okay=True, readable=True),
    help="Directory to write latitude slice results to.",
)
def main(
    labour_cost_JMD_per_hour,
    trade_rerouting_fraction,
    grid_y_index,
    road_splits_path,
    rail_splits_path,
    flow_data_dir,
    output_path
) -> None:
    """
    Calculate economic losses due to removing all road and rail edges within a
    hotspots grid cell from multi-modal transport network.
    """

    logging.info(f"Running transport failure analysis for hotspot cells with {grid_y_index=}")

    logging.info("Read nominal flow data and networks")
    network_data, all_flows, trade_sectors = read_flow_data(flow_data_dir)

    logging.info("Reading splits data")
    splits = pd.concat([gpd.read_parquet(road_splits_path), gpd.read_parquet(rail_splits_path)])

    logging.info("Failing edges within cells, reallocating flows and calculating losses")
    losses_by_cell: list[pd.DataFrame] = []
    grid_x_indices: list[int] = sorted(splits.loc[splits.cell_index_0_y == grid_y_index, "cell_index_0_x"].unique())
    for grid_x_index in grid_x_indices:
        logging.info(f"({grid_y_index}, {int(grid_x_index)})")

        edges_to_remove: list[str] = sorted(
            splits.loc[
                (splits.cell_index_0_y == grid_y_index) & (splits.cell_index_0_x == grid_x_index),
                "edge_id"
            ].unique()
        )

        rerouting_cost: list[dict] = []
        for flow_type, data in network_data.items():
            rerouting_cost += igraph_scenario_edge_failures_premade_network(
                # we will remove edges from the graph, only operate on a copy
                data["network"].copy(),
                edges_to_remove,
                data["flows"],
                data["edge_indexes"],
                "edge_path",
                "time",
            )

        if not rerouting_cost:
            continue
        else:
            losses: pd.DataFrame = economic_losses_from_network_damage(
                pd.DataFrame(rerouting_cost),
                "edge_id",
                "no_access",
                all_flows,
                trade_sectors,
                trade_rerouting_fraction,
                labour_cost_JMD_per_hour,
            )

        losses["cell_index_y"] = grid_y_index
        losses["cell_index_x"] = grid_x_index

        losses_by_cell.append(losses)

    logging.info("Writing results to disk")
    if losses_by_cell:
        output: pd.DataFrame = pd.concat(losses_by_cell)
    else:
        # empty DataFrame
        output: pd.DataFrame = pd.DataFrame(losses_by_cell)
    output.to_parquet(output_path)


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    main()
