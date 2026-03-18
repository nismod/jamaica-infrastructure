"""
Intersect electricity network with hotspots grid for multi-point failure analysis.

Assigns grid IDs to power network nodes and edges. Uses the existing hotspots grid
so that electricity hotspots can be compared with other infrastructure hotspots.
"""

import logging

import click
import geopandas as gpd
import snail.intersection


logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--grid-path",
    "-g",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the hotspots grid TIFF file.",
)
@click.option(
    "--network-path",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the electricity network GPKG file.",
)
@click.option(
    "--output-nodes-path",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Path for output nodes geoparquet file with grid IDs.",
)
@click.option(
    "--output-edges-path",
    "-e",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Path for output edges geoparquet file with grid IDs.",
)
def main(grid_path, network_path, output_nodes_path, output_edges_path):
    """Intersect electricity network with hotspots grid."""
    logging.info("Reading grid and preparing metadata")
    grid = read_grid(grid_path)

    logging.info("Intersecting grid with nodes")
    intersect_grid_with_nodes(network_path, grid, output_nodes_path)

    logging.info("Intersecting grid with edges")
    intersect_grid_with_edges(network_path, grid, output_edges_path)

    logging.info("Complete")


def read_grid(grid_path):
    """Read grid and prepare metadata for snail intersection."""
    import rasterio
    with rasterio.open(grid_path) as src:
        grid = snail.intersection.GridDefinition(
            crs=src.crs,
            width=src.width,
            height=src.height,
            transform=src.transform
        )
    return grid


def intersect_grid_with_nodes(network_path, grid, output_path):
    """Assign grid coordinate indices to network nodes."""
    logging.info("Reading nodes from network")
    nodes = gpd.read_file(network_path, layer="nodes", engine="pyogrio")

    logging.info("Applying grid indices to nodes")
    grid_intersections = snail.intersection.apply_indices(
        nodes, grid, index_i="cell_index_0_x", index_j="cell_index_0_y"
    )

    logging.info(f"Writing {len(grid_intersections)} nodes to {output_path}")
    grid_intersections.to_parquet(output_path)


def intersect_grid_with_edges(network_path, grid, output_path):
    """Split network edges on grid boundaries and assign grid coordinate indices."""
    logging.info("Reading edges from network")
    edges = gpd.read_file(network_path, layer="edges", engine="pyogrio")

    logging.info("Preparing linestrings for splitting")
    edges = snail.intersection.prepare_linestrings(edges)

    logging.info("Splitting linestrings on grid boundaries")
    grid_intersections = snail.intersection.split_linestrings(edges, grid)

    logging.info("Calculating split edge lengths")
    grid_intersections["length_m"] = grid_intersections["geometry"].length

    logging.info("Applying grid indices to split edges")
    grid_intersections = snail.intersection.apply_indices(
        grid_intersections, grid, index_i="cell_index_0_x", index_j="cell_index_0_y"
    )

    logging.info(f"Writing {len(grid_intersections)} split edges to {output_path}")
    grid_intersections.to_parquet(output_path)


if __name__ == "__main__":
    main()
