#!/usr/bin/env python
# coding: utf-8
import logging
import os
import warnings

import click
import fiona
import geopandas
import numpy
import pandas
import rasterio

from shapely.geometry import mapping, shape
from shapely.ops import linemerge, polygonize
from snail.core.intersections import get_cell_indices, split_linestring, split_polygon
from tqdm import tqdm
from jamaica_infrastructure.transform import read_transforms


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the asset layer table",
)
@click.option(
    "--hazard-csv",
    "-h",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the hazard table",
)
@click.option(
    "--data-dir",
    "-d",
    required=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=False, readable=True),
    help="Path to the data directory",
)
@click.option(
    "--asset-gpkg",
    "-g",
    required=True,
    help="asset_gpkg value in the network CSV",
)
@click.option(
    "--asset-layer",
    "-l",
    required=True,
    help="asset_layer value in the network CSV",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to save splits as geoparquet",
)
def split_assets(network_csv, hazard_csv, data_dir, asset_gpkg, asset_layer, output_path):

    # read networks
    networks = pandas.read_csv(network_csv)
    network_path = networks.loc[
        (networks.asset_gpkg == asset_gpkg) & (networks.asset_layer == asset_layer),
        "path"
    ].squeeze()

    hazards = pandas.read_csv(hazard_csv)
    hazard_transforms, transforms = read_transforms(hazards, data_dir)

    fname = os.path.join(data_dir, network_path)
    pq_fname_nodes = output_path.replace(".gpkg", "__nodes.geoparquet")
    pq_fname_edges = output_path.replace(".gpkg", "__edges.geoparquet")
    pq_fname_areas = output_path.replace(".gpkg", "__areas.geoparquet")

    logging.info("Processing %s", os.path.basename(fname))
    layers = fiona.listlayers(fname)
    logging.info("Layers: %s", layers)

    if "nodes" in layers:
        # look up nodes cell index
        nodes = geopandas.read_file(fname, layer="nodes")

        if nodes.empty:
            logging.info("Skipping %s with no data.", os.path.basename(fname))
        else:
            logging.info(os.path.basename(fname))
            logging.info("Node CRS %s", nodes.crs)
            nodes = process_nodes(nodes, transforms, hazard_transforms, data_dir)
            # nodes.to_file(out_fname, driver="GPKG", layer="nodes")
            nodes.to_parquet(pq_fname_nodes)

    if "edges" in layers:
        # split lines
        edges = geopandas.read_file(fname, layer="edges")

        if edges.empty:
            logging.info("Skipping %s with no data.", os.path.basename(fname))
        else:
            logging.info(os.path.basename(fname))
            logging.info("Node CRS %s", nodes.crs)
            logging.info("Edge CRS %s", edges.crs)
            edges = process_edges(edges, transforms, hazard_transforms, data_dir)
            # edges.to_file(out_fname, driver="GPKG", layer="edges")
            edges.to_parquet(pq_fname_edges)

    if "areas" in layers:
        # split polygons
        areas = geopandas.read_file(fname, layer="areas")

        if areas.empty:
            logging.info("Skipping %s with no data.", os.path.basename(fname))
        else:
            logging.info("Area CRS %s", areas.crs)
            areas = explode_multi(areas)
            areas = process_areas(areas, transforms, hazard_transforms, data_dir)
            # areas.to_file(out_fname, driver="GPKG", layer="areas")
            areas.to_parquet(pq_fname_areas)


def associate_raster(df, fname, cell_index_col="cell_index", band_number=1):
    with rasterio.open(fname) as dataset:
        band_data = dataset.read(band_number)
        return df[cell_index_col].apply(lambda i: band_data[i[1], i[0]])


def process_nodes(nodes, transforms, hazard_transforms, data_dir):
    # lookup per transform
    for i, t in enumerate(transforms):
        # transform to grid
        crs_df = nodes.to_crs(t.crs)
        # save cell index for fast lookup of raster values
        crs_df[f"cell_index_{i}"] = crs_df.geometry.progress_apply(
            lambda geom: get_indices(geom, t)
        )
        # transform back
        nodes = crs_df.to_crs(nodes.crs)

    # associate hazard values
    hazard_series: dict[str, pandas.Series] = {}
    for hazard in hazard_transforms.itertuples():
        logging.info("Hazard %s transform %s", hazard.key, hazard.transform_id)
        fname = os.path.join(data_dir, hazard.path)
        cell_index_col = f"cell_index_{hazard.transform_id}"
        hazard_series[hazard.key] = associate_raster(nodes, fname, cell_index_col)
    nodes = pandas.concat([nodes, pandas.DataFrame(hazard_series)], axis=1)

    # split and drop tuple columns so GPKG can save
    for i, t in enumerate(transforms):
        nodes = split_index_column(nodes, f"cell_index_{i}")
        nodes.drop(columns=f"cell_index_{i}", inplace=True)
    return nodes


def try_merge(geom):
    if geom.geom_type == "MultiLineString":
        geom = linemerge(geom)
    return geom


def process_edges(edges, transforms, hazard_transforms, data_dir):
    # handle multilinestrings
    edges.geometry = edges.geometry.apply(try_merge)
    geom_types = edges.geometry.apply(lambda g: g.geom_type)
    logging.info(geom_types.value_counts())
    edges = explode_multi(edges)

    # split edges per transform
    for i, t in enumerate(transforms):
        # transform to grid
        crs_df = edges.to_crs(t.crs)
        crs_df = split_df(crs_df, t)
        # save cell index for fast lookup of raster values
        crs_df[f"cell_index_{i}"] = crs_df.geometry.progress_apply(
            lambda geom: get_indices(geom, t)
        )
        # transform back
        edges = crs_df.to_crs(edges.crs)

    # associate hazard values
    hazard_series: dict[str, pandas.Series] = {}
    for hazard in hazard_transforms.itertuples():
        logging.info("Hazard %s transform %s", hazard.key, hazard.transform_id)
        fname = os.path.join(data_dir, hazard.path)
        cell_index_col = f"cell_index_{hazard.transform_id}"
        hazard_series[hazard.key] = associate_raster(edges, fname, cell_index_col)
    edges = pandas.concat([edges, pandas.DataFrame(hazard_series)], axis=1)

    # split and drop tuple columns so GPKG can save
    for i, t in enumerate(transforms):
        edges = split_index_column(edges, f"cell_index_{i}")
        edges.drop(columns=f"cell_index_{i}", inplace=True)

    return edges


def split_df(df, t):
    # split
    core_splits = []
    for edge in tqdm(df.itertuples(), total=len(df)):
        # split edge
        splits = split_linestring(edge.geometry, t.width, t.height, t.transform)
        # add to collection
        for s in splits:
            s_dict = edge._asdict()
            del s_dict["Index"]
            s_dict["geometry"] = s
            core_splits.append(s_dict)
    logging.info(f"Split {len(df)} edges into {len(core_splits)} pieces")
    sdf = geopandas.GeoDataFrame(core_splits, crs=t.crs, geometry="geometry")
    return sdf


def process_areas(areas, transforms, hazard_transforms, data_dir):
    # split areas per transform
    for i, t in enumerate(transforms):
        # transform to grid
        crs_df = areas.to_crs(t.crs)
        crs_df = split_area_df(crs_df, t)
        # save cell index for fast lookup of raster values
        crs_df[f"cell_index_{i}"] = crs_df.geometry.progress_apply(
            lambda geom: get_indices(geom, t)
        )
        # transform back
        areas = crs_df.to_crs(areas.crs)

    # associate hazard values
    hazard_series: dict[str, pandas.Series] = {}
    for hazard in hazard_transforms.itertuples():
        logging.info("Hazard %s transform %s", hazard.key, hazard.transform_id)
        fname = os.path.join(data_dir, hazard.path)
        cell_index_col = f"cell_index_{hazard.transform_id}"
        hazard_series[hazard.key] = associate_raster(areas, fname, cell_index_col)
    areas = pandas.concat([areas, pandas.DataFrame(hazard_series)], axis=1)

    # split and drop tuple columns so GPKG can save
    for i, t in enumerate(transforms):
        areas = split_index_column(areas, f"cell_index_{i}")
        areas.drop(columns=f"cell_index_{i}", inplace=True)

    return areas


def explode_multi(df):
    items = []
    geoms = []
    for item in df.itertuples(index=False):
        if item.geometry.geom_type in ("MultiPoint", "MultiLineString", "MultiPolygon"):
            for part in item.geometry.geoms:
                items.append(item._asdict())
                geoms.append(part)
        else:
            items.append(item._asdict())
            geoms.append(item.geometry)

    df = geopandas.GeoDataFrame(items, crs=df.crs, geometry=geoms)
    return df


def set_precision(geom, precision):
    """Set geometry precision"""
    geom_mapping = mapping(geom)
    geom_mapping["coordinates"] = numpy.round(
        numpy.array(geom_mapping["coordinates"]), precision
    )
    return shape(geom_mapping)


def split_area_df(df, t):
    # split
    core_splits = []
    for area in tqdm(df.itertuples(), total=len(df)):
        # split area
        try:
            splits = split_polygon(area.geometry, t.width, t.height, t.transform)
        except RuntimeError as e:
            print("SKIPPING GEOMETRY with error:", e)
            print(t)
            print(area.geometry.wkt)
            continue
        # round to high precision (avoid floating point errors)
        splits = [set_precision(s, 9) for s in splits]
        # to polygons
        splits = list(polygonize(splits))
        # add to collection
        for s in splits:
            s_dict = area._asdict()
            del s_dict["Index"]
            s_dict["geometry"] = s
            core_splits.append(s_dict)
    logging.info(f"  Split {len(df)} areas into {len(core_splits)} pieces")
    sdf = geopandas.GeoDataFrame(core_splits)
    sdf.crs = t.crs
    return sdf


def get_indices(geom, t):
    x, y = get_cell_indices(geom, t.width, t.height, t.transform)
    # wrap around to handle edge cases
    x = x % t.width
    y = y % t.height
    return (x, y)


def split_index_column(df, prefix):
    df[f"{prefix}_x"] = df[prefix].apply(lambda i: i[0])
    df[f"{prefix}_y"] = df[prefix].apply(lambda i: i[1])
    return df


if __name__ == "__main__":
    # Ignore writing-to-parquet warnings
    warnings.filterwarnings("ignore", message=".*initial implementation of Parquet.*")
    # Ignore reading-geopackage warnings
    warnings.filterwarnings(
        "ignore", message=".*Sequential read of iterator was interrupted.*"
    )

    # Enable progress_apply
    tqdm.pandas()

    # Enable info logging
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    logging.info("Start.")
    split_assets()
    logging.info("Done.")
