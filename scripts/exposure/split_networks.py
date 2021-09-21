#!/usr/bin/env python
# coding: utf-8
import json
import os
import sys
from collections import namedtuple
from glob import glob

import fiona
import geopandas
import pandas
import rasterio

from shapely.ops import linemerge
from snail.intersections import split_linestring
from snail.intersections import get_cell_indices
from tqdm import tqdm


def load_config():
    """Read config.json"""
    config_path = os.path.join(os.path.dirname(__file__), "..", "..", "config.json")
    with open(config_path, "r") as config_fh:
        config = json.load(config_fh)
    return config


# Enable progress_apply
tqdm.pandas()


def main(data_path, networks_csv, hazards_csv):
    # read transforms, record with hazards
    hazards = pandas.read_csv(hazards_csv)
    hazard_transforms, transforms = read_transforms(hazards, data_path)
    hazard_transforms.to_csv(hazards_csv.replace(".csv", "__with_transforms.csv"), index=False)

    # read networks
    networks = pandas.read_csv(networks_csv)

    for network_path in networks.path:
        fname = os.path.join(data_path, network_path)
        out_fname = fname.replace(".gpkg", "_splits.gpkg")

        # skip if output is there already
        if os.path.exists(out_fname):
            print("Skipping", os.path.basename(fname), ". Already exists:")
            print("-", out_fname)
            continue

        print("Processing", os.path.basename(fname))
        layers = fiona.listlayers(fname)
        print("   layers:", layers)

        if "nodes" in layers:
            # look up nodes cell index
            nodes = geopandas.read_file(fname, layer="nodes")
            print(" nodes", nodes.crs)
            nodes = process_nodes(nodes, transforms, hazard_transforms, data_path)
            nodes.to_file(out_fname, driver="GPKG", layer="nodes")

        if "edges" in layers:
            # split lines
            edges = geopandas.read_file(fname, layer="edges")
            print(" edges", edges.crs)
            edges = process_edges(edges, transforms, hazard_transforms, data_path)
            edges.to_file(out_fname, driver="GPKG", layer="edges")

        if "areas" in layers:
            # split polygons
            pass


# Helper class to store a raster transform and CRS
Transform = namedtuple('Transform', ['crs', 'width', 'height', 'transform'])


def associate_raster(df, key, fname, cell_index_col='cell_index', band_number=1):
    with rasterio.open(fname) as dataset:
        band_data = dataset.read(band_number)
        df[key] = df[cell_index_col].apply(lambda i: band_data[i[1], i[0]])


def read_transforms(hazards, data_path):
    transforms = []
    transform_id = 0
    hazard_transforms = []
    for hazard in hazards.itertuples():
        hazard_path = hazard.path
        with rasterio.open(os.path.join(data_path, hazard_path)) as dataset:
            crs = dataset.crs
            width = dataset.width
            height = dataset.height
            transform = Transform(crs, width, height, tuple(dataset.transform))
        # add transform to list if not present
        if transform not in transforms:
            transforms.append(transform)
            transform_id = transform_id + 1

        # record hazard/transform details
        hazard_transform_id = transforms.index(transform)
        hazard_transform = hazard._asdict()
        del hazard_transform['Index']
        hazard_transform['transform_id'] = hazard_transform_id
        hazard_transform['width'] = transform.width
        hazard_transform['height'] = transform.height
        hazard_transform['crs'] = str(transform.crs)
        hazard_transform['transform_0'] = transform.transform[0]
        hazard_transform['transform_1'] = transform.transform[1]
        hazard_transform['transform_2'] = transform.transform[2]
        hazard_transform['transform_3'] = transform.transform[3]
        hazard_transform['transform_4'] = transform.transform[4]
        hazard_transform['transform_5'] = transform.transform[5]
        hazard_transforms.append(hazard_transform)
    hazard_transforms = pandas.DataFrame(hazard_transforms)

    return hazard_transforms, transforms


def process_nodes(nodes, transforms, hazard_transforms, data_path):
    # lookup per transform
    for i, t in enumerate(transforms):
        # transform to grid
        crs_df = nodes.to_crs(t.crs)
        # save cell index for fast lookup of raster values
        crs_df[f'cell_index_{i}'] = crs_df.geometry.progress_apply(lambda geom: get_indices(geom, t))
        # transform back
        nodes = crs_df.to_crs(nodes.crs)

    # associate hazard values
    for hazard in hazard_transforms.itertuples():
        print(hazard.key, hazard.transform_id)
        fname = os.path.join(data_path, hazard.path)
        cell_index_col = f'cell_index_{hazard.transform_id}'
        associate_raster(nodes, hazard.key, fname, cell_index_col)

    # split and drop tuple columns so GPKG can save
    for i, t in enumerate(transforms):
        nodes = split_index_column(nodes, f'cell_index_{i}')
        nodes.drop(columns=f'cell_index_{i}', inplace=True)
    return nodes


def try_merge(geom):
    if geom.type =='MultiLineString':
        geom = linemerge(geom)
    return geom


def process_edges(edges, transforms, hazard_transforms, data_path):
    # handle multilinestrings
    edges.geometry = edges.geometry.apply(try_merge)
    geom_types = edges.geometry.apply(lambda g: g.type)
    print("   ", (geom_types.value_counts()))
    edges = edges[geom_types == 'LineString']

    # split edges per transform
    for i, t in enumerate(transforms):
        # transform to grid
        crs_df = edges.to_crs(t.crs)
        crs_df = split_df(crs_df, t)
        # save cell index for fast lookup of raster values
        crs_df[f'cell_index_{i}'] = crs_df.geometry.progress_apply(lambda geom: get_indices(geom, t))
        # transform back
        edges = crs_df.to_crs(edges.crs)

    # associate hazard values
    for hazard in hazard_transforms.itertuples():
        print(hazard.key, hazard.transform_id)
        fname = os.path.join(data_path, hazard.path)
        cell_index_col = f'cell_index_{hazard.transform_id}'
        associate_raster(edges, hazard.key, fname, cell_index_col)

    # split and drop tuple columns so GPKG can save
    for i, t in enumerate(transforms):
        edges = split_index_column(edges, f'cell_index_{i}')
        edges.drop(columns=f'cell_index_{i}', inplace=True)

    return edges


def split_df(df, t):
    # split
    core_splits = []
    for edge in tqdm(df.itertuples(), total=len(df)):
        # split edge
        splits = split_linestring(
            edge.geometry,
            t.width,
            t.height,
            t.transform
        )
        # add to collection
        for s in splits:
            s_dict = edge._asdict()
            del s_dict['Index']
            s_dict['geometry'] = s
            core_splits.append(s_dict)
    print("  Split %d edges into %d pieces", len(df), len(core_splits))
    sdf = geopandas.GeoDataFrame(core_splits)
    sdf.crs = t.crs
    return sdf


def get_indices(geom, t):
    x, y = get_cell_indices(
        geom,
        t.width,
        t.height,
        t.transform)
    # wrap around to handle edge cases
    x = x % t.width
    y = y % t.height
    return (x, y)


def split_index_column(df, prefix):
    df[f'{prefix}_x'] = df[prefix].apply(lambda i: i[0])
    df[f'{prefix}_y'] = df[prefix].apply(lambda i: i[1])
    return df


if __name__ == '__main__':
    # Load config
    CONFIG = load_config()
    # Save splits, attribute data later
    data_path = CONFIG["paths"]["data"]
    try:
        networks_csv = sys.argv[1]
        hazards_csv = sys.argv[2]
    except IndexError:
        print("Error. Please provide networks and hazards as CSV.")
        print(f"Usage: python {__file__} networks/network_files.csv hazards/hazard_layers.csv")
    main(data_path, networks_csv, hazards_csv)
