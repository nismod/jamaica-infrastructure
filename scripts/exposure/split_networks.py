#!/usr/bin/env python
# coding: utf-8
import os
from collections import namedtuple
from glob import glob

import fiona
import geopandas
import pandas
import rasterio
import snail

from shapely.ops import linemerge
from snail.intersections import split
from snail.intersections import get_cell_indices
from tqdm import tqdm

# Enable progress_apply
tqdm.pandas()

def main(data_path):
    hazards = pandas.read_csv(os.path.join(data_path, "hazards", "hazard_layers.csv"))

    networks = pandas.read_csv(os.path.join(data_path, "networks", "network_files.csv"))
    for network_path in networks.path:
        fname = os.path.join(data_path, network_path)
        out_fname = fname.replace(".gpkg", "_splits.gpkg")
        print(os.path.basename(fname))
        layers = fiona.listlayers(fname)
        print("   ", layers)
        if "nodes" in layers:
            # look up nodes values
            pass
        if "edges" in layers:
            # split lines
            edges = geopandas.read_file(fname, layer="edges")
            print(" edges", edges.crs)
            edges = process_edges(edges, transforms)
            edges.to_file(out_fname, driver="GPKG", layer="edges")
        if "areas" in layers:
            # split polygons
            pass    


# Helper class to store a raster transform and CRS
Transform = namedtuple('Transform', ['crs', 'width', 'height', 'transform'])


def associate_raster(df, key, fname, band_number=1):
    with rasterio.open(fname) as dataset:
        band_data = dataset.read(band_number)
        df[key] = df.cell_index.apply(lambda i: band_data[i[1], i[0]])


def read_transforms(hazards):
    transforms = set()
    for hazard_path in hazards.path:
        with rasterio.open(os.path.join(data_path, hazard_path)) as dataset:
            crs = dataset.crs
            width = dataset.width
            height = dataset.height
            transforms.add(Transform(crs, width, height, tuple(dataset.transform)))
    transforms = list(transforms)
    return transforms


def process_edges(edges, transforms):
    # handle multilinestrings
    edges.geometry = edges.geometry.apply(linemerge)
    print("   ", (edges.geometry.apply(lambda g: g.type)).value_counts())
    edges = edges[edges.geometry.apply(lambda g: g.type) == 'LineString']

    # split edges per transform
    for i, t in enumerate(transforms):
        print(i, t)
        edges = split_df(edges, t)
        # save cell index for fast lookup of raster values
        edges[f'cell_index_{i}'] = edges.geometry.progress_apply(lambda geom: get_indices(geom, t))
    return edges


def split_df(df, t):
    # transform to grid
    crs_df = df.to_crs(t.crs)
    # split
    core_splits = []
    for edge in tqdm(crs_df.itertuples(), total=len(df)):
        # split edge
        splits = split(
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
    # transform back
    return sdf.to_crs(df.crs)


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

if __name__ == '__main__':
    # Load config

    # Save splits, attribute data later
    pass
