#!/usr/bin/env python
# coding: utf-8
"""Split land use areas into slope squares
"""
import logging
import sys
from glob import glob

import geopandas
import numpy
import rasterio
import shapely.geometry
import shapely.ops
import snail.core.intersections
from tqdm import tqdm

tqdm.pandas()


def associate_raster(fname, df, key, band_number=1, cell_index_col="cell_index"):
    with rasterio.open(fname) as dataset:
        assert df.crs == dataset.crs, "Raster and vector CRS must match"
        band_data = dataset.read(band_number)
        flat_data = band_data.flatten()
        df[key] = df[cell_index_col].apply(lambda i: flat_data[i])
    return df


def set_precision(geom, precision):
    """Set geometry precision"""
    geom_mapping = shapely.geometry.mapping(geom)
    geom_mapping["coordinates"] = numpy.round(
        numpy.array(geom_mapping["coordinates"]), precision
    )
    return shapely.geometry.shape(geom_mapping)


def split_area_df(df, width, height, transform):
    # split
    core_splits = []
    for area in tqdm(df.itertuples(), total=len(df)):
        # split area
        splits = snail.core.intersections.split_polygon(
            area.geometry, width, height, transform
        )
        # round to high precision (avoid floating point errors)
        # 7-9 is appropriate for degrees, 2-3 more appropriate for metres
        splits = [set_precision(s, 3) for s in splits]
        # to polygons
        splits = list(shapely.ops.polygonize(splits))
        # add to collection
        for s in splits:
            s_dict = area._asdict()
            del s_dict["Index"]
            s_dict["geometry"] = s
            core_splits.append(s_dict)
    logging.info(f"  Split {len(df)} areas into {len(core_splits)} pieces")
    sdf = geopandas.GeoDataFrame(core_splits)
    sdf.crs = df.crs
    return sdf


def get_index(geom, width, height, transform):
    x, y = snail.core.intersections.get_cell_indices(
        geom, width, height, transform
    )
    # wrap around to handle edge cases
    x = x % width
    y = y % height
    return x + width * y


def main(slope_fname, landuse_fname):
    with rasterio.open(slope_fname) as slope:
        pass

    # read polygons
    land_use = (
        geopandas.read_file(
            landuse_fname
        )[["Classify", "LU_CODE", "geometry"]]
        .explode()
        .reset_index(drop=True)
    )

    # subset for testing
    #land_use = land_use[:100]

    # split
    land_use_split = split_area_df(land_use, slope.width, slope.height, slope.transform)

    # save cell index for fast lookup of raster values
    land_use_split["cell_index"] = land_use_split.geometry.progress_apply(
        lambda geom: get_index(geom, slope.width, slope.height, slope.transform)
    )

    # associate hazard values
    land_use_slope = associate_raster(slope_fname, land_use_split, "slope_degrees")

    land_use_slope.to_file("processed_data/nbs/landuse_slope.gpkg", driver="GPKG")

if __name__ == '__main__':
    try:
        slope_fname, landuse_fname = sys.argv[1:]
    except:
        print(sys.argv[1:])
        print("Expected usage:")
        print(f"    python {__file__} slope.tif landuse.gpkg")
        exit()

    main(slope_fname, landuse_fname)

