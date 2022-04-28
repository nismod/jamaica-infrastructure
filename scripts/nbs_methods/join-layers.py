#!/usr/bin/env python
# coding: utf-8
import sys
import warnings

import geopandas
import numpy
import pandas
import rasterio
from tqdm import tqdm

# run from `nbs` data root /processed_data/nbs

def main(cell_number):
    # filter geopandas read_file warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    grid = geopandas.read_file('../grid_10km.gpkg').set_index('link_id')
    cell = grid.loc[cell_number]

    df = geopandas.read_file(f'land_use_slope/extracts/jamaica_land_use_combined_with_sectors_and_slope__{cell_number-1}.gpkg')

    df = df[['Classify', 'LU_CODE', 'cell_index', 'slope_degrees', 'geometry']] \
      .rename(columns={'Classify': 'landuse_desc', 'LU_CODE': 'landuse_code'})

    with_elevation = associate_raster('land_use_slope/elevation.tif', df, 'elevation_m')
    with_elevation.elevation_m = numpy.clip(with_elevation.elevation_m, -100, 2256)

    with_vectors = process(with_elevation, cell)
    with_vectors.to_file(f'land_use_slope/extracts/combined__{cell_number-1}.gpkg')


def associate_raster(fname, df, key, band_number=1, cell_index_col="cell_index"):
    with rasterio.open(fname) as dataset:
        assert df.crs == dataset.crs, "Raster and vector CRS must match"
        band_data = dataset.read(band_number)
        flat_data = band_data.flatten()
        df[key] = df[cell_index_col].apply(lambda i: flat_data[i])
    return df


def slice_vector(slice_df, slice_cell):
    slice_df = slice_df \
        .cx[slice_cell.left:slice_cell.right, slice_cell.bottom:slice_cell.top] \
        .copy()
    slice_df.geometry = slice_df.geometry.intersection(slice_cell.geometry)
    return slice_df

def associate_vector(fname, data_df, cell, select_cols):
    vector_df = geopandas.read_file(fname) \
        [list(select_cols.keys()) + ['geometry']] \
        .explode(index_parts=False) \
        .rename(columns=select_cols) \
        .copy()
    vector_df = slice_vector(vector_df, cell)
    if vector_df.empty:
        for col in select_cols.values():
            data_df[col] = None
        return data_df

    chunks = []
    chunk_size = 100
    for lower in tqdm(range(0, len(data_df), chunk_size)):
        chunk = data_df[lower:lower+chunk_size].copy()
        try:
            chunk = chunk.overlay(vector_df, how='identity', keep_geom_type=True)
        except ValueError:
            # probably had nothing to overlay
            for col in select_cols.values():
                chunk[col] = None
        chunks.append(chunk)

    # Alternate idea about clipping by subcell
    # for subcell in tqdm(subcells.itertuples(), total=len(subcells)):
    #     df_chunk = slice_vector(data_df, subcell)
    #     if df_chunk.empty:
    #         continue
    #     vector_chunk = slice_vector(vector_df, subcell)
    #     if vector_chunk.empty:
    #         chunks.append(df_chunk)
    #     else:
    #         chunk = df_chunk.overlay(vector_chunk, how='identity')
    #         chunks.append(chunk)

    return pandas.concat(chunks)

def process(with_elevation, cell):
    with_soils = associate_vector(
        'soils/nsmdb-soils.gpkg', with_elevation, cell,
        {'Type':'soil_type', 'Permeability Code': 'hydrologic_soil_group_code'})

    #with_erosion = associate_vector(
    #    'Terrestrial\Soil erosion susceptibility\Erosion susceptibility by land cover\Soil erosion susceptibility and land use intersect.shp',
    #    with_soils, cell
    #    {'Classes':'erosion_susceptibility'})

    with_bauxite = associate_vector(
        'Terrestrial/Bauxite/nsmdb-bauxite_reserves.gpkg', with_soils, cell,
        {'COLOR': 'within_bauxite_area'})
    with_bauxite.within_bauxite_area = ~with_bauxite.within_bauxite_area.isna()

    with_primary_forest = associate_vector(
        'Terrestrial/Forests/Forests buffered 100m/Primary forest buffered 100m.shp', with_bauxite, cell,
        {'LU_CODE': 'within_forest_100m'})
    with_primary_forest.within_forest_100m = (with_primary_forest.within_forest_100m == 'PF')

    with_forest_reserves = associate_vector(
        'Protected Sites/Forest Reserves/Forest Reserves.shp', with_primary_forest, cell,
        {'NAME':'forest_reserve_name'})
    with_forest_reserves['within_forest_reserve'] = ~with_forest_reserves.forest_reserve_name.isna()

    protected_areas = geopandas.read_file('Protected Sites/Protected Areas/Protected areas.shp')

    with_protected = with_forest_reserves
    for layer in protected_areas.LAYER.unique():
        with_protected = associate_vector(
            f'Protected Sites/Protected Areas/protected_areas_{layer}.gpkg',
            with_protected, cell,
            {'NAME':f'protected_area_{layer}_name'})
    with_protected['is_protected'] = ~(
        with_protected.protected_area_GAME_RESERVES_name.isna() |
        with_protected.protected_area_NATIONAL_PARK_name.isna() |
        with_protected.protected_area_PROTECTED_AREA_name.isna() |
        with_protected.protected_area_MARINE_PARK_name.isna()
    )
    with_protected['is_proposed_protected'] = ~with_protected.protected_area_PROPOSED_PROTECTED_AREA_name.isna()

    with_major = associate_vector(
        'Terrestrial/Riparian NbS/Rivers buffered 50m/Major rivers/Major rivers buffered 50m.gpkg',
        with_protected, cell,
        {'OBJECTID_1': 'within_major_river_50m'})
    with_major.within_major_river_50m = ~with_major.within_major_river_50m.isna()

    with_stream = associate_vector(
        'Terrestrial/Riparian NbS/Rivers buffered 50m/Large streams/Large steams buffered 50m.gpkg',
        with_major, cell,
        {'OBJECTID': 'within_large_stream_50m'})
    with_stream.within_large_stream_50m = ~with_stream.within_large_stream_50m.isna()

    with_headwater = associate_vector(
        'Terrestrial/Riparian NbS/Rivers buffered 50m/Headwater streams/Headwater streams buffered 50m.gpkg',
        with_stream, cell,
        {'OBJECTID': 'within_headwater_stream_50m'})
    with_headwater.within_headwater_stream_50m = ~with_headwater.within_headwater_stream_50m.isna()
    return with_headwater

if __name__ == '__main__':
    cell_number = sys.argv[1]
    main(int(cell_number))
