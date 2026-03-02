"""Assign mining GDP to land use layers in Jamaica"""

import sys
import os
import subprocess

import pandas as pd
import geopandas as gpd
import click

# gpd._compat.USE_PYGEOS = True
# gpd.options.use_pygeos = True
from shapely.geometry import Point
import shapely
import numpy as np
from utils import *
from tqdm import tqdm

tqdm.pandas()


@click.command()
@click.option(
    "--epsg",
    "-e",
    "epsg",
    required=True,
    type=int,
    help="Coordinate system for Jamaica",
)
@click.option(
    "--mining-gdp",
    "-mg",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to mining GDP gpkg input",
)
@click.option(
    "--intermediate-file",
    "-if",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to buildings assigned economic sectors intermediate gpkg",
)
@click.option(
    "--mining-gdp-output",
    "-mgo",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    help="Path to output mining GDP with buildings gpkg",
)
@click.option(
    "--building-mining-gdp",
    "-bmg",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    help="Path to output building mining GDP csv",
)
def main(
    epsg,
    mining_gdp,
    intermediate_file,
    mining_gdp_output,
    building_mining_gdp,
):
    epsg_jamaica = epsg

    mining_areas = gpd.read_file(
        mining_gdp,
        layer="areas",
    )

    tot_gpd = mining_areas["C_GDP"].sum()
    tot_area = mining_areas["area_m2"].sum()
    print("* Estimated GDP", tot_gpd)
    print("* Estimated Areas", tot_area)
    sector_df = gpd.read_file(
        intermediate_file,
        layer="commercial_sectors",
    )
    sector_df = sector_df.to_crs(epsg=epsg_jamaica)
    sector_columns = [c for c in sector_df.columns.values.tolist() if "C_" in c[:2]]
    sector_df["C"] = sector_df[sector_columns].sum(axis=1)
    sector_df["C"] = sector_df.progress_apply(lambda x: 1 if x["C"] > 0 else 0, axis=1)
    mining_buildings = gpd.sjoin(
        mining_areas[["mining_id", "GDP_persqm", "geometry"]],
        sector_df[sector_df["C"] == 1][["osm_id", "geometry"]],
        how="inner",
        predicate="intersects",
    ).reset_index()

    mining_buildings.rename(columns={"geometry": "landuse_geometry"}, inplace=True)
    mining_buildings = pd.merge(
        mining_buildings, sector_df[["osm_id", "geometry"]], how="left", on=["osm_id"]
    )
    mining_buildings["area_sqm"] = mining_buildings.progress_apply(
        lambda x: x["landuse_geometry"].intersection(x["geometry"].buffer(0)).area,
        axis=1,
    )
    print("* Building areas", mining_buildings["area_sqm"].sum())
    print("* Building areas ratio", mining_buildings["area_sqm"].sum() / tot_area)
    mining_buildings["C_GDP_building"] = (
        mining_buildings["GDP_persqm"] * mining_buildings["area_sqm"]
    )
    mining_buildings_areas = (
        mining_buildings.groupby(["mining_id"])["C_GDP_building"].sum().reset_index()
    )
    mining_buildings = (
        mining_buildings.groupby(["osm_id"])["C_GDP_building"].sum().reset_index()
    )
    print("* GDP to buildings", mining_buildings["C_GDP_building"].sum())
    print(
        "* GDP ratio to buildings", mining_buildings["C_GDP_building"].sum() / tot_gpd
    )

    mining_areas = pd.merge(
        mining_areas, mining_buildings_areas, how="left", on=["mining_id"]
    )
    mining_areas["C_GDP_building"] = mining_areas["C_GDP_building"].fillna(0)
    mining_areas["GDP_building_ratio"] = mining_areas.progress_apply(
        lambda x: x["C_GDP_building"] / x["C_GDP"] if x["C_GDP"] > 0 else 0, axis=1
    )

    gpd.GeoDataFrame(
        mining_areas, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
    ).to_file(
        mining_gdp_output,
        layer="areas",
        driver="GPKG",
    )
    mining_buildings.to_csv(
        building_mining_gdp,
        index=False,
    )


if __name__ == "__main__":
    main()
