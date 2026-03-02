"""Assign heights to buildings in Jamaica"""

import os
import glob
import geopandas as gpd
import pandas as pd
import numpy as np
import subprocess
import json
from tqdm import tqdm
import click

from jamaica_infrastructure.geo import LOCAL_PROJ_CRS_EPSG

tqdm.pandas()


@click.command()
@click.version_option("1.0")
@click.option(
    "--data-dir",
    "-d",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to processed data",
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to output gpkg",
)
@click.option(
    "--vector-raster-intersections-script",
    "-vri",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to subprocess script for computing vector raster intersections",
)
@click.option(
    "--building-layers",
    "-bl",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to building layers csv",
)
@click.option(
    "--building-rasters",
    "-br",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to building rasters csv",
)
@click.option(
    "--buildings-econ-activity",
    "-bea",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
def main(
    data_dir,
    output_dir,
    vector_raster_intersections_script,
    building_layers,
    building_rasters,
    buildings_econ_activity,
):
    # Set global paths
    processed_data_path = data_dir
    results_data_path = output_dir

    # Create a path to store intersection outputs
    output_path = os.path.join(results_data_path, "buildings_heights")
    if os.path.exists(output_path) == False:
        os.mkdir(output_path)

    vector_details_csv = building_layers
    raster_details_csv = building_rasters

    run_intersections = True  # Set to True is you want to run this process
    if run_intersections is True:
        args = [
            "python",
            f"{vector_raster_intersections_script}",
            f"{processed_data_path}",
            f"{vector_details_csv}",
            f"{raster_details_csv}",
            f"{output_path}",
        ]
        print("* Start the processing of buildings vector-raster intersections")
        print(args)
        subprocess.run(args, check=True)

    print("* Done with the processing of buildings vector-raster intersections")

    """Post-processing the buildings-raster intersection results
    """
    building_id_column = "osm_id"  # Building ID column
    raster_columns = pd.read_csv(raster_details_csv)["key"].values.tolist()
    # Read in intersection geoparquet
    intersections_file = [
        os.path.join(output_path, file)
        for file in os.listdir(output_path)
        if file.endswith(".geoparquet")
    ][0]
    if intersections_file:
        building_intersections = gpd.read_parquet(intersections_file)
        all_buildings = gpd.read_file(
            buildings_econ_activity,
            layer="areas",
        )
        for col in raster_columns:
            building_values = building_intersections[building_intersections[col] > 0]
            building_values["area_wt"] = building_values.geometry.area
            building_values[col] = building_values[col] * building_values["area_wt"]
            building_values = (
                building_values.groupby(building_id_column)[[col, "area_wt"]]
                .sum()
                .reset_index()
            )
            building_values[col] = building_values[col] / building_values["area_wt"]
            all_buildings = pd.merge(
                all_buildings,
                building_values[[building_id_column, col]],
                how="left",
                on=[building_id_column],
            ).fillna(0)
            del building_values

        del building_intersections

    gpd.GeoDataFrame(
        all_buildings, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    ).to_file(
        os.path.join(output_path, "buildings_heights.gpkg"),
        layer="areas",
        driver="GPKG",
    )
    print("* Done with assigning heights to buildings")


if __name__ == "__main__":
    main()
