"""Process agricultural crop production to areas."""

import os
import subprocess

import click
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

from jamaica_infrastructure.geo import (
    LOCAL_PROJ_CRS_EPSG,
    create_voronoi_layer,
    raster_rewrite,
)


@click.command()
@click.version_option("1.0")
@click.option(
    "--spam-path",
    required=True,
    type=click.Path(exists=True, dir_okay=True, readable=True),
    help="Path to SPAM agriculture rasters directory.",
)
@click.option(
    "--spam-agriculture-outputs-path",
    required=True,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, writable=True),
    help="Path to SPAM agriculture outputs GeoPackage.",
)
def main(
    spam_path,
    spam_agriculture_outputs_path,
):
    """
    Process agricultural crop production to areas.

    Example usage:

        python workflow/1_context/agricultural_crops.py \
            --spam-path incoming_data/agriculture_data \
            --spam-agriculture-outputs-path processed_data/agriculture_data/spam_agriculture_outputs.gpkg
    """
    crop_folders = [
        "spam2010v2r0_global_val_prod_agg.geotiff",
        "spam2010v2r0_global_prod.geotiff",
        "spam2010v2r0_global_yield.geotiff",
    ]
    crop_strings = [
        "spam2010V2r0_global_V_agg_",
        "spam2010V2r0_global_P_",
        "spam2010V2r0_global_Y_",
    ]
    crop_outputs = ["production", "tonnage", "yield"]
    for crop_path, crop_field, crop_layer in zip(
        crop_folders, crop_strings, crop_outputs
    ):
        crop_data_path = os.path.join(spam_path, crop_path, "JAM")
        all_crops = []
        all_fields = []
        for file in os.listdir(crop_data_path):
            if file.endswith(".tif"):
                crop_file = file.replace(".tif", "")
                field_name = crop_file.replace(crop_field, "")
                all_fields.append(field_name)
                crop_raster_in = os.path.join(crop_data_path, f"{crop_file}.tif")
                crop_raster_out = os.path.join(
                    crop_data_path, f"{crop_file}_reproject.tif"
                )
                raster_rewrite(crop_raster_in, crop_raster_out)

                outCSVName = os.path.join(crop_data_path, f"{crop_file}.csv")
                subprocess.run(["gdal2xyz.py", "-csv", crop_raster_out, outCSVName])

                # Load points and convert to geodataframe with coordinates
                load_points = pd.read_csv(
                    outCSVName,
                    header=None,
                    names=["x", "y", field_name],
                    index_col=None,
                )
                # Ensure crop values are numeric and non-negative
                load_points[field_name] = (
                    load_points[field_name].fillna(0).clip(lower=0)
                )

                os.remove(outCSVName)
                os.remove(crop_raster_out)

                if len(all_crops) > 0:
                    all_crops = pd.merge(
                        all_crops, load_points, how="left", on=["x", "y"]
                    )
                else:
                    all_crops = load_points.copy()

                del load_points
                print("* Done with", file)

        all_crops["geometry"] = [Point(xy) for xy in zip(all_crops.x, all_crops.y)]
        crop_points = gpd.GeoDataFrame(
            all_crops, crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}", geometry="geometry"
        )
        crop_points["crop_id"] = crop_points.index.values.tolist()
        del all_crops
        print(crop_points)

        crop_areas = create_voronoi_layer(
            crop_points, "crop_id", epsg=LOCAL_PROJ_CRS_EPSG
        )

        crop_areas = gpd.GeoDataFrame(
            pd.merge(
                crop_areas,
                crop_points[["crop_id"] + all_fields],
                how="left",
                on=["crop_id"],
            ),
            geometry="geometry",
            crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}",
        )

        crop_points.to_file(
            spam_agriculture_outputs_path,
            layer=f"{crop_layer}_value",
            driver="GPKG",
        )
        crop_areas.to_file(
            spam_agriculture_outputs_path,
            layer=f"{crop_layer}_areas",
            driver="GPKG",
        )
        del crop_areas, crop_points


if __name__ == "__main__":
    main()
