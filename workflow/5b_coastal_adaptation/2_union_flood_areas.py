import logging

import click
import geopandas as gpd
from shapely.ops import unary_union


def combine_flood_protection_areas(rcps, rps, flood_areas, flood_layer):

    all_flood_areas = []
    for _ in rcps:
        for _ in rps:
            flood_polygons = gpd.read_file(flood_areas, layer=flood_layer)
            merged_flood_polygon = unary_union(flood_polygons.geometry)
            all_flood_areas.append(merged_flood_polygon)

    # Create a union of all merged flood areas
    final_flood_area = unary_union(all_flood_areas)

    return gpd.GeoDataFrame(geometry=[final_flood_area], crs=flood_polygons.crs)


@click.command()
@click.version_option("1.0")
@click.option(
    "--coastal-adaptation-assets",
    "-c",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to GPKG with coastal protection assets",
)
@click.option(
    "--rcp",
    "-rcp",
    required=True,
    type=float,
    help="Representative Concentration Pathway of flood raster",
)
@click.option(
    "--rp",
    "-rp",
    required=True,
    type=int,
    help="Return period of flood raster",
)
@click.option(
    "--epoch",
    "-ep",
    required=True,
    type=int,
    help="Epoch (rough year) of flood raster",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to output gpkg",
)
def main(coastal_adaptation_assets, rcp, rp, epoch, output_path):

    union_flood_area_gdf = combine_flood_protection_areas(
        [f"{int(rcp*10)}{epoch}"], [f"{rp}"], coastal_adaptation_assets, "areas"
    )

    union_flood_area_gdf.to_file(output_path, driver="GPKG")
    logging.info("Completed combining flood protection areas.")


if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )
    main()
