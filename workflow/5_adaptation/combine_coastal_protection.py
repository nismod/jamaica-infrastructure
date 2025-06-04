import logging
import geopandas as gpd
from shapely.ops import unary_union
import logging

import click

""" Defining Some Helper Functions utilised through the functions
    -------
"""

def join_flood_areas(flood_polygons):
    # Merge all polygons into a single polygon
    merged_polygon = unary_union(flood_polygons.geometry)
    return merged_polygon

""" Defing Main Logic Functions
    -------
"""

def combine_flood_protection_areas(RCP, RP, flood_areas, flood_layer):
    
    all_flood_areas = []

    for rcp in RCP:
        for rp in RP:
            # flood_polygon_layer = f"flood_protection_area_rcp_{rcp}_rp{rp}"
            flood_polygon_layer = flood_layer
            flood_polygons = gpd.read_file(flood_areas, layer=flood_polygon_layer)
            merged_flood_polygon = join_flood_areas(flood_polygons)
            all_flood_areas.append(merged_flood_polygon)
            

    # Create a union of all merged flood areas
    final_flood_area = unary_union(all_flood_areas)

    # Save the final merged flood polygon
    union_flood_area_gdf = gpd.GeoDataFrame(geometry=[final_flood_area], crs=flood_polygons.crs)
    return union_flood_area_gdf



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
    help="RPS value",
)

@click.option(
    "--rp",
    "-rp",
    required=True,
    type=int,
    help="RP value",
)

@click.option(
    "--epoch",
    "-ep",
    required=True,
    type=int,
    help="epcoh value",
)

@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to output gpkg",
)

def main(coastal_adaptation_assets,rcp, rp, epoch, output_dir):

    # RP = ['100']
    # RCP = ['baseline2010', '262050', '262100', '452030', '452050', '452070', '452100', '852030', '852050', '852070', '852100']
    
    RP = [f"{rp}"]

    fl_map_rcp = f"{int(rcp*10)}{epoch}"
    RCP = [fl_map_rcp]

    flood_area_layer = "areas"
    flood_areas = coastal_adaptation_assets

    union_flood_area_gdf = combine_flood_protection_areas(RCP, RP, flood_areas, flood_area_layer)

    union_flood_area_gdf.to_file(f"{output_dir}/coastal_protection_assets/combined_coastal_protection_area.gpkg", driver="GPKG")
    logging.info("Completed combining flood protection areas.")
 



if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )

    main()
