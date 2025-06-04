import os
import math
import random
import numpy as np
import pandas as pd
import fiona
import logging
import fiona

import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPoint, MultiPolygon, LineString
from shapely.ops import unary_union, voronoi_diagram, nearest_points
from shapely import affinity

import logging
import warnings

import click

""" Defining Some Helper Functions utilised through the functions
    -------
"""

def find_intersecting_assets(flood_polygon, layer):
    # Ensure CRS match
    layer = layer.to_crs(flood_polygon.crs)

    # Perform a spatial overlay to get intersecting parts
    intersecting_assets = gpd.overlay(layer, flood_polygon, how="intersection")

    # Save to GPKG
    # intersecting_assets.to_file("test_intersect.gpkg", driver="GPKG")

    return intersecting_assets

def Assign_flood_area_to_asset(asset_network, RCP, RP, path, id_label, flood_areas, cost_col, flood_layer):
    def find_flood_area_asset_intersection(flood_polygons, asset):
        asset_geom = asset.geometry
        if not asset_geom.is_valid:
            asset_geom = asset_geom.buffer(0)
        flood_polygons = flood_polygons.copy()
        flood_polygons["geometry"] = flood_polygons["geometry"].apply(
            lambda geom: geom if geom.is_valid else geom.buffer(0)
        )
        intersecting_floods = flood_polygons[flood_polygons.geometry.intersects(asset_geom)]
        if not intersecting_floods.empty:
            max_flood_polygon = intersecting_floods.loc[intersecting_floods['max_flood_height'].idxmax()]
            flood_polygon_id = max_flood_polygon['id']
            max_flood_height = max_flood_polygon['max_flood_height']
        else:
            flood_polygon_id = None
            max_flood_height = None
        return flood_polygon_id, max_flood_height
    
    df = pd.read_parquet(path)
    all_layers = fiona.listlayers(flood_areas)
    flood_layers = {
        layer: gpd.read_file(flood_areas, layer=layer)
        for layer in all_layers if layer.startswith(flood_layer)
    }
    
    # Create a lookup dictionary for mean_cost from asset_network
    cost_lookup = dict(zip(asset_network[id_label], asset_network[cost_col]))
    
    previous_layer_name = None
    for i, row in df.iterrows():
        asset_id = row.iloc[0]
        
        # Add mean_cost for this asset if not already present
        if 'mean_cost' not in df.columns or pd.isna(df.at[i, 'mean_cost']):
            df.at[i, 'mean_cost'] = cost_lookup.get(asset_id, None)
        
        for rcp in RCP:
            for rp in RP:
                layer_name = flood_layer
                # if layer_name != previous_layer_name:
                # print(f"Current Layer: {layer_name}")
                previous_layer_name = layer_name
                flood_polygons = flood_layers.get(layer_name)
                asset = asset_network[asset_network[id_label] == asset_id].iloc[0]
                flood_id, flood_height = find_flood_area_asset_intersection(flood_polygons, asset)
                flood_id_col, flood_height_col = f"flood_id_rcp_{rcp}_rp_{rp}", f"flood_height_rcp_{rcp}_rp_{rp}"
                df.at[i, flood_id_col] = flood_id
                df.at[i, flood_height_col] = flood_height
        logging.info(f"Completed asset: {asset_id}")
    
    return df

""" Defing Main Logic Functions
    -------
"""

def filter_affected_assets(networks, union_flood_area_gdf, data_path, network_filter, output):
    # Assuming union_flood_area_gdf is already defined
    flood_polygon = union_flood_area_gdf

    # Take just the first (or main) geometry from the flood union
    flood = flood_polygon.loc[flood_polygon.index[0]]
    flood_polygon = gpd.GeoDataFrame([flood], geometry='geometry', crs=flood_polygon.crs)

    for index, n in networks.iterrows():
        fname = os.path.join(data_path, n['path'])
        id_col = n['asset_id_column']
        layer_type = n['asset_layer']
        # ref = n['asset_description']
        # ref = ref.replace(" ", "_")
        gpkg = n['asset_gpkg']
        layer = n['asset_layer']
        ref = f"{gpkg}_{layer}"
        

        if network_filter and ref not in network_filter:
            continue

        # Load the asset layer
        assets = gpd.read_file(fname, layer=layer_type)
        if assets.empty:
            continue

        # Find intersecting assets
        asset_intersections = find_intersecting_assets(flood_polygon, assets)
        
        
        # Keep only the ID column and drop duplicates
        filtered_ids = asset_intersections[[id_col]].drop_duplicates()

        # Save to Parquet
        
        parquet_path = os.path.join(f"{output}/coastal_protection_assets/network_protection_mappings", f"{ref}_coastal_filtered.parquet")
        os.makedirs(os.path.dirname(parquet_path), exist_ok=True)
        filtered_ids.to_parquet(parquet_path, index=False)
        
        logging.info(f"Completed filtering assets for network layer: {ref}")

def map_network_assets_to_protection(RCP, RP, output, networks, data_path, network_filter, flood_areas, flood_layer):
    for index, n in networks.iterrows():
        fname = os.path.join(data_path, n['path'])
        id_col = n['asset_id_column']
        layer_type = n['asset_layer']
        # ref = n['asset_description']
        # ref = ref.replace(" ", "_")
        cost_col = n['asset_mean_cost_column']
        gpkg = n['asset_gpkg']
        layer = n['asset_layer']
        ref = f"{gpkg}_{layer}"
        path = f'{output}/coastal_protection_assets/network_protection_mappings/{ref}_coastal_filtered.parquet'

        if network_filter and ref not in network_filter:
            continue

        logging.info(f"Processing assets for network layer: {ref}")

        # Load the asset layer
        assets = gpd.read_file(fname, layer=layer_type)
        if assets.empty:
            continue

        updated_output = Assign_flood_area_to_asset(assets, RCP, RP, path, id_col, flood_areas, cost_col, flood_layer)
        updated_output.to_parquet(path, index=False)


@click.command()
@click.version_option("1.0")

@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to network layers csv",
)
@click.option(
    "--processed-data-path",
    "-d",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to processed data",
)
@click.option(
    "--coastal-adaptation-assets",
    "-c",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to GPKG with coastal protection assets",
)
@click.option(
    "--combine-coastal-protection",
    "-cb",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to GPKG with combined coastal protectionareas",
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

def main(network_csv,
         processed_data_path,
         coastal_adaptation_assets,
         combine_coastal_protection,
         asset_gpkg,
         asset_layer,
         rcp,
         epoch,
         rp,
         output_dir
         ):
    # network_filter = ["transport_rail_edges"] #filter out certain network layers (for testing)
    network_filter = [f"{asset_gpkg}_{asset_layer}"]

    # RP = ['100']
    # RCP = ['baseline2010', '262050', '262100', '452030', '452050', '452070', '452100', '852030', '852050', '852070', '852100']
    
    RP = [f"{rp}"]

    fl_map_rcp = f"{int(rcp*10)}{epoch}"
    RCP = [fl_map_rcp]
    flood_area_layer = "areas"

    data_path = processed_data_path
    
    networks_csv = network_csv
    networks = pd.read_csv(networks_csv)
    networks = networks[networks["asset_description"] != "buildings"]
    # print (networks)

    flood_areas = coastal_adaptation_assets
    output_path = output_dir

    union_flood_area_gdf = gpd.read_file(combine_coastal_protection)
    filter_affected_assets(networks, union_flood_area_gdf, data_path, network_filter, output_path)
    map_network_assets_to_protection(RCP, RP, output_path, networks, data_path, network_filter, flood_areas, flood_area_layer)



if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )

    main()
