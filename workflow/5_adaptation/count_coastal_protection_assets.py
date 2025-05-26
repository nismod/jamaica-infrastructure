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

def process_network_assets(network_info, flood_polygons, data_path):
    """Process a single network and count assets intersecting with flood polygons."""
    fname = os.path.join(data_path, network_info['path'])
    layer_type = network_info['gpkg_layer']
    ref = network_info['ref']
    
    # Initialize results dictionary
    results = {str(pid): 0 for pid in flood_polygons['id']}
    
    logging.info(f"Processing network '{ref}' (layer: {layer_type})")
    
    try:
        # Load the asset layer
        assets = gpd.read_file(fname, layer=layer_type)
        if assets.empty:
            logging.warning(f"No assets found in layer {layer_type}")
            return results
            
        logging.info(f"Found {len(assets)} assets in network '{ref}'")
        
        # Fix invalid geometries
        assets['geometry'] = assets['geometry'].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))
        flood_copy = flood_polygons.copy()
        flood_copy["geometry"] = flood_copy["geometry"].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))
        
        # Process each polygon
        for _, poly in flood_copy.iterrows():
            polygon_id = str(poly['id'])
            count = len(assets[assets.geometry.intersects(poly['geometry'])])
            results[polygon_id] = count
            
            if count > 0:
                logging.debug(f"Found {count} intersections for polygon {polygon_id}")
        
    except Exception as e:
        logging.error(f"Error processing network {ref}: {str(e)}")
        # logging.error(traceback.format_exc())
    
    total = sum(results.values())
    logging.info(f"Network '{ref}': found {total} total intersections")
    return results


def count_assets(RCP, RP, output, networks, data_path, flood_areas):
    """Count assets from each network intersecting with flood polygons."""
    output_dir = f'{output}/coastal_protection_assets/networks_to_protection_asset_overiew'
    os.makedirs(output_dir, exist_ok=True)
    
    for rcp in RCP:
        for rp in RP:
            flood_layer = f"flood_protection_area_rcp_{rcp}_rp{rp}"
            logging.info(f"Processing layer: {flood_layer}")
            
            try:
                # Load flood polygons
                flood_polygons = gpd.read_file(flood_areas, layer=flood_layer)
                
                if flood_polygons.empty or 'id' not in flood_polygons.columns:
                    logging.warning(f"No valid polygons or missing 'id' column in {flood_layer}")
                    continue
                
                # Initialize results
                all_results = {str(pid): {} for pid in flood_polygons['id']}
                
                # Process each network
                for _, network in networks.iterrows():
                    ref = network['ref']
                    network_results = process_network_assets(network, flood_polygons, data_path)
                    
                    # Add results to main dictionary
                    for pid, count in network_results.items():
                        all_results[pid][ref] = count
                
                # Create dataframe
                result_df = pd.DataFrame.from_dict(all_results, orient='index')
                result_df.index.name = 'polygon_id'
                result_df.reset_index(inplace=True)
                
                # Ensure all network columns exist
                for ref in networks['ref'].unique():
                    if ref not in result_df.columns:
                        result_df[ref] = 0
                
                # Save to CSV
                output_file = f'{output_dir}/coastal_protection_assets_breakdown_rcp_{rcp}_rp_{rp}.csv'
                result_df.to_csv(output_file, index=False)
                
                total_count = result_df.drop('polygon_id', axis=1).sum().sum()
                logging.info(f"Saved to {os.path.basename(output_file)} with {total_count} total intersections")
                
            except Exception as e:
                logging.error(f"Error processing layer {flood_layer}: {str(e)}")
                # logging.error(traceback.format_exc())

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
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to output gpkg",
)

def main(network_csv,processed_data_path,coastal_adaptation_assets,output_dir):
    RP = ['100']
    RCP = ['baseline2010', '262050', '262100', '452030', '452050', '452070', '452100', '852030', '852050', '852070', '852100']
    # RCP = ['baseline2010']

    data_path = processed_data_path
    
    networks_csv = network_csv
    networks = pd.read_csv(networks_csv)
    networks = networks[networks["ref"] != "buildings"]
    # print (networks)

    flood_areas = coastal_adaptation_assets
    output_path = output_dir

    count_assets(RCP, RP, output_path, networks, data_path, flood_areas)
    

if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )

    main()