import connectivity
from pathlib import Path
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from rasterio.warp import reproject
from rasterio.enums import Resampling
import connectivity
import tifffile
import rioxarray as rxr
import os
# from osgeo import gdal
import rioxarray
import Robynlibrary


def rasterize_condition(geodf, out_file, resolution=10):
    # Determine bounds and resolution
    minx, miny, maxx, maxy = geodf.total_bounds
    width = int((maxx - minx) / resolution)
    height = int((maxy - miny) / resolution)
    transform = from_origin(minx, maxy, resolution, resolution)
    
    # Prepare shapes: tuple of (geometry, condition value)
    shapes = ((geom, value) for geom, value in zip(geodf.geometry, geodf['condition_factor_p4'])) #condition_factor_p4 used to be Condition
    
    # Rasterize the shapes
    condition_raster = rasterize(
        shapes=shapes,
        out_shape=(height, width),
        fill=0,  # No data value
        transform=transform,
        dtype='float32'
    )
    
    # Write the raster to a GeoTIFF file
    with rasterio.open(
        out_file,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype='float32',
        crs=geodf.crs,
        transform=transform,
    ) as dst:
        dst.write(condition_raster, 1)
    print(f"Raster saved at: {out_file}")
    
    # Return the raster array and its transform for further processing
    return condition_raster, transform

# Define a function to resample a raster and save it

def resample_and_save(source_raster, src_transform, src_crs, bounds, out_file, new_resolution=100):
    minx, miny, maxx, maxy = bounds
    new_width = int((maxx - minx) / new_resolution)
    new_height = int((maxy - miny) / new_resolution)
    new_transform = from_origin(minx, maxy, new_resolution, new_resolution)
    
    # Create an empty array to hold the resampled data
    resampled_raster = np.empty((new_height, new_width), dtype='float32')
    
    # Resample using the average method
    reproject(
        source=source_raster,
        destination=resampled_raster,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=new_transform,
        dst_crs=src_crs,
        resampling=Resampling.average
    )
    
    # Save the resampled raster to a GeoTIFF
    with rasterio.open(
        out_file,
        "w",
        driver="GTiff",
        height=new_height,
        width=new_width,
        count=1,
        dtype='float32',
        crs=src_crs,
        transform=new_transform,
    ) as dst:
        dst.write(resampled_raster, 1)
    print(f"100 m resolution raster saved at: {out_file}")




def open_raster_as_array(path):
    """
    Opens a TIFF raster file as a numpy array.
    
    Parameters:
        path (Path or str): Path to the TIFF raster.
        
    Returns:
        np.array: Raster data as a numpy array.
    """
    raster = tifffile.imread(path)
    array = np.array(raster)
    return array



def compute_connectivity(raster_path):
    """
    Opens a resampled raster, creates a binary land array,
    and computes the landscape connectivity.
    
    Parameters:
        raster_path (Path or str): Path to the resampled TIFF raster.
        
    Returns:
        connectivity_value: The computed connectivity value.
    """
    # Open the raster as an array
    condition_layer = open_raster_as_array(raster_path)
    
    # Create a binary array: 1 where condition_layer != 0, 0 elsewhere
    land_array = np.zeros_like(condition_layer, dtype='uint8')
    land_array[condition_layer != 0] = 1
    
    # Setup parameters for connectivity analysis
    n_processes = os.cpu_count()
    lambda_parameter = 5
    generations_mode = "one_generation"
    number_of_species_generations = 1
    
    # Compute connectivity (using your connectivity module)
    connectivity_value = connectivity.landscape_connectivity(
        condition_layer, 
        n_processes, 
        land_array, 
        lambda_parameter, 
        generations_mode, 
        number_of_species_generations
    )
    return connectivity_value


# --- Define a function to compute forest area ---
def compute_forest_area(gdf, classes, mixed=None):
    """
    Computes the total forest area from a GeoDataFrame based on the selected classes.
    For classes present in the 'mixed' dictionary, uses the provided forest fraction;
    otherwise, assumes a full forest (1.0) if the class is in the list.
    
    Parameters:
        gdf (GeoDataFrame): Must have a 'Classify' column and geometry.
        classes (list): Landcover classes to consider as forest.
        mixed (dict): Keys are mixed classes, values are the forest fraction (e.g., 0.25).
    
    Returns:
        float: Total forest area in m².
    """
    if mixed is None:
        mixed = {}
    # Create a series that assigns a forest fraction to each row:
    # - If the class is in the provided list and is in the mixed dict, use the mixed fraction.
    # - If it is in the list but not in the mixed dict, assign 1.
    # - Otherwise, assign 0.
    fraction_series = gdf['Classify'].apply(
        lambda x: mixed.get(x, 1.0) if x in classes else 0.0
    )
    # Calculate the forest area by multiplying each polygon's area by its fraction and summing.
    forest_area = (gdf.geometry.area * fraction_series).sum()
    return forest_area


def calc_connectivity_normalized(baseline, max_val, min_val):
    """
    Calculate the connectivity percentage based on a baseline value, 
    where max_val is 100% connectivity and min_val is 0% connectivity.
    
    Parameters:
        baseline (float): The baseline connectivity value.
        max_val (float): The maximum connectivity value (test ones).
        min_val (float): The minimum connectivity value (test zeros).
    
    Returns:
        float: The connectivity percentage of the baseline value.
    """
    normalized_connectivity = ((baseline - min_val) / (max_val - min_val)) * 100
    return normalized_connectivity




