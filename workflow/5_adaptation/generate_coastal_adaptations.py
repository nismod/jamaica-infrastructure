import os
import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPoint, MultiPolygon, LineString, MultiLineString, GeometryCollection
from shapely.ops import unary_union, voronoi_diagram, nearest_points, linemerge
from shapely import affinity
from scipy.ndimage import label
from scipy.spatial import Voronoi
from scipy.optimize import minimize
from skimage import measure
from pathlib import Path

import rasterio
from rasterio.transform import from_origin
from rasterstats import zonal_stats
from shapely.strtree import STRtree

import sklearn.cluster

import logging
import os
import warnings

import click

""" Defining Some Helper Functions utilised through the functions
    -------
"""
def add_Layer_to_File(data, output_file, layer_name, driver):
    data.to_file(output_file, layer=layer_name, driver=driver)

# Adapted and modified from https://github.com/thomas-fred/jam-coastal-protection
def process_raster_to_clusters(threshold_m, eps, raster_file, minPts):
    """
    Processes a raster file to identify clusters of depth values using DBSCAN, 
    converts them to polygons, and returns a GeoDataFrame for further processing.

    Args:
        threshold_m (float): Minimum depth value to include pixels in the analysis.
        eps (float): DBSCAN epsilon parameter for cluster proximity.
        raster_file (str): File path to the input raster file.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing minimum enclosing polygons for clusters.
    """
    # Open the raster file and read the first band (depth values)
    raster = rasterio.open(raster_file)
    depth_m = raster.read(1)
    i, j = np.indices(depth_m.shape)  # Create row (i) and column (j) indices

    # Extract the transformation and CRS from the raster
    transform = from_origin(raster.bounds.left, raster.bounds.top, raster.res[0], raster.res[1])
    original_crs = raster.crs

    # Create a DataFrame and filter out pixels below the depth threshold
    df = pd.DataFrame(data={"i": i.ravel(), "j": j.ravel(), "depth_m": depth_m.ravel()})
    df = df[df["depth_m"] > threshold_m]  # Keep only pixels above the threshold

    if not df.empty:
        # Initialize the DBSCAN clustering algorithm and fit the filtered data
        db = sklearn.cluster.DBSCAN(eps=eps, min_samples=minPts)
        X = df.loc[:, ["i", "j"]].to_numpy()  # Extract pixel coordinates
        db.fit(X)

        # Retrieve cluster labels and count the number of clusters
        labels = db.labels_
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)  # Exclude noise points (-1)

        if n_clusters >= 2:
            # Create polygons from clusters and convert to a GeoDataFrame
            gdf = create_polygons_from_clusters(labels, depth_m.shape, transform, original_crs, df)

            # Reproject GeoDataFrame to EPSG:3448 coordinate reference system
            gdf = gdf.to_crs("EPSG:3448")

            # Generate minimum enclosing polygons for each cluster
            gdf = create_min_enclosing_polygons(gdf)
            return gdf

def create_polygons_from_clusters(labels, shape, transform, crs, df):
    """
    Converts DBSCAN cluster labels into polygons and organizes them in a GeoDataFrame.

    Args:
        labels (np.ndarray): Cluster labels assigned by DBSCAN (-1 for noise).
        shape (tuple): Shape of the raster (rows, cols).
        transform (Affine): Transformation to convert pixel indices to coordinates.
        crs (str): Coordinate reference system of the raster.
        df (pd.DataFrame): DataFrame containing pixel indices and other relevant data.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame containing polygons representing clusters.
    """
    # Create an array to store cluster labels mapped to raster pixels
    cluster_raster = np.full(shape, -1, dtype=np.int32)
    cluster_raster[df["i"], df["j"]] = labels

    polygons = []
    for cluster_label in np.unique(labels):
        if cluster_label == -1:
            continue  # Skip noise points

        # Identify the mask for the current cluster
        cluster_mask = cluster_raster == cluster_label
        contours = measure.find_contours(cluster_mask, level=0.5)

        # Convert contours to polygons
        for contour in contours:
            polygon_coords = [
                transform * (col, row) for row, col in contour
            ]
            polygon = Polygon(polygon_coords)
            if polygon.is_valid:
                polygons.append((cluster_label, polygon))

    # Create a GeoDataFrame from the polygons
    gdf = gpd.GeoDataFrame(polygons, columns=["cluster_label", "geometry"], crs=crs)
    return gdf

def create_min_enclosing_polygons(gdf):
    """
    Generates minimum enclosing polygons for each cluster in the GeoDataFrame.

    Args:
        gdf (gpd.GeoDataFrame): GeoDataFrame containing cluster polygons.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with updated minimum enclosing polygons.
    """
    enclosing_polygons = []

    for cluster_label in gdf["cluster_label"].unique():
        # Combine all polygons for the same cluster label into a single geometry
        cluster_polygons = gdf[gdf["cluster_label"] == cluster_label]
        combined_geometry = unary_union(cluster_polygons.geometry)

        # Create the minimum enclosing polygon (remove self-intersections with buffer)
        enclosing_polygon = combined_geometry.buffer(0)
        enclosing_polygons.append((cluster_label, enclosing_polygon))

    # Create a new GeoDataFrame with the enclosing polygons
    gdf_enclosing = gpd.GeoDataFrame(enclosing_polygons, columns=["cluster_label", "geometry"], crs=gdf.crs)
    return gdf_enclosing

def add_max_zonal_stats(geopkg, raster_path, new_column_name):
    """
    Add a new column to a GeoDataFrame with the maximum raster value for each polygon,
    mimicking the functionality of the Zonal Statistics Tool in QGIS. 
    The GeoDataFrame is reprojected to match the raster CRS during the calculation 
    and reprojected back to its original CRS afterward.

    Parameters:
    - geopkg (GeoDataFrame): GeoDataFrame containing polygons for which zonal statistics will be calculated.
    - raster_path (str): Path to the raster file used for calculating zonal statistics.
    - new_column_name (str): Name of the new column to store the maximum raster value for each polygon.

    Returns:
    - GeoDataFrame: The updated GeoDataFrame with the new column containing maximum raster values.
    """
    # Save the original CRS (Coordinate Reference System) of the GeoDataFrame for re-projection later
    original_crs = geopkg.crs

    # Open the raster file to retrieve its CRS
    with rasterio.open(raster_path) as src:
        raster_crs = src.crs

    # Reproject the GeoDataFrame to match the raster CRS if their CRS do not align
    if geopkg.crs != raster_crs:
        geopkg = geopkg.to_crs(raster_crs)


    geopkg_temp = geopkg.copy()
    # Ensure all geometries are valid (use a buffer of 0 as a fix for invalid geometries)
    geopkg_temp["geometry"] = geopkg_temp["geometry"].buffer(50)

    # Calculate zonal statistics to determine the maximum raster value within each polygon
    stats = zonal_stats(
        geopkg_temp,          # GeoDataFrame containing polygons
        raster_path,     # Path to the raster file
        stats="max",     # Calculate the maximum value for each zone
        geojson_out=False,  # Do not output results in GeoJSON format
        nodata=-9999     # Value to treat as NoData in the raster (adjust as necessary)
    )

    # Extract the maximum raster values from the zonal statistics results and add to a new column
    geopkg[new_column_name] = [
        max(stat["max"], 0) if stat["max"] is not None else None
        for stat in stats
    ]

    # Reproject the GeoDataFrame back to its original CRS
    if geopkg.crs != original_crs:
        geopkg = geopkg.to_crs(original_crs)

    return geopkg  # Return the updated GeoDataFrame with the new column



""" Generating Flood Zones Polygons via S-CAPE Method
    -------
"""
def generate_polygons_along_edge(edge_layer, node_layer, buffer):
    """
    Generate polygons (bounding boxes) aligned along edges in a GeoDataFrame.

    For each edge, this function creates a rectangular bounding box centered on the edge's midpoint, 
    rotated to align with the edge's direction, and expanded vertically by the specified buffer.

    Parameters:
    - edge_layer (GeoDataFrame): GeoDataFrame containing edges (LineString) with 'from_id' and 'to_id' attributes.
    - node_layer (GeoDataFrame): GeoDataFrame containing nodes (Point) with 'node_id' attributes.
    - buffer (float): Additional height (vertical buffer) to be added to the bounding box, in the km.

    Returns:
    - GeoDataFrame: A GeoDataFrame containing rotated and translated bounding box polygons for each coastline segment.
    """
    # List to store the resulting bounding boxes for each edge
    bounding_boxes = []

    # Lists to store additional attributes for the resulting GeoDataFrame
    distances = []  # Store the length of each edge
    midpoints = []  # Store the midpoint of each edge
    angles = []     # Store the rotation angle for each edge

    # Iterate through each edge in the edge layer
    for _, edge in edge_layer.iterrows():
        # Retrieve the 'from_id' and 'to_id' for the current edge
        from_node_id = edge['from_id']
        to_node_id = edge['to_id']
        
        # Get the coordinates of the start and end nodes using their IDs
        from_node = node_layer[node_layer['node_id'] == from_node_id].geometry.iloc[0]
        to_node = node_layer[node_layer['node_id'] == to_node_id].geometry.iloc[0]

        # Skip edges where the start and end nodes are at the same location
        if (round(from_node.x, 5) == round(to_node.x, 5)) and (round(from_node.y, 5) == round(to_node.y, 5)):
            continue

        # Calculate the straight-line distance (diagonal) between the two nodes
        distance = from_node.distance(to_node)
        distances.append(distance)

        # Calculate the midpoint of the edge
        midpoint = calculate_midpoint(from_node, to_node)
        midpoints.append(midpoint)

        # Calculate the angle of the edge relative to the horizontal axis
        angle = calculate_angle(from_node, to_node)
        angles.append(angle)

        # Define the height of the bounding box (including buffer)
        min_y = min(from_node.y, to_node.y)
        max_y = max(from_node.y, to_node.y)
        height = max_y - min_y + (buffer * 1000)  # Apply buffer to the height

        # Define the width of the bounding box (distance between the two nodes)
        min_x = min(from_node.x, to_node.x)
        max_x = min_x + distance  # Adjust width to match edge length

        # Create an initial rectangular bounding box
        bounding_box = Polygon([
            (min_x, min_y),
            (min_x, min_y + height),  # Top-left corner
            (max_x, min_y + height),  # Top-right corner
            (max_x, min_y),           # Bottom-right corner
            (min_x, min_y)            # Close back to bottom-left
        ])

        # Translate the bounding box to center it on the edge's midpoint
        translated_bounding_box = translate_bounding_box(bounding_box, midpoint)

        # Rotate the bounding box to align with the edge's orientation
        rotated_bounding_box = rotate_bounding_box(translated_bounding_box, midpoint, angle)

        # Append the final bounding box to the list
        bounding_boxes.append(rotated_bounding_box)

    # Convert the list of bounding boxes into a GeoDataFrame
    bounding_boxes_gdf = gpd.GeoDataFrame(geometry=bounding_boxes, crs=edge_layer.crs)

    # Add additional attributes (distance, angle) to the resulting GeoDataFrame
    bounding_boxes_gdf['id'] = range(len(bounding_boxes_gdf))  # Unique ID for each polygon
    edge_layer['rectangle_id'] = range(len(bounding_boxes_gdf))
    bounding_boxes_gdf['distance'] = distances  # Distance between the nodes
    bounding_boxes_gdf['angle'] = angles        # Angle of the edge in degrees

    return bounding_boxes_gdf, edge_layer

def calculate_midpoint(from_node, to_node):
    """
    Calculate the midpoint between two nodes.

    Parameters:
    - from_node (Point): The starting node as a Shapely Point.
    - to_node (Point): The ending node as a Shapely Point.

    Returns:
    - Point: The midpoint between the two nodes as a Shapely Point.
    """
    midpoint_x = (from_node.x + to_node.x) / 2
    midpoint_y = (from_node.y + to_node.y) / 2
    return Point(midpoint_x, midpoint_y)

def calculate_angle(from_node, to_node):
    """
    Calculate the angle of the line connecting two nodes relative to the horizontal axis.

    Parameters:
    - from_node (Point): The starting node as a Shapely Point.
    - to_node (Point): The ending node as a Shapely Point.

    Returns:
    - float: The angle in degrees, measured counterclockwise from the positive x-axis.
    """
    delta_x = to_node.x - from_node.x
    delta_y = to_node.y - from_node.y
    angle_rad = math.atan2(delta_y, delta_x)
    return math.degrees(angle_rad)

def translate_bounding_box(bounding_box, midpoint):
    """
    Translate a bounding box to center it on a specified midpoint.

    Parameters:
    - bounding_box (Polygon): The bounding box polygon to translate.
    - midpoint (Point): The target center point for translation.

    Returns:
    - Polygon: The translated bounding box polygon.
    """
    current_center_x = (bounding_box.bounds[0] + bounding_box.bounds[2]) / 2
    current_center_y = (bounding_box.bounds[1] + bounding_box.bounds[3]) / 2
    translate_x = midpoint.x - current_center_x
    translate_y = midpoint.y - current_center_y
    return affinity.translate(bounding_box, xoff=translate_x, yoff=translate_y)

def rotate_bounding_box(bounding_box, midpoint, angle):
    """
    Rotate a bounding box around a specified midpoint by a given angle.

    Parameters:
    - bounding_box (Polygon): The bounding box polygon to rotate.
    - midpoint (Point): The center point for rotation.
    - angle (float): The angle in degrees by which to rotate the bounding box.

    Returns:
    - Polygon: The rotated bounding box polygon.
    """
    return affinity.rotate(bounding_box, angle, origin=(midpoint.x, midpoint.y))

def join_polygons_by_flood_height(gdf, flood_height_threshold, length_limit, group_size, full_coastline_layer):
    """
    Groups polygons in a GeoDataFrame based on their maximum flood height, with additional grouping 
    conditions based on total distance/width and maximum group size. 
    
    Polygons are joined together if their flood height exceeds a specified threshold, but groups
    are finalized when the total distance exceeds the specified limit or when switching between
    above/below threshold states.

    Parameters:
    - gdf (GeoDataFrame): Input GeoDataFrame containing polygons with "distance" column for width
    - flood_height_threshold (float): The flood height threshold above which polygons are grouped together
    - length_limit (float): The total distance limit for groups; groups exceeding this limit are finalized
    - group_size (int): The minimum number of polygons that can make up a group
    - full_coastline_layer: The coastline layer (kept for compatibility, not used in distance calculation)

    Returns:
    - GeoDataFrame: A new GeoDataFrame containing the grouped polygons
    """
    
    # Initialize lists to store results for each group
    groups = []  # List of combined polygons (grouped geometries)
    heights = []  # List of maximum flood heights for each group
    ids = []  # List of original polygon IDs for each group
    total_distances = []  # List of total distances for each group
    
    # Initialize variables to track the current group being processed
    current_group = []  # List of polygons in the current group
    current_heights = []  # List of max flood heights for the current group
    current_ids = []  # List of original IDs for the current group
    current_distances = []  # List of distances for the current group
    is_below_threshold = None  # Flag indicating whether the current polygon is below the threshold

    def finalize_current_group():
        """Helper function to finalize the current group"""
        if not current_group:
            return
            
        # Calculate total distance for the current group
        total_distance = sum(current_distances)
            
        if len(current_group) < group_size and groups:
            # Add the current polygons to the previous group if current group is too small
            groups[-1] = unary_union([groups[-1]] + current_group).convex_hull
            heights[-1] = max(heights[-1], max(current_heights))
            ids[-1].extend(current_ids)
            # Add distances to the previous group's total
            total_distances[-1] += total_distance
        else:
            # Finalize the group normally by combining geometries
            combined_geometry = unary_union(current_group)
            groups.append(combined_geometry.convex_hull)
            heights.append(max(current_heights))
            ids.append(current_ids)
            total_distances.append(total_distance)

    # Iterate over each polygon in the GeoDataFrame
    for index, row in gdf.iterrows():
        polygon_above_threshold = row["max_flood_height"] > flood_height_threshold
        
        # Check if threshold state is changing - if so, finalize current group
        if is_below_threshold is not None and is_below_threshold != (not polygon_above_threshold):
            finalize_current_group()
            # Reset the current group and associated data
            current_group = []
            current_heights = []
            current_ids = []
            current_distances = []
        
        # Check if adding this polygon would exceed the length limit
        if current_group:
            test_total_distance = sum(current_distances) + row["distance"]
            
            if test_total_distance > length_limit:
                # Finalize current group before adding this polygon
                finalize_current_group()
                # Reset the current group and associated data
                current_group = []
                current_heights = []
                current_ids = []
                current_distances = []

        # Update threshold state and add polygon to current group
        is_below_threshold = not polygon_above_threshold
        current_group.append(row["geometry"])
        current_heights.append(row["max_flood_height"])
        current_ids.append(row["id"])
        current_distances.append(row["distance"])

    # Finalize the last group if any polygons remain
    finalize_current_group()

    # Convert the grouped data back into a GeoDataFrame
    result_gdf = gpd.GeoDataFrame(
        {
            "geometry": gpd.GeoSeries(groups),  # Combined geometries for each group
            "max_flood_height": heights,  # Max flood heights for each group
            "original_ids": ids,  # Original IDs for each group
            "total_distance": total_distances,  # Total distances for each group
        }
    )

    # Add a unique ID field to the result GeoDataFrame
    result_gdf["id"] = range(len(result_gdf))

    # Set the CRS (Coordinate Reference System) to match the original GeoDataFrame
    if not result_gdf.empty:
        result_gdf.set_crs(gdf.crs, inplace=True)

    return result_gdf

def order_edges_around_island(edge_layer):
    """
    Order edges around an island, starting with edge_0 and following the
    connections based on from_id and to_id.

    Parameters:
    - edge_layer (GeoDataFrame): GeoDataFrame containing edges with 'id', 'from_id', and 'to_id' columns

    Returns:
    - List[str]: List of associated bounding box IDs in the order the edges are traced around the island
    """
    # Initialize the ordered list with the first edge (edge_0)
    ordered_edges = []

    # Create a dictionary to lookup edges by their 'from_id' for efficient access
    edges_by_from_id = {
        edge['from_id']: edge
        for _, edge in edge_layer.iterrows()
    }

    # Start with the edge having 'id' equal to 'edge_0'
    current_edge = edge_layer[edge_layer['id'] == 'edge_0'].iloc[0]
    ordered_edges.append(current_edge['rectangle_id'])  # Append the first edge's ID

    # Track the starting node to detect when we've completed a loop
    start_node = current_edge['from_id']

    # Continue looping through the edges, following the 'to_id' until we return to the starting node
    while True:
        # Find the next edge based on the current edge's 'to_id'
        next_from_id = current_edge['to_id']

        # Break the loop if we return to the starting node (complete the loop)
        if next_from_id == start_node:
            break

        # Look up the next edge using the 'from_id' as the key in the dictionary
        current_edge = edges_by_from_id[next_from_id]
        ordered_edges.append(current_edge['rectangle_id'])  # Add the next edge's ID to the ordered list

    return ordered_edges  # Return the ordered list of edge IDs

def find_coast_segment_for_polygon(coastline, polygon):
    # if coastline.geometry.iloc[0].intersects(polygon.geometry):
    intersection = coastline.geometry.iloc[0].intersection(polygon)
    return intersection

def find_max_coast_length(coast_segment):
    # Check if the geometry is empty
    if coast_segment.is_empty:
        return 0  # or return some other value if appropriate

    if isinstance(coast_segment, MultiLineString):
        # Use .geoms to iterate through individual LineStrings within the MultiLineString
        max_length = max(coast.length for coast in coast_segment.geoms)
    elif isinstance(coast_segment, LineString):
        # If it's a LineString, return its length
        max_length = coast_segment.length
    elif isinstance(coast_segment, GeometryCollection):
        # Handle GeometryCollection case (can be multiple different geometries)
        max_length = max(find_max_coast_length(geom) for geom in coast_segment)
    else:
        raise TypeError("Input must be a LineString, MultiLineString, or GeometryCollection.")

    return max_length



""" Refining Flood Zones
    -------
"""
def merge_overlapping_polygons(gdf, overlap_threshold):
    """
    Iteratively merges overlapping polygons that overlap by more than the specified threshold.
    
    Args:
        gdf (GeoDataFrame): GeoDataFrame containing polygons with 'geometry' and 'max_flood_height' fields.
        overlap_threshold (float): The minimum percentage of overlap required to trigger a merge (e.g., 0.35 for 35% overlap).
    
    Returns:
        GeoDataFrame: A new GeoDataFrame with merged polygons until no pair overlaps by more than the threshold.
    """
    
    # Ensure that the CRS (Coordinate Reference System) is defined for the input data
    if gdf.crs is None:
        gdf.set_crs("EPSG:3097", inplace=True)  # Set CRS to Jamaica Metric Grid, adjust if necessary

    # Helper function to perform one pass of merging overlapping polygons
    def merge_once(gdf, overlap_threshold):
        """
        Performs one pass over the GeoDataFrame to merge overlapping polygons.
        Merges polygons whose intersection area exceeds the given overlap threshold.
        
        Args:
            gdf (GeoDataFrame): GeoDataFrame to process.
            overlap_threshold (float): The overlap threshold to trigger a merge.
        
        Returns:
            GeoDataFrame: A new GeoDataFrame with merged polygons from this pass.
        """
        merged = []  # List to hold merged geometries
        indices_to_merge = set()  # Set to keep track of merged polygons

        # Loop through all polygons in the GeoDataFrame to check for overlaps
        for i, geom1 in enumerate(gdf.geometry):
            if i in indices_to_merge:
                continue  # Skip geometries that are already merged

            for j, geom2 in enumerate(gdf.geometry):
                if i >= j or j in indices_to_merge:
                    continue  # Avoid redundant checks and already-merged geometries

                # Check if the two polygons intersect
                if geom1.intersects(geom2):
                    intersection = geom1.intersection(geom2)
                    if intersection.is_empty:
                        continue  # Skip if the intersection is empty

                    # Calculate the area of the intersection and compare with the polygons' areas
                    area1 = geom1.area
                    area2 = geom2.area
                    intersection_area = intersection.area

                    # Check if the overlap exceeds the threshold
                    if (intersection_area / min(area1, area2)) > overlap_threshold:
                        # Merge the two polygons by creating a convex hull around the union
                        new_geom = unary_union([geom1, geom2]).convex_hull
                        new_height = max(gdf.loc[i, "max_flood_height"], gdf.loc[j, "max_flood_height"])

                        # Store the merged geometry and its max flood height
                        merged.append({
                            "geometry": new_geom,
                            "max_flood_height": new_height
                        })
                        indices_to_merge.update([i, j])  # Mark the merged polygons

                        break  # Stop checking other polygons for the current geometry once merged

        # Retain the geometries that were not merged
        non_merged = gdf.loc[~gdf.index.isin(indices_to_merge)]

        # Combine the merged geometries with the non-overlapping ones
        merged_gdf = gpd.GeoDataFrame(
            {
                "geometry": [m["geometry"] for m in merged],
                "max_flood_height": [m["max_flood_height"] for m in merged],
            }
        )

        # Ensure CRS is explicitly set for the merged geometries
        merged_gdf.set_crs(gdf.crs, inplace=True)

        # Concatenate the non-overlapping geometries with the merged ones
        final_gdf = gpd.GeoDataFrame(pd.concat([non_merged, merged_gdf], ignore_index=True))

        # Set CRS again after concatenation to ensure consistency
        final_gdf.set_crs(gdf.crs, inplace=True)

        return final_gdf

    # Iteratively merge polygons until no more overlaps are found
    prev_gdf = gdf
    while True:
        new_gdf = merge_once(prev_gdf, overlap_threshold)
        if len(new_gdf) == len(prev_gdf):  # No changes (no more overlaps)
            break
        prev_gdf = new_gdf  # Update the GeoDataFrame for the next iteration

    # Ensure final CRS consistency
    prev_gdf.set_crs(gdf.crs, inplace=True)

    # Reassign unique IDs starting from 1 for the merged polygons
    prev_gdf = prev_gdf.reset_index(drop=True)  # Reset the index
    prev_gdf['id'] = prev_gdf.index + 1  # Assign new unique IDs based on the index

    prev_gdf["length"] = prev_gdf.geometry.apply(lambda geom: longest_bounding_box_side(geom, prev_gdf.crs))

    return prev_gdf

def get_intersections_by_id(gdf, id_column="id"):
    """
    Creates a dictionary where the key is the polygon ID, and the value is a list of IDs
    of polygons that intersect with it.

    Args:
        gdf (GeoDataFrame): GeoDataFrame containing polygons with 'geometry'.
        id_column (str): The column name containing unique polygon IDs.

    Returns:
        dict: A dictionary where keys are polygon IDs and values are lists of intersecting polygon IDs.
    """
    intersections = {}

    # Ensure the CRS is defined
    if gdf.crs is None:
        raise ValueError("GeoDataFrame must have a defined CRS.")

    # Iterate over all pairs of polygons
    for i, geom1 in gdf.iterrows():
        polygon_id = geom1[id_column]
        intersecting_ids = []

        for j, geom2 in gdf.iterrows():
            other_id = geom2[id_column]
            
            # Skip self-comparison
            if polygon_id == other_id:
                continue
            
            # Check if the polygons intersect
            if geom1.geometry.intersects(geom2.geometry):
                intersecting_ids.append(other_id)

        # Update the dictionary
        intersections[polygon_id] = intersecting_ids

    return intersections

def clip_larger_polygons(gdf, id_column="id"):
    """
    Clips larger polygons to remove intersections with smaller polygons based on the intersection dictionary.

    Args:
        gdf (GeoDataFrame): GeoDataFrame containing polygons with 'geometry'.
        id_column (str): The column name containing unique polygon IDs.

    Returns:
        GeoDataFrame: Updated GeoDataFrame with clipped polygons.
    """
    # Ensure the CRS is defined
    if gdf.crs is None:
        raise ValueError("GeoDataFrame must have a defined CRS.")

    # Create the intersection dictionary
    intersections = get_intersections_by_id(gdf, id_column=id_column)
    
    # Iterate through the dictionary
    for polygon_id, intersecting_ids in intersections.items():
        # Get the geometry of the current polygon
        geom1 = gdf.loc[gdf[id_column] == polygon_id, 'geometry'].values[0]

        # Check each intersecting polygon
        for other_id in intersecting_ids:
            # Get the geometry of the other polygon
            geom2 = gdf.loc[gdf[id_column] == other_id, 'geometry'].values[0]

            # Check if they still intersect
            if geom1.intersects(geom2):
                # Calculate the areas of the polygons
                area1 = geom1.area
                area2 = geom2.area

                # Clip the larger polygon by subtracting the intersection
                if area1 >= area2:
                    clipped_geom = geom1.difference(geom2)
                    geom1 = extract_polygons(clipped_geom)  # Ensure only polygons
                else:
                    clipped_geom = geom2.difference(geom1)
                    geom2 = extract_polygons(clipped_geom)  # Ensure only polygons

                # Update the geometries in the GeoDataFrame
                gdf.loc[gdf[id_column] == polygon_id, 'geometry'] = geom1
                gdf.loc[gdf[id_column] == other_id, 'geometry'] = geom2

    # Reset index for consistency
    gdf = gdf.reset_index(drop=True)
    return gdf

def extract_polygons(geometry):
    """
    Extracts only polygonal geometries from a GeometryCollection or MultiPolygon.

    Args:
        geometry (shapely.geometry): Input geometry.

    Returns:
        shapely.geometry.Polygon or MultiPolygon: Cleaned geometry.
    """
    if isinstance(geometry, Polygon) or isinstance(geometry, MultiPolygon):
        return geometry
    elif geometry.is_empty:
        return None  # Return None for empty geometries
    elif hasattr(geometry, 'geoms'):  # Handle GeometryCollection
        polygons = [geom for geom in geometry.geoms if isinstance(geom, (Polygon, MultiPolygon))]
        if len(polygons) == 1:
            return polygons[0]  # Return as a single Polygon
        elif len(polygons) > 1:
            return MultiPolygon(polygons)  # Return as MultiPolygon
    return None  # In case no valid polygon is found

def longest_bounding_box_side(polygon, crs):
    """
    Computes the longest side of the minimum bounding box in meters.
    Ensures the bounding box covers the entire MultiPolygon.

    Parameters:
    - polygon: shapely.geometry.Polygon or MultiPolygon
    - crs: Coordinate Reference System (from the source GeoDataFrame)

    Returns:
    - The longest side of the bounding box in meters.
    """

    # Ensure the CRS is projected
    if crs is None or not crs.is_projected:
        raise ValueError("The provided CRS must be projected with units in meters.")

    # If it's a MultiPolygon, merge into a single shape
    if isinstance(polygon, MultiPolygon):
        polygon = polygon.convex_hull  # Create a convex hull around all parts

    # Get the minimum rotated rectangle
    min_rect = polygon.minimum_rotated_rectangle
    coords = list(min_rect.exterior.coords)  # Extract rectangle coordinates

    # Compute side lengths
    side_lengths = [((coords[i][0] - coords[i+1][0])**2 + (coords[i][1] - coords[i+1][1])**2) ** 0.5 
                    for i in range(4)]  # Only 4 sides in a rectangle

    return max(side_lengths)

# Adapted and modified from https://github.com/thomas-fred/jam-coastal-protection
def combine_floods_within_scape(voronoi_polygons: gpd.GeoDataFrame,flood_polygons: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    combined_polygons = []  # List to store resulting geometries
    voronoi_ids = []  # List to store Voronoi IDs (optional, for tracking)

    # Loop through each Voronoi polygon
    for voronoi_index, voronoi in voronoi_polygons.iterrows():
        voronoi_geom = voronoi.geometry

        # Find all flood polygons that intersect with the current Voronoi polygon
        intersecting_floods = flood_polygons[flood_polygons.intersects(voronoi_geom)]

        if not intersecting_floods.empty:
            # Combine the intersecting flood polygons within the Voronoi polygon
            combined_geom = intersecting_floods.intersection(voronoi_geom).union_all()

            if combined_geom:
                # Append the filtered geometry and its Voronoi index
                combined_polygons.append(combined_geom)
                voronoi_ids.append(voronoi['id'])
        else:
            combined_polygons.append(None)
            voronoi_ids.append(voronoi['id'])

    # Create a GeoDataFrame for the combined polygons
    result_gdf = gpd.GeoDataFrame({
        "id": voronoi_ids,  # Include the Voronoi IDs
        "geometry": combined_polygons  # Include the geometries
    }, crs=voronoi_polygons.crs)  # Set CRS to match Voronoi polygons

    return result_gdf  # Return the combined flood polygons GeoDataFrame

""" 
    -------
"""

def find_midpoints(geometry, num_splits):
    """Finds the midpoints of `num_splits` segments along a LineString or MultiLineString."""
    if not isinstance(geometry, (LineString, MultiLineString)):
        raise ValueError("Input must be a LineString or MultiLineString")
    
    total_length = geometry.length
    if total_length == 0:
        return []

    segment_length = total_length / num_splits
    return [geometry.interpolate((i + 0.5) * segment_length, normalized=False) for i in range(num_splits)]

def generate_extra_points(polygon, coastline, num_pts, scale_factor):
    """Generates points at segment midpoints from a split coastline within a polygon."""
    coast_intersection = coastline.intersection(polygon)
    if num_pts < 2:
        num_splits = 2
    else:
        num_splits = num_pts
    
        
    if coast_intersection.empty or coast_intersection.isna().all():
        return []
    
    if isinstance(coast_intersection, gpd.GeoSeries):
        coast_intersection = coast_intersection.dropna().reset_index(drop=True)
        if coast_intersection.empty:
            return []
        coast_segment = coast_intersection.iloc[0]  # Take the first valid geometry
    else:
        coast_segment = coast_intersection
    
    if not isinstance(coast_segment, (LineString, MultiLineString)) or coast_segment.length == 0:
        return []
    
    # t_scale = 10_000_000
    # print(polygon.area/scale_factor)
    # coast_segment = coast_segment.simplify(polygon.area/scale_factor)
    

    # If the coast_segment is a MultiPolygon
    if isinstance(coast_segment, MultiLineString):
        # Calculate the total length of the MultiLineString
        total_length = coast_segment.length

        # Filter out LineStrings with less than 10% of the total length
        large_lines = [line for line in coast_segment.geoms if line.length >= total_length * 0.1]

        # Create a new MultiLineString with the remaining lines
        coast_segment = MultiLineString(large_lines)
    
    # segments = split_geometry(coast_segment, num_splits)
    midpoints = find_midpoints(coast_segment, num_splits)
    
    layers = {
        "test_line": [coast_segment],
        # "test_coast": [segments],
        "test_points": midpoints,
        "test_poly": [polygon]
    }
    
    # for name, geometries in layers.items():
    #     try:
    #         gdf = gpd.GeoDataFrame(geometry=geometries, crs=jamaica_polygon_revised.crs)
    #         add_Layer_to_File(gdf, processing_file, name, "GPKG")
    #     except:
    #         pass
    
    return midpoints

def get_larger_polygon(jamaica_polygon):
    """Returns the larger polygon (Jamaica's convex hull)."""
    convex_gdf = jamaica_polygon
    return convex_gdf.iloc[0].geometry

def compute_centroids(intersections, smaller_polygons):
    """Compute centroids for the smaller polygons (coastline polygons)."""
    def get_centroid(row, smaller_polygons):
        if row.geometry and not row.geometry.is_empty:
            return row.geometry.centroid
        match = smaller_polygons[smaller_polygons["id"] == row["id"]]
        return match.geometry.iloc[0].centroid if not match.empty else None
    
    return intersections.apply(lambda row: get_centroid(row, smaller_polygons), axis=1)

def create_centroid_gdf(centroids, smaller_polygons):
    """Create a GeoDataFrame with centroids and associated centroid IDs."""
    if 'id' in smaller_polygons.columns and not smaller_polygons['id'].isnull().any():
        centroid_ids = smaller_polygons['id'].values
    else:
        centroid_ids = range(len(centroids))

    return gpd.GeoDataFrame({'c_id': centroid_ids, 's_id': centroid_ids }, geometry=centroids, crs=smaller_polygons.crs)

def adjust_centroids_to_polygon(centroids, larger_polygon):
    """Adjust centroids to be within or touching the larger polygon's boundary."""
    adjusted_centroids = [adjust_point_to_polygon(pt, larger_polygon) for pt in centroids]
    return adjusted_centroids

def adjust_point_to_polygon(point, polygon):
    """Adjust a point to be within or touching the boundary of a polygon."""
    if polygon.contains(point) or polygon.touches(point):
        return point
    else:
        nearest_boundary_point = nearest_points(point, polygon.boundary)[1]
        adjusted_point = Point(
            point.x + (nearest_boundary_point.x - point.x),
            point.y + (nearest_boundary_point.y - point.y)
        )
        if polygon.contains(adjusted_point) or polygon.touches(adjusted_point):
            return adjusted_point
        return nearest_boundary_point

def create_voronoi_tessellation(centroids):
    """Create Voronoi tessellation based on centroids."""
    seed_points = MultiPoint(centroids)
    return voronoi_diagram(seed_points)

def clip_voronoi_to_polygon(voronoi, larger_polygon):
    """Clip Voronoi regions to the boundary of a larger polygon."""
    return [region.intersection(larger_polygon) for region in voronoi.geoms]

def associate_centroids_to_regions(valid_regions, centroids, centroid_gdf):
    """Associate centroids to the valid Voronoi regions based on proximity."""
    associated_centroid_ids = []
    associated_polygon_ids = []
    for region in valid_regions:
        region_centroid = region.centroid
        nearest_centroid = centroids.distance(region_centroid).idxmin()
        associated_centroid_ids.append(centroid_gdf.loc[nearest_centroid, 'c_id'])
        associated_polygon_ids.append(centroid_gdf.loc[nearest_centroid, 's_id'])
    return associated_centroid_ids, associated_polygon_ids

def create_clipped_voronoi_gdf(valid_regions, associated_centroid_ids, associated_polygon_ids, crs):
    """Create a GeoDataFrame for the clipped Voronoi regions with associated centroid IDs."""
    return gpd.GeoDataFrame(
        {'id': range(len(valid_regions)),
         'flood_zone_id': associated_polygon_ids,
         'centroid_id': associated_centroid_ids},
        geometry=valid_regions,
        crs=crs
    )

def find_coastline_length(coastline, polygon):
    """Find the length of the coastline segment that intersects with the polygon."""
    intersection = coastline.geometry.iloc[0].intersection(polygon)
    return intersection.length

def compute_bbox_length(coasline_polygon, voronoi_polygon, crs):
    intersection = coasline_polygon.intersection(voronoi_polygon)
    
    if intersection.is_empty:
        return 0  # No intersection found
    
    bounding_box = intersection.convex_hull.minimum_rotated_rectangle

    # Extract the coordinates of the bounding box
    coords = list(bounding_box.exterior.coords)

    # Compute distances between consecutive points
    distances = [np.sqrt((coords[i][0] - coords[i+1][0])**2 + (coords[i][1] - coords[i+1][1])**2) 
                 for i in range(len(coords) - 1)]
    
    longest_side = max(distances)  # Longest side of the bounding box
    
    # Save to file for debugging
    # add_Layer_to_File(gpd.GeoDataFrame({'geometry': [bounding_box]}, crs=crs), processing_file, "bbox", "GPKG")

    return longest_side
    
def create_voronoi(jamaica_polygon, jam_contour, coastline_polygons, intersections, length_limit, coastline):
    """
    Create Voronoi tessellation based on the centroids of smaller polygons (coastline),
    clipped to the boundary of a larger polygon (Jamaica's convex hull).
    Uses a combined approach with multiple strategies to guarantee length constraint compliance.
    
    Args:
        jamaica_polygon (GeoDataFrame): The larger polygon (e.g., Jamaica's convex hull).
        jam_contour (GeoDataFrame): A GeoDataFrame containing the exact contour of the Jamaica.
        coastline_polygons (GeoDataFrame): The smaller polygons (e.g., coastline polygons) for which Voronoi regions will be calculated.
        intersections (GeoDataFrame): Intersections for computing centroids.
        length_limit (float): Maximum allowed coastal length for any Voronoi region.
        coastline (GeoDataFrame): The coastline geometry.

    Returns:
        GeoDataFrame: Clipped Voronoi polygons inside the larger polygon with associated centroid IDs.
        GeoDataFrame: Centroids of the smaller polygons used for Voronoi tessellation with a 'c_id' column.
    """
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point, Polygon
    from shapely.ops import unary_union, voronoi_diagram
    from shapely.geometry import MultiPoint
    
    logging.info("   [4.0] Starting Combined Voronoi Refinement Strategy")
    
    # Initial setup
    combined_coastline_polygons = unary_union(coastline_polygons.geometry)
    larger_polygon = get_larger_polygon(jamaica_polygon)
    smaller_polygons = coastline_polygons
    smaller_polygons.crs = coastline_polygons.crs
    
    # Calculate scale factor for point generation
    areas = coastline_polygons.geometry.area.values
    scale_factor = np.mean(areas) / len(areas)
    
    # Helper function: process_voronoi (from your original code)
    def process_voronoi(centroids, centroid_gdf, larger_polygon, smaller_polygons):
        """Process Voronoi tessellation with the given centroids"""
        # Adjust centroids to be within the larger polygon
        adjusted_centroids = adjust_centroids_to_polygon(centroids, jam_contour.iloc[0].geometry)
        
        # Convert the adjusted centroids back to a GeoSeries
        adjusted_centroids_gs = gpd.GeoSeries(adjusted_centroids, crs=centroid_gdf.crs)
        
        # Update the GeoDataFrame with adjusted centroids
        centroid_gdf.geometry = adjusted_centroids_gs
        
        # Create Voronoi tessellation based on centroids
        voronoi = create_voronoi_tessellation(adjusted_centroids)
        
        # Clip Voronoi regions to the boundary of the larger polygon
        clipped_regions = clip_voronoi_to_polygon(voronoi, larger_polygon)
        
        # Filter out invalid or empty regions, keeping only valid polygons
        valid_regions = [region for region in clipped_regions if not region.is_empty and isinstance(region, Polygon)]
        
        # Match centroids to Voronoi regions based on proximity
        associated_centroid_ids, associated_polygon_ids = associate_centroids_to_regions(valid_regions, centroids, centroid_gdf)
        
        # Create GeoDataFrame for clipped Voronoi regions with associated centroid IDs
        clipped_gdf = create_clipped_voronoi_gdf(valid_regions, associated_centroid_ids, associated_polygon_ids, smaller_polygons.crs)
        
        return clipped_gdf, centroid_gdf
    
    # Strategy 1: Pre-calculation Strategy
    def calculate_required_points_strategy(coastline_polygons, coastline, length_limit):
        """Pre-calculate exactly how many points each coastline polygon needs"""
        required_points = {}
        
        for idx, row in coastline_polygons.iterrows():
            polygon_id = row['id']
            polygon = row.geometry
            
            # Get coastline intersection with this polygon
            try:
                coast_intersection = coastline.intersection(polygon)
                if hasattr(coast_intersection, 'length') and coast_intersection.length > 0:
                    # Calculate required number of segments to meet length limit
                    required_segments = max(2, int(np.ceil(coast_intersection.length / length_limit)))
                    required_points[polygon_id] = required_segments
                else:
                    required_points[polygon_id] = 2
            except:
                required_points[polygon_id] = 2
        
        return required_points

    def generate_all_points_upfront(coastline_polygons, coastline, length_limit, scale_factor):
        """Generate all required seed points upfront based on pre-calculated requirements"""
        required_points = calculate_required_points_strategy(coastline_polygons, coastline, length_limit)
        
        all_centroids = []
        centroid_data = []
        c_id_counter = 0
        
        for polygon_id, num_points in required_points.items():
            polygon_row = coastline_polygons[coastline_polygons['id'] == polygon_id]
            if polygon_row.empty:
                continue
                
            polygon = polygon_row.iloc[0].geometry
            
            # Generate the required number of points
            new_points = generate_extra_points(polygon, coastline, num_points, scale_factor)
            
            if new_points:
                for point in new_points:
                    centroid_data.append({
                        'c_id': c_id_counter,
                        's_id': polygon_id,
                        'geometry': point
                    })
                    c_id_counter += 1
                all_centroids.extend(new_points)
            else:
                # Fallback to centroid if no points generated
                centroid_point = polygon.centroid
                centroid_data.append({
                    'c_id': c_id_counter,
                    's_id': polygon_id,
                    'geometry': centroid_point
                })
                c_id_counter += 1
                all_centroids.append(centroid_point)
        
        # Create GeoDataFrame with all centroids
        if centroid_data:
            centroid_gdf = gpd.GeoDataFrame(centroid_data, crs=coastline_polygons.crs)
        else:
            # Fallback to original centroids if something went wrong
            centroids = compute_centroids(intersections, smaller_polygons)
            centroid_gdf = create_centroid_gdf(centroids, smaller_polygons)
        
        return centroid_gdf, all_centroids

    # Strategy 2: Hierarchical Subdivision
    def hierarchical_subdivision_strategy(clipped_gdf, centroid_gdf, coastline_polygons, 
                                        coastline, length_limit, larger_polygon, 
                                        smaller_polygons, scale_factor):
        """Use hierarchical approach to subdivide all problematic regions simultaneously"""
        max_depth = 10
        current_depth = 0
        
        logging.info("   [4.2] Applying Hierarchical Subdivision")
        
        while current_depth < max_depth:
            # Calculate coastal lengths for all regions
            clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(
                lambda poly: find_coastline_length(coastline, poly)
            )
            
            # Find ALL regions that exceed the limit
            problematic_regions = clipped_gdf[clipped_gdf["coastal_length"] > length_limit]
            
            if problematic_regions.empty:
                logging.info(f"   Hierarchical subdivision complete at depth {current_depth}")
                break
            
            logging.info(f"   Depth {current_depth}: Subdividing {len(problematic_regions)} regions")
            
            # Group by flood_zone_id to process each coastline polygon
            problematic_by_zone = problematic_regions.groupby('flood_zone_id')
            
            new_points_batch = []
            zones_to_remove = set()
            
            for zone_id, zone_regions in problematic_by_zone:
                zones_to_remove.add(zone_id)
                
                # Get the coastline polygon
                polygon_row = coastline_polygons[coastline_polygons['id'] == zone_id]
                if polygon_row.empty:
                    continue
                    
                polygon = polygon_row.iloc[0].geometry
                
                # Calculate how many MORE points we need
                max_length_in_zone = zone_regions["coastal_length"].max()
                current_points_in_zone = len(zone_regions)
                
                # Calculate multiplication factor needed
                subdivision_factor = max(2, int(np.ceil(max_length_in_zone / length_limit)))
                new_total_points = max(current_points_in_zone * 2, subdivision_factor)
                
                # Generate new points for this zone
                new_points = generate_extra_points(polygon, coastline, new_total_points, scale_factor)
                
                if new_points:
                    for point in new_points:
                        new_points_batch.append({
                            'c_id': len(centroid_gdf) + len(new_points_batch),
                            's_id': zone_id,
                            'geometry': point
                        })
            
            # Remove old centroids for problematic zones
            centroid_gdf = centroid_gdf[~centroid_gdf['s_id'].isin(zones_to_remove)].reset_index(drop=True)
            
            # Add all new points at once
            if new_points_batch:
                new_centroid_gdf = gpd.GeoDataFrame(new_points_batch, crs=centroid_gdf.crs)
                centroid_gdf = pd.concat([centroid_gdf, new_centroid_gdf], ignore_index=True)
            
            # Regenerate Voronoi with all new points
            centroids = centroid_gdf.geometry
            clipped_gdf, centroid_gdf = process_voronoi(
                centroids, centroid_gdf, larger_polygon, smaller_polygons
            )
            
            current_depth += 1
        
        return clipped_gdf, centroid_gdf

    # Strategy 3: Constraint-Based Refinement
    def constraint_based_refinement(clipped_gdf, centroid_gdf, coastline_polygons, 
                                   coastline, length_limit, larger_polygon, 
                                   smaller_polygons, scale_factor):
        """Constraint-based approach with progress monitoring"""
        previous_max_length = float('inf')
        stagnation_count = 0
        max_stagnation = 3
        iteration_count = 0
        
        logging.info("   [4.3] Applying Constraint-Based Refinement")
        
        while True:
            # Calculate current state
            clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(
                lambda poly: find_coastline_length(coastline, poly)
            )
            
            current_max_length = clipped_gdf["coastal_length"].max()
            regions_exceeding = len(clipped_gdf[clipped_gdf["coastal_length"] > length_limit])
            
            # Check if we're done
            if regions_exceeding == 0:
                logging.info(f"   Constraint satisfaction achieved after {iteration_count} iterations")
                break
            
            # Check for progress
            if current_max_length >= previous_max_length:
                stagnation_count += 1
                if stagnation_count >= max_stagnation:
                    logging.info("   Stagnation detected, applying aggressive subdivision")
                    clipped_gdf, centroid_gdf = apply_aggressive_subdivision(
                        clipped_gdf, centroid_gdf, coastline_polygons, coastline,
                        length_limit, larger_polygon, smaller_polygons, scale_factor
                    )
                    stagnation_count = 0
            else:
                stagnation_count = 0
            
            # Apply standard refinement step
            clipped_gdf, centroid_gdf = apply_standard_refinement_step(
                clipped_gdf, centroid_gdf, coastline_polygons, coastline,
                length_limit, larger_polygon, smaller_polygons, scale_factor
            )
            
            previous_max_length = current_max_length
            iteration_count += 1
            
            logging.info(f"   Iteration {iteration_count}: Max length = {current_max_length:.1f}, "
                  f"Violations = {regions_exceeding}")
            
            # Safety check - if we've done too many iterations, apply final aggressive subdivision
            if iteration_count > 20:
                logging.info("   Applying final aggressive subdivision")
                clipped_gdf, centroid_gdf = apply_final_aggressive_subdivision(
                    clipped_gdf, centroid_gdf, coastline_polygons, coastline,
                    length_limit, larger_polygon, smaller_polygons, scale_factor
                )
                break
        
        return clipped_gdf, centroid_gdf

    def apply_aggressive_subdivision(clipped_gdf, centroid_gdf, coastline_polygons, 
                                   coastline, length_limit, larger_polygon, 
                                   smaller_polygons, scale_factor):
        """Apply aggressive subdivision to break through stagnation"""
        
        problematic_regions = clipped_gdf[clipped_gdf["coastal_length"] > length_limit]
        
        for zone_id in problematic_regions['flood_zone_id'].unique():
            zone_regions = problematic_regions[problematic_regions['flood_zone_id'] == zone_id]
            max_length = zone_regions["coastal_length"].max()
            
            # Calculate aggressive subdivision factor
            subdivision_factor = max(5, int(np.ceil(max_length / (length_limit * 0.3))))
            
            # Remove old points
            centroid_gdf = centroid_gdf[centroid_gdf['s_id'] != zone_id].reset_index(drop=True)
            
            # Add many new points
            polygon_row = coastline_polygons[coastline_polygons['id'] == zone_id]
            if not polygon_row.empty:
                polygon = polygon_row.iloc[0].geometry
                new_points = generate_extra_points(polygon, coastline, subdivision_factor, scale_factor)
                
                if new_points:
                    new_data = []
                    for point in new_points:
                        new_data.append({
                            'c_id': len(centroid_gdf) + len(new_data),
                            's_id': zone_id,
                            'geometry': point
                        })
                    
                    new_centroid_gdf = gpd.GeoDataFrame(new_data, crs=centroid_gdf.crs)
                    centroid_gdf = pd.concat([centroid_gdf, new_centroid_gdf], ignore_index=True)
        
        # Regenerate Voronoi
        centroids = centroid_gdf.geometry
        clipped_gdf, centroid_gdf = process_voronoi(
            centroids, centroid_gdf, larger_polygon, smaller_polygons
        )
        
        return clipped_gdf, centroid_gdf

    def apply_standard_refinement_step(clipped_gdf, centroid_gdf, coastline_polygons, 
                                     coastline, length_limit, larger_polygon, 
                                     smaller_polygons, scale_factor):
        """Apply one step of standard refinement"""
        
        # Find the longest region
        longest_polygon = clipped_gdf.loc[clipped_gdf["coastal_length"].idxmax()]
        
        if longest_polygon['coastal_length'] <= length_limit:
            return clipped_gdf, centroid_gdf
        
        # Get the problematic zone
        scape_id = longest_polygon['flood_zone_id']
        polygon_row = coastline_polygons[coastline_polygons['id'] == scape_id]
        
        if polygon_row.empty:
            return clipped_gdf, centroid_gdf
            
        polygon = polygon_row.iloc[0].geometry
        
        # Remove existing points for this zone
        centroid_gdf = centroid_gdf[centroid_gdf['s_id'] != scape_id].reset_index(drop=True)
        
        # Calculate new points needed
        current_regions = clipped_gdf[clipped_gdf["flood_zone_id"] == scape_id]
        num_pts = len(current_regions) + 2  # Add more points
        
        # Generate new points
        new_points = generate_extra_points(polygon, coastline, num_pts, scale_factor)
        
        if new_points:
            new_data = []
            for point in new_points:
                new_data.append({
                    'c_id': len(centroid_gdf) + len(new_data),
                    's_id': scape_id,
                    'geometry': point
                })
            
            new_centroid_gdf = gpd.GeoDataFrame(new_data, crs=centroid_gdf.crs)
            centroid_gdf = pd.concat([centroid_gdf, new_centroid_gdf], ignore_index=True)
            
            # Regenerate Voronoi
            centroids = centroid_gdf.geometry
            clipped_gdf, centroid_gdf = process_voronoi(
                centroids, centroid_gdf, larger_polygon, smaller_polygons
            )
        
        return clipped_gdf, centroid_gdf

    def apply_final_aggressive_subdivision(clipped_gdf, centroid_gdf, coastline_polygons, 
                                         coastline, length_limit, larger_polygon, 
                                         smaller_polygons, scale_factor):
        """Final aggressive subdivision to guarantee constraint satisfaction"""
        
        logging.info("   Applying final aggressive subdivision to guarantee constraints")
        
        # Find all violating regions
        clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(
            lambda poly: find_coastline_length(coastline, poly)
        )
        
        problematic_regions = clipped_gdf[clipped_gdf["coastal_length"] > length_limit]
        
        for zone_id in problematic_regions['flood_zone_id'].unique():
            zone_regions = problematic_regions[problematic_regions['flood_zone_id'] == zone_id]
            max_length = zone_regions["coastal_length"].max()
            
            # Very aggressive subdivision - guarantee we have enough points
            subdivision_factor = max(10, int(np.ceil(max_length / (length_limit * 0.2))))
            
            # Remove old points
            centroid_gdf = centroid_gdf[centroid_gdf['s_id'] != zone_id].reset_index(drop=True)
            
            # Get polygon and generate many points
            polygon_row = coastline_polygons[coastline_polygons['id'] == zone_id]
            if not polygon_row.empty:
                polygon = polygon_row.iloc[0].geometry
                new_points = generate_extra_points(polygon, coastline, subdivision_factor, scale_factor)
                
                if new_points:
                    new_data = []
                    for point in new_points:
                        new_data.append({
                            'c_id': len(centroid_gdf) + len(new_data),
                            's_id': zone_id,
                            'geometry': point
                        })
                    
                    new_centroid_gdf = gpd.GeoDataFrame(new_data, crs=centroid_gdf.crs)
                    centroid_gdf = pd.concat([centroid_gdf, new_centroid_gdf], ignore_index=True)
        
        # Final Voronoi generation
        centroids = centroid_gdf.geometry
        clipped_gdf, centroid_gdf = process_voronoi(
            centroids, centroid_gdf, larger_polygon, smaller_polygons
        )
        
        # Final check
        clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(
            lambda poly: find_coastline_length(coastline, poly)
        )
        
        return clipped_gdf, centroid_gdf

    # Combined Strategy Implementation
    try:
        # Step 1: Try pre-calculation strategy
        logging.info("   [4.1] Attempting Pre-calculation Strategy")
        centroid_gdf, all_centroids = generate_all_points_upfront(
            coastline_polygons, coastline, length_limit, scale_factor
        )
        
        # Generate initial Voronoi
        clipped_gdf, centroid_gdf = process_voronoi(
            centroid_gdf.geometry, centroid_gdf, larger_polygon, smaller_polygons
        )
        
        # Check if this solved the problem
        clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(
            lambda poly: find_coastline_length(coastline, poly)
        )
        
        violations = len(clipped_gdf[clipped_gdf["coastal_length"] > length_limit])
        
        if violations == 0:
            logging.info("   [4.1] Pre-calculation strategy successful!")
            return clipped_gdf, centroid_gdf
        else:
            logging.info(f"   [4.1] Pre-calculation left {violations} violations, trying hierarchical")
            
    except Exception as e:
        logging.info(f"   [4.1] Pre-calculation failed: {e}, falling back to original method")
        # Fallback to original centroid computation
        centroids = compute_centroids(intersections, smaller_polygons)
        centroid_gdf = create_centroid_gdf(centroids, smaller_polygons)
        clipped_gdf, centroid_gdf = process_voronoi(centroids, centroid_gdf, larger_polygon, smaller_polygons)
        clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(lambda poly: find_coastline_length(coastline, poly))

    # Step 2: Try hierarchical subdivision
    try:
        clipped_gdf, centroid_gdf = hierarchical_subdivision_strategy(
            clipped_gdf, centroid_gdf, coastline_polygons, coastline, 
            length_limit, larger_polygon, smaller_polygons, scale_factor
        )
        
        # Check if hierarchical solved it
        clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(
            lambda poly: find_coastline_length(coastline, poly)
        )
        violations = len(clipped_gdf[clipped_gdf["coastal_length"] > length_limit])
        
        if violations == 0:
            logging.info("   [4.2] Hierarchical subdivision successful!")
            return clipped_gdf, centroid_gdf
        else:
            logging.info(f"   [4.2] Hierarchical left {violations} violations, using constraint-based")
            
    except Exception as e:
        logging.info(f"   [4.2] Hierarchical subdivision failed: {e}")

    # Step 3: Final fallback to constraint-based approach
    logging.info("   [4.3] Using constraint-based approach as final guarantee")
    clipped_gdf, centroid_gdf = constraint_based_refinement(
        clipped_gdf, centroid_gdf, coastline_polygons, coastline,
        length_limit, larger_polygon, smaller_polygons, scale_factor
    )
    
    # Final validation and reporting
    clipped_gdf["coastal_length"] = clipped_gdf.geometry.apply(
        lambda poly: find_coastline_length(coastline, poly)
    )
    
    final_violations = len(clipped_gdf[clipped_gdf["coastal_length"] > length_limit])
    max_length = clipped_gdf["coastal_length"].max()
    mean_length = clipped_gdf["coastal_length"].mean()
    
    logging.info(f"   [4.4] Final Results:")
    logging.info(f"   Total regions: {len(clipped_gdf)}")
    logging.info(f"   Regions exceeding limit: {final_violations}")
    logging.info(f"   Max coastal length: {max_length:.2f}")
    logging.info(f"   Mean coastal length: {mean_length:.2f}")
    logging.info(f"   Length limit: {length_limit}")
    
    if final_violations > 0:
        logging.info(f"   WARNING: {final_violations} regions still exceed the length limit!")
        violating_regions = clipped_gdf[clipped_gdf["coastal_length"] > length_limit]
        logging.info(f"   Violating lengths: {violating_regions['coastal_length'].tolist()}")
    else:
        logging.info("   SUCCESS: All regions comply with length limit!")
    
    return clipped_gdf, centroid_gdf

# Adapted and modified from https://github.com/thomas-fred/jam-coastal-protection
def combine_floods_within_voronoi(voronoi_polygons: gpd.GeoDataFrame,flood_polygons: gpd.GeoDataFrame,voronoi_points: gpd.GeoDataFrame,pts_id_field:str ,max_distance: float) -> gpd.GeoDataFrame:
    """
    Combines all parts of flood polygons that intersect with each Voronoi polygon,
    but remove any parts whose centroids are farther than a given distance from
    the corresponding Voronoi point.

    Args:
        voronoi_polygons (gpd.GeoDataFrame): GeoDataFrame containing the Voronoi polygons.
        flood_polygons (gpd.GeoDataFrame): GeoDataFrame containing the flood polygons to combine.
        voronoi_points (gpd.GeoDataFrame): GeoDataFrame containing the Voronoi points with a 'fid' column that matches the 'centroid_id' column.
        max_distance (float): Maximum allowed distance (in the same CRS units) between a polygon's centroid and the corresponding Voronoi point.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the combined flood polygons for each Voronoi polygon,
                           with the flood areas merged based on their proximity to the Voronoi points.
    """
    combined_polygons = []  # List to store resulting geometries
    voronoi_ids = []  # List to store Voronoi IDs (optional, for tracking)

    # Loop through each Voronoi polygon
    for voronoi_index, voronoi in voronoi_polygons.iterrows():
        voronoi_geom = voronoi.geometry
        centroid_id = voronoi["centroid_id"]  # Get the Voronoi polygon's centroid ID

        # Get the corresponding Voronoi point based on the centroid ID
        voronoi_point = voronoi_points[voronoi_points[pts_id_field] == centroid_id]
        if voronoi_point.empty:
            continue  # Skip if no matching Voronoi point is found

        voronoi_point_geom = voronoi_point.geometry.iloc[0]

        # Find all flood polygons that intersect with the current Voronoi polygon
        intersecting_floods = flood_polygons[flood_polygons.intersects(voronoi_geom)]

        if not intersecting_floods.empty:
            # Combine the intersecting flood polygons within the Voronoi polygon
            combined_geom = intersecting_floods.intersection(voronoi_geom).union_all()

            # Filter out parts of the combined geometry that are too far from the Voronoi point
            if isinstance(combined_geom, MultiPolygon):
                # If the combined geometry is a MultiPolygon, check each part
                filtered_parts = [
                    part for part in combined_geom.geoms
                    if part.centroid.distance(voronoi_point_geom) <= max_distance
                ]
                combined_geom = MultiPolygon(filtered_parts) if filtered_parts else None
            elif combined_geom.centroid.distance(voronoi_point_geom) > max_distance:
                # If the combined geometry is a single polygon, check its centroid distance
                combined_geom = None

            if combined_geom:
                # Append the filtered geometry and its Voronoi index
                combined_polygons.append(combined_geom)
                voronoi_ids.append(voronoi_index)

    # Create a GeoDataFrame for the combined polygons
    result_gdf = gpd.GeoDataFrame({
        "id": voronoi_ids,  # Include the Voronoi IDs
        "geometry": combined_polygons  # Include the geometries
    }, crs=voronoi_polygons.crs)  # Set CRS to match Voronoi polygons

    return result_gdf  # Return the combined flood polygons GeoDataFrame

def linestring_intersect_polygons(default_inland_distance: float, coast: gpd.GeoDataFrame, polygons: gpd.GeoDataFrame, flood_polygons: gpd.GeoDataFrame) -> (gpd.GeoDataFrame, gpd.GeoDataFrame):
    """
    Sections the coastline based on intersections with each Voronoi polygon and individually 
    buffers the new linestring segments until they completely enclose the flood polygon 
    within the corresponding Voronoi section.

    Also returns a GeoDataFrame of intersections between the linestring and Voronoi polygons before any buffering.

    Args:
        default_inland_distance (float): Buffer radius in kilometers.
        coast (gpd.GeoDataFrame): GeoDataFrame containing the coastline geometry to be buffered.
        polygons (gpd.GeoDataFrame): GeoDataFrame of Voronoi polygons to check intersections.
        flood_polygons (gpd.GeoDataFrame): GeoDataFrame of flood polygons with a 'voronoi_id' column to associate with Voronoi polygons.

    Returns:
        tuple: 
            - GeoDataFrame containing the intersections of the linestring and Voronoi polygons before buffering.
            - GeoDataFrame containing the final intersected and buffered geometries, with Voronoi IDs.
    """
    # Prepare lists to store intersections and final buffered geometries
    initial_intersections = []
    final_intersections = []

    # Create a GeoDataFrame containing the union of all coastline geometries
    linestring = gpd.GeoDataFrame({"id": [0], "geometry": [coast.geometry.union_all()]})
    linestring.crs = coast.crs  # Set CRS to match the coast's CRS

    # Iterate over each Voronoi polygon
    for idx, poly in polygons.iterrows():
        # Check if the linestring intersects with the Voronoi polygon
        if linestring.geometry.iloc[0].intersects(poly.geometry):
            # Perform the intersection between the linestring and the Voronoi polygon
            intersection = linestring.geometry.iloc[0].intersection(poly.geometry)
            # print(poly['id'])
            
            # Append the raw intersection with the associated Voronoi ID
            initial_intersections.append({
                "geometry": intersection,
                "id": poly['id']
            })

            # Buffer the intersection piece separately based on the provided inland distance
            buffered_intersection = intersection.buffer(default_inland_distance * 1_000)
            
            # Find the corresponding flood polygon that has the same 'id'
            flood_polygon = flood_polygons[flood_polygons['id'] == idx]
            
            if not flood_polygon.empty:
                flood_geometry = flood_polygon.geometry.iloc[0]
                
                # Start iteratively buffering the intersection from 0.5 km and increase by 0.1 km
                current_buffer_radius = default_inland_distance * 1_000  # Start with a buffer radius of 0.5 km
                while not buffered_intersection.contains(flood_geometry):
                    # Increase buffer by 0.1 km at each iteration until it contains the flood geometry
                    current_buffer_radius += 0.1 * 1_000
                    buffered_intersection = intersection.buffer(current_buffer_radius)
            
            # Perform the final intersection with the flood polygon
            final_intersection = buffered_intersection.intersection(poly.geometry)
            
            # Append the final intersection with its associated Voronoi ID to the list
            final_intersections.append({
                "geometry": final_intersection,
                "id": poly['id']
            })
    
    # Create GeoDataFrames for initial and final intersections
    initial_intersections_gdf = gpd.GeoDataFrame(initial_intersections, crs=polygons.crs)
    final_intersections_gdf = gpd.GeoDataFrame(final_intersections, crs=polygons.crs)

    return initial_intersections_gdf, final_intersections_gdf

def split_multipart_polygons(gdf):
    new_features = []
    separated_parts = []  # Track only the separated small parts
    small_sections = []
    
    for _, row in gdf.iterrows():
        geom = row.geometry
        if isinstance(geom, MultiPolygon):  # Check if it's a multi-part polygon
            parts = list(geom.geoms)
            largest_part = max(parts, key=lambda p: p.area)  # Identify the largest part
            for part in parts:
                new_row = row.copy()
                new_row.geometry = part
                separated_parts.append(new_row)
                if part != largest_part:  # Ignore the largest part
                    new_row = row.copy()
                    new_row.geometry = part
                    small_sections.append(new_row)
        else:
            new_features.append(row)
    
    new_gdf = gpd.GeoDataFrame(new_features + separated_parts, columns=gdf.columns, crs=gdf.crs)
    new_gdf["id"] = range(1, len(new_gdf) + 1)  # Reassign IDs sequentially
    
    
    
    # Use a more robust method to identify small sections
    small_sections = new_gdf[new_gdf.geometry.apply(lambda geom: any(geom.equals(part.geometry) for part in small_sections))]
    touching_pairs = set()
    
    # Identify touching polygons and store as ordered pairs
    for _, poly in small_sections.iterrows():
        for _, other_poly in new_gdf.iterrows():
            if poly["id"] != other_poly["id"] and poly.geometry.buffer(0.5).intersects(other_poly.geometry.buffer(0.5)):
                pair = (poly["id"], other_poly["id"]) \
                    if poly.geometry.area < other_poly.geometry.area else \
                    (other_poly["id"], poly["id"])
                touching_pairs.add(pair)
    
    touching_pairs = list(touching_pairs)
    area_dict = {row["id"]: row.geometry.area for _, row in new_gdf.iterrows()}
    
    # Sort touching pairs by polygon area
    touching_pairs.sort(key=lambda x: area_dict[x[0]])

    buffer_distance = 0.5
    sorted_pairs = []

    best_pairs = {}  # Dictionary to store the best (first_id, second_id) pair

    for first_id, second_id in touching_pairs:
        first_polygon = new_gdf.loc[new_gdf['id'] == first_id, 'geometry'].values[0]
        second_polygon = new_gdf.loc[new_gdf['id'] == second_id, 'geometry'].values[0]
        
        shared_border_length = first_polygon.buffer(buffer_distance).intersection(second_polygon).length  

        # If first_id is not in best_pairs OR this pair has a longer shared border, update it
        if first_id not in best_pairs or shared_border_length > best_pairs[first_id][1]:
            best_pairs[first_id] = (second_id, shared_border_length)

    # Convert dictionary back to a list of tuples
    sorted_pairs = [(first_id, second_id) for first_id, (second_id, _) in best_pairs.items()]
    # print (sorted_pairs)

    merged_gdf = new_gdf.copy()
    srtdp = sorted_pairs.copy()
    merged_polygons = []
    
    i = 0
    while i < len(srtdp):
        first_id, second_id = srtdp[i]
        
        poly1 = merged_gdf.loc[merged_gdf["id"] == first_id, "geometry"].values[0]
        poly2 = merged_gdf.loc[merged_gdf["id"] == second_id, "geometry"].values[0]
        
        new_polygon = poly1.union(poly2.buffer(buffer_distance))
        new_id = first_id
        
        merged_gdf.loc[merged_gdf["id"] == first_id, "geometry"] = new_polygon
        merged_polygons.append(second_id)
        
        srtdp = [(new_id if x == second_id else x, new_id if y == second_id else y) for x, y in srtdp]
        i += 1
    
    merged_gdf = merged_gdf[~merged_gdf['id'].isin(merged_polygons)].reset_index(drop=True)

    return merged_gdf

def merge_adjacent_polygons_by_coastline(gdf, coastline, length_limit):
    # Ensure coastline is a single (multi)line geometry
    if isinstance(coastline, gpd.GeoDataFrame) or isinstance(coastline, gpd.GeoSeries):
        coastline = coastline.union_all()

    buffer_distance = 0.5
    gdf = gdf.copy()
    gdf["buffered"] = gdf.geometry.buffer(buffer_distance)
    gdf["area"] = gdf.geometry.area
    area_dict = dict(zip(gdf["id"], gdf["area"]))
    geom_dict = dict(zip(gdf["id"], gdf.geometry))
    buffered_dict = dict(zip(gdf["id"], gdf["buffered"]))

    sindex = gdf.sindex
    touching_pairs = set()

    for idx, row in gdf.iterrows():
        possible_matches_index = list(sindex.intersection(row["buffered"].bounds))
        for other_idx in possible_matches_index:
            other_row = gdf.iloc[other_idx]
            if row["id"] != other_row["id"] and row["buffered"].intersects(other_row["buffered"]):
                pair = (row["id"], other_row["id"]) if row["area"] < other_row["area"] else (other_row["id"], row["id"])
                touching_pairs.add(pair)

    touching_pairs = list(touching_pairs)
    max_area = np.mean(gdf["area"].values) / 2

    touching_pairs.sort(key=lambda x: area_dict[x[0]])
    sorted_pairs = []

    for first_id, second_id in touching_pairs:
        if area_dict[first_id] > max_area:
            continue

        poly1 = geom_dict[first_id]
        poly2 = geom_dict[second_id]
        buf1 = buffered_dict[first_id]
        buf2 = buffered_dict[second_id]

        shared_border_length = buf1.intersection(poly2).length
        total_coastline_length = (
            buf1.intersection(coastline).length +
            buf2.intersection(coastline).length
        )

        if total_coastline_length <= length_limit:
            sorted_pairs.append((first_id, second_id, shared_border_length))

    sorted_pairs.sort(key=lambda x: (area_dict[x[0]], -x[2]))
    srtdp = [(x[0], x[1]) for x in sorted_pairs]

    merged_ids = set()
    i = 0
    while i < len(srtdp):
        first_id, second_id = srtdp[i]
        if first_id in merged_ids or second_id in merged_ids:
            i += 1
            continue

        poly1 = geom_dict[first_id]
        poly2 = geom_dict[second_id]
        buf1 = poly1.buffer(buffer_distance)
        buf2 = poly2.buffer(buffer_distance)

        length1 = buf1.intersection(coastline).length
        length2 = buf2.intersection(coastline).length

        new_length = length1 + length2

        if new_length < length_limit:
            new_geom = poly1.union(buf2)
            geom_dict[first_id] = new_geom
            buffered_dict[first_id] = new_geom.buffer(buffer_distance)
            area_dict[first_id] = new_geom.area
            merged_ids.add(second_id)

            srtdp = [(first_id if x == second_id else x, first_id if y == second_id else y) for x, y in srtdp]

        i += 1

    merged_gdf = gdf[~gdf["id"].isin(merged_ids)].copy()
    merged_gdf["geometry"] = merged_gdf["id"].map(geom_dict)
    merged_gdf = merged_gdf.drop(columns=["buffered", "area"]).reset_index(drop=True)
    return merged_gdf   


""" Main Function that calls the processes
    -------
"""
def Generate_Coastal_Flood_Protection_Areas(
        RCP,
        RP,
        input_file, 
        raster_file,
        processing_file, 
        output_file,
        initial_inland_buffer_distance,
        flood_depth_threshold,
        eps_threshold,
        minPts,
        max_coast_segement_length,
        minGrpSize,
        overlap_threshold,
    ):

    coastline_edge_layer = gpd.read_file(input_file, layer="edges")  # Load edges as GeoDataFrame
    coastline_node_layer = gpd.read_file(input_file, layer="nodes")  # Load nodes as GeoDataFrame

    jamaica_polygon_revised = gpd.read_file(input_file, layer="jam") 
    jamaica_polygon_revised = jamaica_polygon_revised.to_crs(3448)

    jamaica_polygon_buffered = gpd.read_file(input_file, layer="jam_buffered") 
    jamaica_polygon_buffered = jamaica_polygon_revised.to_crs(3448)

    full_coastline_layer = gpd.read_file(input_file, layer="Jamaica_Coastline_Layer")

    #------------------------------------------------------------------------------------------------#
    logging.info ("[1] Creating DBSCAN Cluster")
    # Calling the function with appropriate parameters
    dbscan_flood_areas = process_raster_to_clusters(
        threshold_m=flood_depth_threshold,  # Minimum depth threshold
        eps=eps_threshold,            # DBSCAN epsilon value
        raster_file=raster_file,    # Input raster file path
        minPts = minPts
    )

    # Saving the output to a GeoPackage file
    layer_name = f"dbscan_flood_areas_RP_{RP}_eps_{eps_threshold}_thresh_{flood_depth_threshold}"
    add_Layer_to_File(dbscan_flood_areas, processing_file, layer_name, "GPKG")
    #------------------------------------------------------------------------------------------------#


    #-----------------------------------------SCAPE METHOD-------------------------------------------#
    logging.info ("[2] Generating SCAPE cluster polygons")
    edge_polygons, coastline_edge_layer = generate_polygons_along_edge(
        edge_layer=coastline_edge_layer, 
        node_layer=coastline_node_layer,
        buffer=initial_inland_buffer_distance,
    )

    # Save the output to a file
    layer_name = f"edge_polygons_depth_{initial_inland_buffer_distance}"
    add_Layer_to_File(edge_polygons, processing_file, layer_name, "GPKG")

    layer_name = "edges"
    add_Layer_to_File(coastline_edge_layer, input_file, layer_name, "GPKG")
    
    # Call the function to calculate and add the maximum flood height from the raster
    edge_polygons = add_max_zonal_stats(
        geopkg=edge_polygons,
        raster_path=raster_file,
        new_column_name="max_flood_height"
    )

    # Save the updated GeoDataFrame back to the file, ensuring CRS consistency
    layer_name = f"edge_polygons_depth_{initial_inland_buffer_distance}"
    add_Layer_to_File(edge_polygons, processing_file, layer_name, "GPKG")

    # Call the function to get the ordered list of edge IDs
    ordered_edge_ids = order_edges_around_island(
        edge_layer = coastline_edge_layer
    )

    # Filter the original polygon layer to include only the ordered edge IDs
    ordered_edge_ids = [eid for eid in ordered_edge_ids if eid in edge_polygons["id"].values]
    input_polygons = edge_polygons.set_index("id").loc[ordered_edge_ids].reset_index()

    # Call the function to join polygons based on flood height and area conditions
    joined_polygons = join_polygons_by_flood_height(
        gdf = input_polygons, 
        flood_height_threshold = flood_depth_threshold, 
        length_limit= max_coast_segement_length, 
        group_size = minGrpSize,
        full_coastline_layer=full_coastline_layer
    )

    # print (len(joined_polygons))
    # joined_polygons = gpd.clip(joined_polygons, jamaica_polygon_revised)


    # Define layer name for the final output
    layer_name = f"joined_edge_polygons_flood_{flood_depth_threshold}_area_{max_coast_segement_length}_group_{minGrpSize}"

    # Save the result back to a file, ensuring CRS consistency
    add_Layer_to_File(joined_polygons, processing_file, layer_name, "GPKG")
    
    s_cape_flood_areas = merge_overlapping_polygons(
        gdf = joined_polygons, 
        overlap_threshold = overlap_threshold,
        # max_coast_segement_length = max_coast_segement_length
    )

    s_cape_flood_areas = clip_larger_polygons(s_cape_flood_areas, "id")
    # s_cape_flood_areas = gpd.clip(s_cape_flood_areas, jamaica_polygon_revised)  

    # s_cape_flood_areas["geometry"] = s_cape_flood_areas["geometry"].apply(keep_largest_polygon)

    # Create a layer name dynamically for the merged flood areas
    layer_name = f"scape_flood_areas_flood_{flood_depth_threshold}_area_{max_coast_segement_length}_group_{minGrpSize}_overlap_{overlap_threshold}"

    # Add the merged polygons layer to the file
    add_Layer_to_File(s_cape_flood_areas, processing_file, layer_name, "GPKG")
    #------------------------------------------------------------------------------------------------#


    #------------------------------------------------------------------------------------------------#
    logging.info ("[3] Generating Island Boundary")
    # Generate the convex hull for the flood polygons
    jamaica_convex_hull = gpd.GeoDataFrame({
        'id': [1],  # Assign an ID to the new hull
        'geometry': jamaica_polygon_revised.buffer(10000)  # The geometry of the convex hull
    }, crs=coastline_edge_layer.crs)

    # Define the layer name for the convex hull and save it to the GeoPackage
    layer_name = "jamaica_convex"
    add_Layer_to_File(jamaica_convex_hull, processing_file, layer_name, "GPKG")
    #------------------------------------------------------------------------------------------------#
    

    #------------------------------------------------------------------------------------------------#
    # logging.info ("    -  Generating Voronoi Polygons")
    # Load contour and coastline polygons from the input files
    jam_contour = gpd.read_file(input_file, layer="jam") 

    dbscan_flood_areas_all = process_raster_to_clusters(
        threshold_m=flood_depth_threshold,  # Minimum depth threshold
        eps=eps_threshold,            # DBSCAN epsilon value
        raster_file=raster_file,    # Input raster file path
        minPts=5
    )

    # Combine flood polygons for each Voronoi polygon using the function above
    scape_intersection = combine_floods_within_scape(
        voronoi_polygons = s_cape_flood_areas,  # Voronoi polygons GeoDataFrame
        flood_polygons = dbscan_flood_areas_all  # Flood polygons GeoDataFrame
    )

    # Add the resulting dbscan_voronoi_intersection layer to a processing file (GeoPackage format)
    layer_name = "scape_voronoi_intersection"
    add_Layer_to_File(scape_intersection, processing_file, layer_name, "GPKG")
    #------------------------------------------------------------------------------------------------#

    #------------------------------------------------------------------------------------------------#
    logging.info ("[4] Generating Voronoi Polygons")
    # Create the Voronoi polygons clipped to Jamaica's convex hull
    voronoi_polygons, voronoi_centroid_points = create_voronoi(
        jamaica_polygon = jamaica_convex_hull, 
        jam_contour = jam_contour,
        coastline_polygons = s_cape_flood_areas,
        intersections= scape_intersection,
        length_limit = max_coast_segement_length,
        coastline= full_coastline_layer
    )

    # Save the Voronoi polygons and centroid points to a GeoPackage
    layer_name = "voronoi_polygons"
    add_Layer_to_File(voronoi_polygons, processing_file, layer_name, "GPKG")

    layer_name = "voronoi_centroid_points"
    add_Layer_to_File(voronoi_centroid_points, processing_file, layer_name, "GPKG")

    #------------------------------------------------------------------------------------------------#

    #------------------------------------------------------------------------------------------------#
    logging.info ("[5] Extend Voronoi to enclose flooding")
    # Combine flood polygons for each Voronoi polygon using the function above
    dbscan_voronoi_intersection = combine_floods_within_voronoi(
        voronoi_polygons = voronoi_polygons,  # Voronoi polygons GeoDataFrame
        flood_polygons = dbscan_flood_areas,  # Flood polygons GeoDataFrame
        voronoi_points = voronoi_centroid_points,  # Voronoi centroids GeoDataFrame
        pts_id_field= "c_id",
        max_distance = 80_000_000  # Maximum distance (in km) for filtering based on proximity to Voronoi points
    )

    # Add the resulting dbscan_voronoi_intersection layer to a processing file (GeoPackage format)
    layer_name = "dbscan_voronoi_intersection"
    add_Layer_to_File(dbscan_voronoi_intersection, processing_file, layer_name, "GPKG")
    #------------------------------------------------------------------------------------------------#

    #------------------------------------------------------------------------------------------------#
    # Calculate initial and final intersections
    flood_protection_coastline, flood_protection_areas = linestring_intersect_polygons(
        default_inland_distance=0.4,  # Buffer radius in kilometers
        coast=coastline_edge_layer,  # Coastline GeoDataFrame
        polygons=voronoi_polygons,  # Voronoi polygons GeoDataFrame
        flood_polygons=dbscan_voronoi_intersection,  # Flood polygons GeoDataFrame from previous step
    )

    flood_protection_coastline["geometry"] = flood_protection_coastline["geometry"].apply(lambda geom: linemerge(geom) if geom.geom_type == "MultiLineString" else geom)
    flood_protection_areas = gpd.clip(flood_protection_areas, jamaica_polygon_buffered)

    # Save the raw intersections to the processing file
    # add_Layer_to_File(flood_protection_coastline, processing_file, "flood_protection_coast", "GPKG")

    # Save the final flood protection areas to the processing file
    add_Layer_to_File(flood_protection_areas, processing_file, "flood_protection_areas", "GPKG")
    #------------------------------------------------------------------------------------------------#

    #------------------------------------------------------------------------------------------------#
    logging.info ("[6] Refining Flood Areas - stage 1")
    flood_protection_areas_v2 = split_multipart_polygons(flood_protection_areas)
    add_Layer_to_File(flood_protection_areas_v2, processing_file, "flood_protection_areas_v2", "GPKG")
    #------------------------------------------------------------------------------------------------#

    #------------------------------------------------------------------------------------------------#

    logging.info ("[7] Refining Flood Areas - stage 2")
    flood_protection_areas_v3 = merge_adjacent_polygons_by_coastline(flood_protection_areas_v2, full_coastline_layer, max_coast_segement_length)
    add_Layer_to_File(flood_protection_areas_v3, processing_file, "flood_protection_areas_v3", "GPKG")

    #------------------------------------------------------------------------------------------------#

    #------------------------------------------------------------------------------------------------#
    # Process the cleaned GeoDataFrame by adding the maximum zonal statistics from the raster data
    logging.info ("[8] Adding Flood Height Information")
    final_coastal_protection_area = add_max_zonal_stats(
        geopkg = flood_protection_areas_v3,  # The cleaned GeoDataFrame
        raster_path = raster_file,  # Path to the raster file containing flood height data
        new_column_name= "max_flood_height"  # New column to store the maximum flood height value
    )
    #------------------------------------------------------------------------------------------------#
    
    #------------------------------------------------------------------------------------------------#
    logging.info ("[9] Find Coastline Segments")
    def find_coast_segment_for_polygon(coastline, polygon):
        # if coastline.geometry.iloc[0].intersects(polygon.geometry):
        intersection = coastline.geometry.iloc[0].intersection(polygon.geometry.buffer(0.1))
        return intersection

    flood_protection_coastline = []

    for idx, poly in final_coastal_protection_area.iterrows():
        intersection = find_coast_segment_for_polygon(full_coastline_layer, poly)
        # print (intersection.length)
        flood_protection_coastline.append({
            "geometry": intersection,
            "id": poly["id"],
            "length": intersection.length
        })

    flood_protection_coastline = gpd.GeoDataFrame(flood_protection_coastline, crs=flood_protection_areas_v2.crs)
    flood_protection_coastline = flood_protection_coastline.merge(final_coastal_protection_area[['id', 'max_flood_height']], on='id', how='left')

    add_Layer_to_File(flood_protection_coastline, processing_file, "flood_protection_coast", "GPKG")
    #------------------------------------------------------------------------------------------------#
    
    #------------------------------------------------------------------------------------------------#
    # Define the layer name to be used for storing the final protection area in a GeoPackage
    layer_name = "final_protection_area"
    add_Layer_to_File(final_coastal_protection_area, processing_file, layer_name, "GPKG")

    # Also save the final protection area to the output file in GeoPackage format
    layer_name = f"areas"
    add_Layer_to_File(final_coastal_protection_area, output_file, layer_name, "GPKG")

    final_coastal_protection_coastline = flood_protection_coastline
    layer_name = f"edges"
    add_Layer_to_File(final_coastal_protection_coastline, output_file, layer_name, "GPKG")
    #------------------------------------------------------------------------------------------------#

@click.command()
@click.version_option("1.0")
@click.option(
    "--island-inputs",
    "-i",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the foundation data assets",
)

@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to output gpkg",
)

@click.option(
    "--data-dir",
    "-d",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=False, readable=True),
    help="Path to processed data",
)

@click.option(
    "--inland-buffer",
    "-ib",
    default=2,
    required=True,
    type=float,
    help="Initial inland buffer distance for coastal polygons",
)

@click.option(
    "--flood-threshold",
    "-ft",
    default=0.1,
    required=True,
    type=float,
    help="flood depth threshold",
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
    "--eps",
    "-e",
    default=3,
    required=True,
    type=int,
    help="DBSCAN epsilon value",
)

@click.option(
    "--minpts",
    "-m",
    default=5,
    required=True,
    type=int,
    help="DBSCAN minPts value",
)

@click.option(
    "--coast-length",
    "-cl",
    default=5000,
    required=True,
    type=float,
    help="max coastal segment length for final coastal protection assets",
)

@click.option(
    "--overlap",
    "-ov",
    default=0.5,
    required=True,
    type=float,
    help="Minimum % groups can overlap to trigger merging (Used in SCAPE)",
)

def main(island_inputs, output_dir, data_dir, inland_buffer, flood_threshold, rcp, rp, epoch, eps, minpts, coast_length, overlap):

    #------------------ FILE PATHS -----------------------
    input_file = island_inputs  # GeoPackage file with coastline edges and nodes

    output_dir = Path(output_dir) / "coastal_protection_assets"
    output_dir.mkdir(parents=True, exist_ok=True)

    processing_file = f"{output_dir}/coastal_protection_processing_layers.gpkg"  # GeoPackage file to store intermediate processing layers
    output_file = f"{output_dir}/Jamaica_coastal_protection_areas.gpkg"

    # print (processing_file)
    #------------------------------------------------------

    #------------------ PARAMS -----------------------
    initial_inland_buffer_distance = inland_buffer  # Buffer distance (in kilometers) to expand the inital bounding boxes & coastline buffeering
    flood_depth_threshold = flood_threshold  # Flood height threshold above which flood pixels will be considered
    eps_threshold = eps # DBSCAN epsilon value
    minPts = minpts #DBSCAN minPts value & Minimum number of polygons to make up a group (Used in SCAPE)
    minGrpSize = 2
    max_coast_segement_length = coast_length
    overlap_threshold = overlap #Minimum % groups can overlap to trigger merging (Used in SCAPE)
    # max_separation_distance = 8_000 #Max Distacne flooding can be from coastline to be included in flood area
    
    # print (initial_inland_buffer_distance, flood_depth_threshold, eps_threshold, minPts, max_coast_segement_length, overlap_threshold)

    RP = [f"{rp}"]

    fl_map_rcp = f"{int(rcp*10)}{epoch}"
    RCP = [fl_map_rcp]
    # RCP = ['262100']
    #-------------------------------------------------

    #-------------------------------------------------
    folder = f"{data_dir}/hazards/Coastal_flood_data"
    for rcp in RCP:
        if rcp == 'baseline2010':
            direc = f'{folder}/Flood_maps_current_climate'
        else:
            pt1 = rcp[:2]
            pt2 = rcp[2:]
            direc = f'{folder}/Flood_maps_future_climate/RCP{pt1}_{pt2}'

        for rp in RP:
            raster_file = f"{direc}/JamaicaJAM001RCP{rcp}_epsg_32618_RP_{rp}.tif"
            logging.info (f"    <-      Processing For RCP - {rcp} for RP - {rp}      ->")
            Generate_Coastal_Flood_Protection_Areas(
                RCP = rcp,
                RP = rp,
                input_file = input_file, 
                raster_file = raster_file,
                processing_file = processing_file, 
                output_file = output_file,
                initial_inland_buffer_distance = initial_inland_buffer_distance,
                flood_depth_threshold = flood_depth_threshold,
                eps_threshold = eps_threshold,
                minPts = minPts,
                max_coast_segement_length = max_coast_segement_length,
                minGrpSize = minGrpSize,
                overlap_threshold = overlap_threshold,
                # max_separation_distance = max_separation_distance
            )
    #-------------------------------------------------


if __name__ == "__main__":

    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )

    main()
    