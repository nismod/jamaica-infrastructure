"""Create buffers around NBS = mangrove, seagrass and coral layers for Jamaica
    Find the intersections of different buffer boundaries
"""

import sys
import os
import pandas as pd
import geopandas as gpd
import numpy as np
from utils import *
from tqdm import tqdm

tqdm.pandas()
# The projection system for Jamaica
epsg_jamaica = 3448


def intersect_layers_and_find_area(layer_i, layer_j, layer_j_id, all_ids):
    common_ij = gpd.sjoin(
        layer_i, layer_j, how="inner", predicate="intersects"
    ).reset_index()
    common_ij.rename(columns={"geometry": "layer_i_geometry"}, inplace=True)
    common_ij = pd.merge(common_ij, layer_j, how="left", on=[layer_j_id])
    common_ij.rename(columns={"geometry": "layer_j_geometry"}, inplace=True)
    common_ij["geometry"] = common_ij.progress_apply(
        lambda x: x.layer_i_geometry.intersection(x.layer_j_geometry), axis=1
    )
    common_ij["area_m2"] = common_ij.progress_apply(lambda x: x.geometry.area, axis=1)
    common_ij["area_hectares"] = 0.0001 * common_ij["area_m2"]
    common_ij.drop(["layer_i_geometry", "layer_j_geometry"], axis=1, inplace=True)
    common_ij = gpd.GeoDataFrame(
        common_ij[all_ids + ["area_m2", "geometry"]],
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    return common_ij


def main(config):
    # Specify the data folder paths where all datasets as stored
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]

    nbs_input_data_path = os.path.join(incoming_data_path, "nbs")

    # Create a for loop over the buffer distances and the layer
    # Add the buffer to the geometry of the layer
    # Write the result to a new geopackage file
    natural_layers = [
        "nsmdb-mangroves.gpkg",
        "nsmdb-seagrass.gpkg",
        "nsmdb-reefs_jan06_NEPA.gpkg",
    ]
    output_layer_name = ["mangroves", "seagrass", "corals"]
    for b in buffer_distances:
        for i, (input_layer, output_layer) in enumerate(
            list(zip(natural_layers, output_layer_name))
        ):
            # Read a Geopackage layer or shapefile
            get_layer = gpd.read_file(os.path.join(nbs_input_data_path, input_layer))
            # Reproject to Jamaica projection system
            get_layer = get_layer.to_crs(epsg=epsg_jamaica)
            # Create buffer around geometry
            get_layer["geometry"] = get_layer.progress_apply(
                lambda x: x.geometry.buffer(b), axis=1
            )
            get_layer.to_file(
                os.path.join(processed_data_path, "nbs", f"{output_layer}_buffer.gpkg"),
                layer=f"{b}m",
                driver="GPKG",
            )

    # Intersect the different buffer layers in pairs to find common areas
    # Buffer distances in meters
    buffer_distances = [250, 500, 1000]
    intersection_layers = [
        {
            "layer": "mangroves",
            "id": "ID",
        },
        {"layer": "seagrass", "id": "ID"},
        {"layer": "corals", "id": "OBJECTID"},
    ]
    for b in buffer_distances:
        for i in range(len(intersection_layers) - 1):
            name_i = intersection_layers[i]["layer"]
            id_i = intersection_layers[i]["id"]
            for j in intersection_layers[i + 1 :]:
                name_j = j["layer"]
                id_j = j["id"]
                layer_i = gpd.read_file(
                    os.path.join(processed_data_path, "nbs", f"{name_i}_buffer.gpkg"),
                    layer=f"{b}m",
                )
                layer_i.rename(columns={id_i: f"{name_i}_id"}, inplace=True)
                print(layer_i)
                layer_j = gpd.read_file(
                    os.path.join(processed_data_path, "nbs", f"{name_j}_buffer.gpkg"),
                    layer=f"{b}m",
                )
                layer_j.rename(columns={id_j: f"{name_j}_id"}, inplace=True)
                print(layer_j)
                common_ij = intersect_layers_and_find_area(
                    layer_i[[f"{name_i}_id", "geometry"]],
                    layer_j[[f"{name_j}_id", "geometry"]],
                    f"{name_j}_id",
                    [f"{name_i}_id", f"{name_j}_id"],
                )

                print(common_ij)
                common_ij.to_file(
                    os.path.join(
                        processed_data_path,
                        "nbs",
                        f"{name_i}_{name_j}_intersections.gpkg",
                    ),
                    layer=f"{b}m",
                    driver="GPKG",
                )

    # Intersect the different buffer layers together to find common areas
    # Buffer distances in meters
    buffer_distances = [250, 500, 1000]
    for b in buffer_distances:
        mangrove_layer = gpd.read_file(
            os.path.join(processed_data_path, "nbs", f"mangroves_buffer.gpkg"),
            layer=f"{b}m",
        )
        mangrove_layer.rename(columns={"ID": "mangroves_id"}, inplace=True)
        seagrass_layer = gpd.read_file(
            os.path.join(processed_data_path, "nbs", f"seagrass_buffer.gpkg"),
            layer=f"{b}m",
        )
        seagrass_layer.rename(columns={"ID": "seagrass_id"}, inplace=True)
        corals_layer = gpd.read_file(
            os.path.join(processed_data_path, "nbs", f"corals_buffer.gpkg"),
            layer=f"{b}m",
        )
        corals_layer.rename(columns={"OBJECTID": "corals_id"}, inplace=True)
        mangroves_seagrass = intersect_layers_and_find_area(
            mangrove_layer[["mangroves_id", "geometry"]],
            seagrass_layer[["seagrass_id", "geometry"]],
            "seagrass_id",
            ["mangroves_id", "seagrass_id"],
        )

        mangroves_seagrass_corals = intersect_layers_and_find_area(
            mangroves_seagrass[["mangroves_id", "seagrass_id", "geometry"]],
            corals_layer[["corals_id", "geometry"]],
            "corals_id",
            ["mangroves_id", "seagrass_id", "corals_id"],
        )
        mangroves_seagrass_corals.to_file(
            os.path.join(
                processed_data_path,
                "nbs",
                f"mangroves_seagrass_corals_intersections.gpkg",
            ),
            layer=f"{b}m",
            driver="GPKG",
        )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
