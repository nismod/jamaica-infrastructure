"""Find areas of different layers of land use and NBS value for Jamaica
Find the areas by different classes as totals and percentages
"""

import sys
import os
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

from utils import load_config

tqdm.pandas()

# The projection system for Jamaica
epsg_jamaica = 3448
m2_to_hectare = 0.0001


def get_areas(gdf):
    gdf = gdf.to_crs(epsg=epsg_jamaica)
    gdf["area_m2"] = gdf.progress_apply(lambda x: x.geometry.area, axis=1)
    gdf["area_hectares"] = m2_to_hectare * gdf["area_m2"]
    gdf["percentage"] = 100.0 * gdf["area_m2"] / gdf["area_m2"].sum()

    return gdf


def round_and_convert_to_thousands(gdf, columns_rounding):
    for i, (col_name, value_precision) in enumerate(columns_rounding):
        gdf[col_name] = gdf.apply(
            lambda x: "{:,}".format(round(x[col_name], value_precision)), axis=1
        )
    return gdf


def group_areas_and_find_totals(gdf, groupby_column, columns_rename, columns_rounding):
    gdf = get_areas(gdf)
    gdf = (
        gdf.groupby([groupby_column])["area_m2", "area_hectares", "percentage"]
        .sum()
        .reset_index()
    )

    gdf.columns = columns_rename
    gdf = round_and_convert_to_thousands(gdf, columns_rounding)
    return gdf


def main(config):
    # Specify the data folder paths where all datasets as stored
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]

    nbs_input_data_path = os.path.join(incoming_data_path, "nbs")

    jamaica_land_area = gpd.read_file(
        os.path.join(processed_data_path, "boundaries", "admin_boundaries.gpkg"),
        layer="admin0",
    )
    jamaica_land_area = jamaica_land_area.to_crs(epsg=epsg_jamaica)
    total_land_area = jamaica_land_area.geometry.area.sum()
    print("Total Jamaica admin area in kilometer square", 1.0e-6 * total_land_area)

    forest_land_use = gpd.read_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "input_land_use_layers.gpkg"
        ),
        layer="forest_land_use",
    )
    forest_land_use["classes_lower"] = forest_land_use.progress_apply(
        lambda x: str(x["Classify"]).lower(), axis=1
    )
    forest_land_use = get_areas(forest_land_use)
    total_forest_land_use_area = forest_land_use["area_m2"].sum()
    print(
        "Total Jamaica land use area in kilometer square",
        1.0e-6 * total_forest_land_use_area,
    )

    excel_writer = pd.ExcelWriter(
        os.path.join(processed_data_path, "nbs", "land_use_layers_with_areas.xlsx")
    )

    forest_grouping_column = "Classify"
    forest_land_use_by_class = group_areas_and_find_totals(
        forest_land_use,
        forest_grouping_column,
        [
            "Land-Use Class",
            "Area (m2)",
            "Area (Hectares)",
            "Area (% of total mainland area)",
        ],
        [
            ("Area (m2)", 0),
            ("Area (Hectares)", 0),
            ("Area (% of total mainland area)", 2),
        ],
    )
    print(forest_land_use_by_class)

    forest_land_use_by_class.to_excel(excel_writer, "all_landuse", index=False)

    forest_only_classes = pd.read_excel(
        os.path.join(incoming_data_path, "nbs", "NbS information.xlsx"),
        sheet_name="forest_layers",
    )
    forest_only_classes = [
        str(f).lower() for f in forest_only_classes["Forest sublayers"].values.tolist()
    ]
    print(forest_only_classes)
    print(list(set(forest_land_use["classes_lower"].values.tolist())))
    forest_only_layers = forest_land_use[
        forest_land_use["classes_lower"].isin(forest_only_classes)
    ]
    print(forest_only_layers)

    forest_only_use = group_areas_and_find_totals(
        forest_only_layers,
        forest_grouping_column,
        [
            "Land-Use Class",
            "Area (m2)",
            "Area (Hectares)",
            "Area (% of total forest area)",
        ],
        [
            ("Area (m2)", 0),
            ("Area (Hectares)", 0),
            ("Area (% of total forest area)", 2),
        ],
    )
    forest_only_use.to_excel(excel_writer, "forest_layers", index=False)

    agri_only_classes = pd.read_excel(
        os.path.join(incoming_data_path, "nbs", "NbS information.xlsx"),
        sheet_name="agriculture_forest_layers",
    )
    agri_only_classes = [
        str(a).lower()
        for a in agri_only_classes["Agriculture sublayers"].values.tolist()
    ]
    agri_only_layers = forest_land_use[
        forest_land_use["classes_lower"].isin(agri_only_classes)
    ]
    agri_only_use = group_areas_and_find_totals(
        agri_only_layers,
        forest_grouping_column,
        [
            "Land-Use Class",
            "Area (m2)",
            "Area (Hectares)",
            "Area (% of total agriculture area)",
        ],
        [
            ("Area (m2)", 0),
            ("Area (Hectares)", 0),
            ("Area (% of total agriculture area)", 2),
        ],
    )
    print(agri_only_use)
    agri_only_use.to_excel(excel_writer, "agriculture_layers", index=False)

    natural_layers = [
        "nsmdb-mangroves.gpkg",
        "nsmdb-seagrass.gpkg",
        "nsmdb-reefs_jan06_NEPA.gpkg",
    ]
    layer_name = ["mangroves", "seagrass", "corals"]
    areas_and_percentages = []
    for i, (input_layer, output_layer) in enumerate(
        list(zip(natural_layers, layer_name))
    ):
        get_layer = gpd.read_file(os.path.join(nbs_input_data_path, input_layer))
        get_layer = get_areas(get_layer)
        total_ = get_layer["area_m2"].sum()
        areas_and_percentages.append(
            (
                output_layer,
                get_layer["area_m2"].sum(),
                get_layer["area_hectares"].sum(),
                100.0 * get_layer["area_m2"].sum() / total_forest_land_use_area,
            )
        )

    areas_and_percentages = pd.DataFrame(
        areas_and_percentages,
        columns=[
            "Coastal Nbs",
            "Area (m2)",
            "Area (Hectares)",
            "Area (% of total mainland area)",
        ],
    )
    areas_and_percentages = round_and_convert_to_thousands(
        areas_and_percentages,
        [
            ("Area (m2)", 0),
            ("Area (Hectares)", 0),
            ("Area (% of total mainland area)", 2),
        ],
    )
    print(areas_and_percentages)
    areas_and_percentages.to_excel(excel_writer, "nbs_coastal_layers", index=False)

    excel_writer.save()


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
