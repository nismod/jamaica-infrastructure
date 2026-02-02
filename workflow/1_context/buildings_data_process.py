"""Take the buildings footprints from OSM and add attributes to them economic sectors
    Write the final buildings footprints into a geopackage
"""

import sys
import os

import geopandas as gpd
import pandas as pd
import numpy as np
import fiona
from collections import OrderedDict, defaultdict
from shapely.geometry import Point
from preprocess_utils import *
from tqdm import tqdm

tqdm.pandas()


def quantiles(dataframe, grouping_by_columns, grouped_columns):
    quantiles_list = ["mean", "min", "max", "median", "q5", "q95"]
    df_list = []
    for quant in quantiles_list:
        if quant == "mean":
            # print (dataframe)
            df = dataframe.groupby(grouping_by_columns, dropna=False)[
                grouped_columns
            ].mean()
        elif quant == "min":
            df = dataframe.groupby(grouping_by_columns, dropna=False)[
                grouped_columns
            ].min()
        elif quant == "max":
            df = dataframe.groupby(grouping_by_columns, dropna=False)[
                grouped_columns
            ].max()
        elif quant == "median":
            df = dataframe.groupby(grouping_by_columns, dropna=False)[
                grouped_columns
            ].quantile(0.5)
        elif quant == "q5":
            df = dataframe.groupby(grouping_by_columns, dropna=False)[
                grouped_columns
            ].quantile(0.05)
        elif quant == "q95":
            df = dataframe.groupby(grouping_by_columns, dropna=False)[
                grouped_columns
            ].quantile(0.95)

        df.rename(
            columns=dict((g, "{}_{}".format(g, quant)) for g in grouped_columns),
            inplace=True,
        )
        df_list.append(df)
    return pd.concat(df_list, axis=1).reset_index()


def find_polygon_dimensions(x):
    # create example polygon
    poly = x.geometry

    # get minimum bounding box around polygon
    box = poly.minimum_rotated_rectangle

    # get coordinates of polygon vertices
    x, y = box.exterior.coords.xy

    # get length of bounding box edges
    edge_length = (
        Point(x[0], y[0]).distance(Point(x[1], y[1])),
        Point(x[1], y[1]).distance(Point(x[2], y[2])),
    )

    # get length of polygon as the longest edge of the bounding box
    length = max(edge_length)

    # get width of polygon as the shortest edge of the bounding box
    width = min(edge_length)

    return length, width, length / width


def get_sector_subsector(x, column_pairs):
    """This function filters out the most likely sector and subsector codes assigned to a building
    In an OSM dataset, based on a mapping of building attributes to macroeconomic sectors

    Based on our mapping a building might be either assigned a sector based on some attributes
        Or not assigned any sector, and instead assigned a code 'X'
    If a building is assigned a code 'X' in addition to a known macroeconomic sector code
        Then the code 'X' is removed and the sector code is retained

    The output of the function results in 1 code assigned to a building
    """
    p0 = []
    p1 = []
    for cl in column_pairs:
        p0 += x[cl[0]].split(",")
        p1 += x[cl[1]].split(",")

    vals = list(OrderedDict.fromkeys(zip(p0, p1)))
    # vals = list(set(zip(p0,p1)))
    if len(vals) > 1:
        sector_val = ",".join([v[0] for v in vals if v[0] != "X"])
        subsector_val = ",".join([v[1] for v in vals if v[1] != "X"])
    else:
        sector_val = vals[0][0]
        subsector_val = vals[0][1]

    return sector_val, subsector_val


def join_values(x, column_types):
    """This function filters out the most likely sector and subsector codes assigned to a building
    In an OSM dataset, based on a mapping of building attributes to macroeconomic sectors

    Based on our mapping a building might be either assigned a sector based on some attributes
        Or not assigned any sector, and instead assigned a code 'X'
    If a building is assigned a code 'X' in addition to a known macroeconomic sector code
        Then the code 'X' is removed and the sector code is retained

    The output of the function results in 1 code assigned to a building
    """
    p0 = []
    for cl in column_types:
        p0 += str(x[cl]).split(",")

    vals = list(OrderedDict.fromkeys(p0))
    # vals = list(set(p0))
    # print (p0,vals)
    if len(vals) > 1:
        return ",".join(
            [
                v
                for v in vals
                if v not in ["", "X", "nan", "None", "none", "yes", "Yes", "Unknown"]
            ]
        )
    else:
        return vals[0]


def join_values_in_column(x, column_merge):
    """This function filters out the most likely sector and subsector codes assigned to a building
    In an OSM dataset, based on a mapping of building attributes to macroeconomic sectors

    Based on our mapping a building might be either assigned a sector based on some attributes
        Or not assigned any sector, and instead assigned a code 'X'
    If a building is assigned a code 'X' in addition to a known macroeconomic sector code
        Then the code 'X' is removed and the sector code is retained

    The output of the function results in 1 code assigned to a building
    """
    p0 = str(x[column_merge]).split(",")
    vals = list(OrderedDict.fromkeys(p0))

    if len(vals) > 1:
        return ",".join(
            [v for v in vals if v not in ["", "X", "nan", "None", "none", "yes", "Yes"]]
        )
    else:
        return vals[0]


def create_sector_subsector_attribute_columns(
    building_dataframe,
    building_id_column,
    merge_columns,
    new_columns,
    attribute_columns,
    infra_columns,
):
    for col in new_columns + attribute_columns + infra_columns:
        building_dataframe[col] = building_dataframe[col].apply(str)
        building_dataframe[col] = building_dataframe[col].replace(
            [r"^\s+$", "", "nan", "None", "none"], "X", regex=True
        )

    building_dataframe["sector_subsector"] = building_dataframe.progress_apply(
        lambda x: get_sector_subsector(x, merge_columns), axis=1
    )
    building_dataframe[["sector_code", "subsector_code"]] = building_dataframe[
        "sector_subsector"
    ].apply(pd.Series)
    attribute_columns = [c for c in attribute_columns if c != "other_tags"]
    building_dataframe["assigned_attribute"] = building_dataframe.progress_apply(
        lambda x: join_values(x, attribute_columns), axis=1
    )
    building_dataframe["infra_type"] = building_dataframe.progress_apply(
        lambda x: join_values(x, infra_columns), axis=1
    )
    building_dataframe.drop(
        new_columns + attribute_columns + infra_columns + ["sector_subsector"],
        axis=1,
        inplace=True,
    )
    building_dataframe = building_dataframe.drop_duplicates(
        subset=[building_id_column], keep="first"
    )

    return building_dataframe


def sector_subsector_asset_mapping(
    building_dataframe,
    building_id_column,
    building_sector_mapping_file,
    attribute_columns,
):
    merge_columns = []
    new_columns = []
    infra_columns = []
    for column in attribute_columns:
        sector_map = pd.read_excel(building_sector_mapping_file, sheet_name=column)
        sector_map = sector_map.drop_duplicates(
            [column, "sector_code", "subsector_code", "infra_sector"], keep="first"
        )
        sector_map = sector_map[
            [column, "sector_code", "subsector_code", "infra_sector"]
        ]
        sector_map.rename(
            columns={
                "sector_code": f"{column}_sector_code",
                "subsector_code": f"{column}_subsector_code",
                "infra_sector": f"{column}_infra_sector",
            },
            inplace=True,
        )
        merge_columns += [(f"{column}_sector_code", f"{column}_subsector_code")]
        new_columns += [f"{column}_sector_code", f"{column}_subsector_code"]
        infra_columns += [f"{column}_infra_sector"]
        building_dataframe = pd.merge(
            building_dataframe, sector_map, how="left", on=[column]
        )

    return create_sector_subsector_attribute_columns(
        building_dataframe,
        building_id_column,
        merge_columns,
        new_columns,
        attribute_columns,
        infra_columns,
    )


def match_points_to_nearest_buildings(points_of_interest, points_file, sector_codes):
    """Match points to their nearest buildings corresponding to a set of sector codes"""
    buildings_select = points_of_interest[
        points_of_interest["sector_code"].isin(sector_codes)
    ].reset_index()
    sindex_buildings = buildings_select.sindex
    points_file["osm_id"] = points_file.geometry.progress_apply(
        lambda x: get_nearest_node(x, sindex_buildings, buildings_select, "osm_id")
    )
    points_file["geom"] = points_file.geometry.progress_apply(
        lambda x: get_nearest_node(x, sindex_buildings, buildings_select, "geometry")
    )
    points_file["distance"] = points_file.progress_apply(
        lambda x: x.geometry.distance(x.geom), axis=1
    )
    points_file.drop(["geom"], axis=1, inplace=True)

    return points_file


def identify_specific_areas(x, area_string):
    if (
        area_string in str(x["Classify"]).lower()
        or area_string in str(x["NAME"]).lower()
        or area_string in str(x["global_LU_type"]).lower()
    ):
        return 1
    else:
        return 0


def column_merging(
    value_dataframe, building_dataframe, building_id_column, dataframe_merge_columns
):
    for column in ["sector_code", "subsector_code"]:
        value_dataframe.rename(columns={column: f"{column}_match"}, inplace=True)
        building_dataframe.rename(columns={column: f"{column}_building"}, inplace=True)

    merge_columns = []
    new_columns = []
    # columns_types = []
    for column in ["match", "building"]:
        merge_columns += [(f"sector_code_{column}", f"subsector_code_{column}")]
        new_columns += [f"sector_code_{column}", f"subsector_code_{column}"]
        # columns_types.append(f'assigned_attribute_{column}')

    building_dataframe = pd.merge(
        building_dataframe, value_dataframe, how="left", on=dataframe_merge_columns
    )
    building_dataframe[new_columns] = building_dataframe[new_columns].astype("str")
    building_dataframe[new_columns] = building_dataframe[new_columns].replace(
        [r"^\s+$", "nan", "None", "none"], "X", regex=True
    )
    building_dataframe["sector_subsector"] = building_dataframe.progress_apply(
        lambda x: get_sector_subsector(x, merge_columns), axis=1
    )
    building_dataframe[["sector_code", "subsector_code"]] = building_dataframe[
        "sector_subsector"
    ].apply(pd.Series)
    building_dataframe.drop(new_columns + ["sector_subsector"], axis=1, inplace=True)
    building_dataframe = building_dataframe.drop_duplicates(
        subset=[building_id_column], keep="first"
    )

    return building_dataframe


def match_buildings_to_landplanning(buildings, gdf, gdf_list):
    matches = gpd.sjoin(
        buildings[["osm_id", "geometry"]],
        gdf[["layer_id", "geometry"]],
        how="inner",
        predicate="intersects",
    ).reset_index()
    # print (matches)
    matches.rename(columns={"geometry": "building_geometry"}, inplace=True)
    matches = pd.merge(
        matches, gdf[["layer_id", "geometry"]], how="left", on=["layer_id"]
    )
    matches["area_match"] = matches.progress_apply(
        lambda x: x["building_geometry"].intersection(x["geometry"].buffer(0)).area,
        axis=1,
    )
    matches = matches.sort_values(by=["area_match"], ascending=False)
    matches = matches.drop_duplicates(subset=["layer_id"], keep="first")

    gdf = pd.merge(gdf, matches[["layer_id", "PARISH"]], how="left", on=["layer_id"])
    parish = gdf["PARISH"].mode().values[0]
    gdf["PARISH"] = gdf["PARISH"].astype(str)
    gdf["PARISH"] = gdf["PARISH"].replace(
        [r"^\s+$", "nan", "None", "none"], parish, regex=True
    )
    gdf.drop("layer_id", axis=1, inplace=True)
    gdf_list.append(gdf)

    return gdf_list


def match_buildings_to_polygon_dataframe(buildings, gdf, epsg=4326, spatial_join=True):
    if spatial_join is True:
        # Remove some columns from the polygon layer
        if "osm_id" in gdf.columns.values.tolist():
            gdf.drop("osm_id", axis=1, inplace=True)
        gdf = gdf.to_crs(epsg=epsg)
        buildings = buildings.to_crs(epsg=epsg)
        gdf_intersections = gpd.sjoin(
            gdf, buildings[["osm_id", "geometry"]], how="inner", predicate="intersects"
        ).reset_index()
        df = pd.DataFrame(
            list(set(gdf_intersections["osm_id"].values.tolist())), columns=["osm_id"]
        )

        for column in [
            "sector_code",
            "subsector_code",
            "assigned_attribute",
            "infra_type",
        ]:
            gdf_intersections[column] = gdf_intersections[column].apply(str)
            gdf_intersections[column] = gdf_intersections[column].replace(
                [r"^\s+$", "", "nan", "None", "none"], "X", regex=True
            )

            df = pd.merge(
                df,
                gdf_intersections[
                    [
                        "osm_id",
                        column,
                    ]
                ]
                .groupby(["osm_id"])[column]
                .apply(",".join)
                .reset_index(),
                how="left",
                on=["osm_id"],
            )
            df.rename(columns={column: f"{column}_match"}, inplace=True)
            buildings.rename(columns={column: f"{column}_building"}, inplace=True)
    else:
        df = gdf.copy()

    merge_columns = []
    new_columns = []
    columns_types = []
    infra_columns = []
    for column in ["match", "building"]:
        merge_columns += [(f"sector_code_{column}", f"subsector_code_{column}")]
        new_columns += [f"sector_code_{column}", f"subsector_code_{column}"]
        columns_types.append(f"assigned_attribute_{column}")
        infra_columns.append(f"infra_type_{column}")

    buildings = pd.merge(buildings, df, how="left", on=["osm_id"])
    buildings = create_sector_subsector_attribute_columns(
        buildings, "osm_id", merge_columns, new_columns, columns_types, infra_columns
    )

    return buildings


def possible_nonbuilding(x, land_use_assigned):
    non_building_score = 0
    non_building_type = ""
    b = str(x["assigned_attribute"]).split(",")
    for l in land_use_assigned:
        non_building_type += ",".join(
            [a for a in b if str(a).lower() == str(l).lower()]
        )
        non_building_score += len([a for a in b if str(a).lower() == str(l).lower()])

    return non_building_score, non_building_type


def assign_building_type(x):
    if x.sector_code == "RES":
        return "Residential"
    elif x.sector_code == "H":
        return "Resort"
    elif x.sector_code == "G":
        return "Commercial"
    elif x.sector_code in ("A", "B", "C", "D"):
        return "Industrial"
    elif x.sector_code in ("L", "M", "N"):
        return "Institutional"
    elif x.sector_code == "O":
        return "Recreation"
    elif len(list(set(str(x.sector_code).split(",")))) >= 2:
        return "Mixed Use"
    else:
        return "Other"


def main(config):
    # Set the file paths on your machine for the data
    # All input data is stored in a path ../incoming_data/..
    # All output data will be written in a path ../processed_data/..

    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    building_data_path = os.path.join(incoming_data_path, "buildings")
    points_of_interest_data_path = os.path.join(
        incoming_data_path, "hotosm_data", "hotosm_jam_points_of_interest_gpkg"
    )
    buildings_input_mapping = os.path.join(
        building_data_path, "osm_building_column_mapping_economic_sectors.xlsx"
    )

    epsg_jamaica = 3448  # Jamaica projection system
    residential_min_area = (
        11.15  # Minimum area of residential dwelling allowed in Jamaica
    )

    """Step 1: Extract the different types of entries in useful columns 
        To infer what kind of building attributes are listed in the OSM building layer
        Using these attributes we want to map a macroeconomic sector to a building based on this attribute 
        The columns in the layer which have useful information are - 'building','amenity','office','shop','other_tags'
    """
    matched_buildings = []
    buildings_input = gpd.read_file(
        os.path.join(building_data_path, "jamaica-buildings.gpkg")
    )
    buildings_input = buildings_input.to_crs(epsg=epsg_jamaica)
    buildings_input[["osm_id", "osm_way_id"]] = buildings_input[
        ["osm_id", "osm_way_id"]
    ].astype("str")
    buildings_input[["osm_id", "osm_way_id"]] = buildings_input[
        ["osm_id", "osm_way_id"]
    ].replace([r"^\s+$", "nan", "None", "none"], "", regex=True)
    buildings_input["osm_id"] = buildings_input["osm_id"].astype(
        "str"
    ) + buildings_input["osm_way_id"].astype("str")

    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        buildings_input.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="original",
            driver="GPKG",
        )

    points_of_interest = gpd.read_file(
        os.path.join(points_of_interest_data_path, "hotosm_jam_points_of_interest.gpkg")
    )
    points_of_interest = points_of_interest.to_crs(epsg=epsg_jamaica)
    points_of_interest["shop"] = points_of_interest["shop"].replace(
        "yes", "wholesale", regex=True
    )

    # """
    # Excel outputs generated for post-processing and
    # Creating a mapping of known building tags in the datasets with economic sectors
    # The rest of the code uses the outputs files from these tags

    # # Extract categories within columns

    # output_excel = os.path.join(incoming_data_path,'buildings','unique_column_mapping.xlsx')
    # output_wrtr = pd.ExcelWriter(output_excel)
    # for column in columns_types:
    #     df = buildings_input[column].value_counts().reset_index()
    #     df.columns = [column,'count']
    #     df.to_excel(output_wrtr,sheet_name=column, index=False)

    # output_wrtr.save()

    # # Extract categories within columns

    # columns_types = ['building','amenity','office','shop']
    # sector_maps = []
    # for column in columns_types:
    #     sector_map = pd.read_excel(os.path.join(
    #                                 building_data_path,
    #                                 'osm_building_column_mapping_economic_sectors.xlsx'),
    #                                 sheet_name=column)
    #     sector_map = sector_map[[column,'sector_code','subsector_code']]
    #     sector_map.rename(columns = {column:'asset_type'},inplace=True)
    #     sector_maps.append(sector_map)

    # sector_maps = pd.concat(sector_maps,axis=0,ignore_index=True)

    # columns_types = ['amenity','man_made','shop','tourism']
    # output_excel = os.path.join(points_of_interest_data_path,
    #                                 'unique_column_mapping.xlsx')
    # output_wrtr = pd.ExcelWriter(output_excel)
    # for column in columns_types:
    #     df = points_of_interest[column].value_counts().reset_index()
    #     df.columns = [column,'count']
    #     df = pd.merge(df,sector_maps,how='left',left_on=[column],right_on='asset_type')
    #     df.drop('asset_type',axis=1,inplace=True)
    #     df.to_excel(output_wrtr,sheet_name=column, index=False)

    # output_wrtr.save()
    # """

    """Step 2: Use a preprocessed file that maps buildings to economic sectors and subsectors
       The preprocessed file builds on the file output of Step 1 where now we have added two new columns
       sector_code & subsector_code - which are the codes for the macroeconomic sectors and subsectors respectively   
       The output of this step is to map the codes onto the spatial building layer to create two new columns in the layer 

       The default values for sector_code and subsector_code is 'X' - which signifies we cannot assign an economic attributes to the building 

       NOTE: From this step onwards we filter out buildings that have been matched to a given land use sector != X
    """
    points_of_interest_mapping = os.path.join(
        points_of_interest_data_path, "hotosm_mapping_economic_sectors.xlsx"
    )
    points_of_interest = sector_subsector_asset_mapping(
        points_of_interest,
        "osm_id",
        points_of_interest_mapping,
        ["amenity", "man_made", "shop", "tourism"],
    )
    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        points_of_interest.to_file(
            os.path.join(
                points_of_interest_data_path,
                "points_of_interest_assigned_economic_sectors_intermediate.gpkg",
            ),
            driver="GPKG",
        )
    buildings_input = sector_subsector_asset_mapping(
        buildings_input,
        "osm_id",
        buildings_input_mapping,
        ["building", "amenity", "office", "shop", "other_tags"],
    )
    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        buildings_input.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="known_attributes",
            driver="GPKG",
        )

    matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    buildings_input = match_buildings_to_polygon_dataframe(
        buildings_input, points_of_interest, epsg=epsg_jamaica
    )
    matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    buildings_input = buildings_input[buildings_input["sector_code"] == "X"]
    write_output = True
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        write_df = gpd.GeoDataFrame(
            pd.concat(matched_buildings + [buildings_input], axis=0, ignore_index=True),
            geometry="geometry",
            crs=f"EPSG:{epsg_jamaica}",
        )
        write_df.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="points_merged",
            driver="GPKG",
        )
        del write_df

    """Step 3: Identify buildings by the planning layers data, which gives further classification of land use types
        We mapped these land use types to economic sectors and subsectors as well 
    """
    landuse_layers = [
        "existing_proposed_landuse",
        "existing_proposed_landuse",
        "existing_landuse",
        "clarendon_landuse",
        "manchester_landuse",
    ]
    landuse_types = ["existing", "proposals", "existing", "existing", "existing"]

    for j, (layer, layer_type) in enumerate(list(zip(landuse_layers, landuse_types))):
        planning_layer = gpd.read_file(
            os.path.join(
                incoming_data_path,
                "buildings",
                "landuse_planning_layers_with_sectors.gpkg",
            ),
            layer=layer,
        )
        if layer_type == "proposals":
            planning_layer.rename(
                columns={
                    "Proposals": "assigned_attribute",
                    "sector_code_proposals": "sector_code",
                    "subsector_code_proposals": "subsector_code",
                    "infra_sector_proposals": "infra_type",
                },
                inplace=True,
            )
        else:
            planning_layer.rename(
                columns={
                    "LU_Zone": "assigned_attribute",
                    "sector_code_existing": "sector_code",
                    "subsector_code_existing": "subsector_code",
                    "infra_sector_existing": "infra_type",
                    "infra_sector": "infra_type",
                },
                inplace=True,
            )
        buildings_input = match_buildings_to_polygon_dataframe(
            buildings_input, planning_layer, epsg=epsg_jamaica
        )
        matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
        buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        write_df = gpd.GeoDataFrame(
            pd.concat(matched_buildings + [buildings_input], axis=0, ignore_index=True),
            geometry="geometry",
            crs=f"EPSG:{epsg_jamaica}",
        )
        write_df.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="land_planning_merged",
            driver="GPKG",
        )
        del write_df

    """Step 4: Identify the buildings within mining and quarrying areas
    """

    land_use_types = gpd.read_file(
        os.path.join(
            processed_data_path, "land_type_and_use", "jamaica_land_use_combined.gpkg"
        ),
        layer="areas",
    )

    if "index_left" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_left", axis=1, inplace=True)
    elif "index_right" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_right", axis=1, inplace=True)

    mining_types = [
        ("bauxite", "C", "132", "Bauxite mining and alumina"),
        ("quarry", "C", "141", "Quarry"),
    ]
    for i, (mtype, msect, msubsect, mname) in enumerate(mining_types):
        land_use_types[f"{mtype}_identify"] = land_use_types.progress_apply(
            lambda x: identify_specific_areas(x, mtype), axis=1
        )
        mining_quarry_areas = land_use_types[land_use_types[f"{mtype}_identify"] == 1]

        mining_quarry_areas["sector_code"] = msect
        mining_quarry_areas["subsector_code"] = msubsect
        mining_quarry_areas["assigned_attribute"] = mname
        mining_quarry_areas["infra_type"] = "X"

        buildings_input = match_buildings_to_polygon_dataframe(
            buildings_input, mining_quarry_areas, epsg=epsg_jamaica
        )
        matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
        buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        write_df = gpd.GeoDataFrame(
            pd.concat(matched_buildings + [buildings_input], axis=0, ignore_index=True),
            geometry="geometry",
            crs=f"EPSG:{epsg_jamaica}",
        )
        write_df.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="mines_added",
            driver="GPKG",
        )
        del write_df

    """Step 5: Add the commercial buildings mapping
    """
    # matched_buildings = []
    # buildings_input = gpd.read_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'),
    #                                 layer="mines_added")
    # matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    # buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    commercial_buildings = gpd.read_file(
        os.path.join(
            incoming_data_path,
            "JAM_classified_buildings_industries",
            "JAM_classified_buildings.shp",
        )
    )[
        [
            "osm_id",
            "layer",
            "Name_of_Pl",
            "SECTOR1",
            "SUB_SECTOR",
            "PRODUCTS",
            "HS_CODE",
            "COMPANY_NA",
            "geometry",
        ]
    ]
    commercial_buildings.rename(
        columns={"SECTOR1": "SECTOR", "Name_of_Pl": "Name_of_Plant"}, inplace=True
    )
    commercial_buildings["layer"] = commercial_buildings["layer"].astype("str")
    commercial_buildings["layer"] = commercial_buildings["layer"].replace(
        [r"^\s+$", "nan", "None", "none", "points_of_interest", " "], "", regex=True
    )
    commercial_buildings["layer"] = commercial_buildings["layer"].replace(
        "quarries", "", regex=True
    )
    commercial_buildings["sector_code"] = "X"
    commercial_buildings["subsector_code"] = "X"

    industry_export_trade_buildings = commercial_buildings[
        commercial_buildings["layer"].isin(["trade_plants", "exporters", "industry"])
    ]
    industry_export_trade_ids = list(
        set(industry_export_trade_buildings["osm_id"].values.tolist())
    )
    social_buildings = commercial_buildings[
        ~commercial_buildings["osm_id"].isin(industry_export_trade_ids)
    ]

    commercial_layers = []
    """Step 6: Do similar point layer mapping for three specific layers of 
        Trade plants, Exporters, and Industry layers
    """

    layer_dict = [
        {
            "layer_name": "trade_plants",
            "excel_file": "nsdmb_tradeplants_layer_economic_sectors.xlsx",
            "layer_columns": ["Name_of_Plant"],
        },
        {
            "layer_name": "exporters",
            "excel_file": "nsdmb_exporters_layer_economic_sectors.xlsx",
            "layer_columns": ["SECTOR", "SUB_SECTOR", "PRODUCTS", "HS_CODE"],
        },
        {
            "layer_name": "industry",
            "excel_file": "nsdmb_industry_layer_economic_sectors.xlsx",
            "excel_sheets": ["manufacturing", "mining"],
            "layer_columns": [["SECTOR"], ["SECTOR", "COMPANY_NA"]],
        },
    ]

    for layer_properties in layer_dict:
        commercial_layer = industry_export_trade_buildings[
            industry_export_trade_buildings["layer"] == layer_properties["layer_name"]
        ]
        if layer_properties["layer_name"] == "industry":
            # points_file['SECTOR'] = points_file['SECTOR'].replace(' ', 'Unknown')
            for sh in range(len(layer_properties["excel_sheets"])):
                nsdmb_layers = pd.read_excel(
                    os.path.join(building_data_path, layer_properties["excel_file"]),
                    sheet_name=layer_properties["excel_sheets"][sh],
                )
                commercial_layer = column_merging(
                    nsdmb_layers,
                    commercial_layer,
                    "osm_id",
                    layer_properties["layer_columns"][sh],
                )
        else:
            nsdmb_layers = pd.read_excel(
                os.path.join(building_data_path, layer_properties["excel_file"]),
                sheet_name="Sheet1",
            )
            commercial_layer = column_merging(
                nsdmb_layers,
                commercial_layer,
                "osm_id",
                layer_properties["layer_columns"],
            )

        commercial_layers.append(commercial_layer)
        print(f"* Done with points layer - {layer_properties['layer_name']}")

    nsdmb_layers = pd.read_excel(
        os.path.join(building_data_path, "nsdmb_layers_economic_sectors.xlsx"),
        sheet_name="Sheet1",
    )
    for layer_properties in nsdmb_layers.itertuples():
        layer_name = layer_properties.layer
        sector_code = layer_properties.sector_code
        subsector_code = layer_properties.subsector_code

        commercial_layer = social_buildings[social_buildings["layer"] == layer_name]
        commercial_layer["sector_code"] = sector_code
        commercial_layer["subsector_code"] = subsector_code
        commercial_layers.append(commercial_layer)

        print(f"* Done with points layer - {layer_name}")

    commercial_layers = pd.concat(commercial_layers, axis=0, ignore_index=True)
    commercial_layers_ids = list(set(commercial_layers["osm_id"].values.tolist()))
    commercial_buildings = commercial_buildings[
        ~commercial_buildings["osm_id"].isin(commercial_layers_ids)
    ]
    commercial_buildings = pd.concat(
        [commercial_layers, commercial_buildings], axis=0, ignore_index=True
    )
    with_geom = commercial_buildings[["osm_id", "geometry"]].drop_duplicates(
        subset=["osm_id"], keep="first"
    )
    without_geom = commercial_buildings.drop("geometry", axis=1)
    without_geom = (
        without_geom.set_index("osm_id")
        .astype(str)
        .groupby("osm_id")
        .agg(",".join)
        .reset_index()
    )

    mcs = [c for c in without_geom.columns.values.tolist() if c != "osm_id"]
    for mc in mcs:
        without_geom[mc] = without_geom.progress_apply(
            lambda x: join_values_in_column(x, mc), axis=1
        )

    commercial_buildings = gpd.GeoDataFrame(
        pd.merge(without_geom, with_geom, how="left", on=["osm_id"]),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        commercial_buildings.to_file(
            os.path.join(
                building_data_path,
                "commercial_buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            driver="GPKG",
        )

    commercial_buildings["assigned_attribute"] = commercial_buildings["layer"]
    commercial_buildings["infra_type"] = "X"
    commercial_buildings = commercial_buildings[
        ["osm_id", "sector_code", "subsector_code", "assigned_attribute", "infra_type"]
    ]

    for column in ["sector_code", "subsector_code", "assigned_attribute", "infra_type"]:
        commercial_buildings.rename(columns={column: f"{column}_match"}, inplace=True)
        buildings_input.rename(columns={column: f"{column}_building"}, inplace=True)

    buildings_input = match_buildings_to_polygon_dataframe(
        buildings_input, commercial_buildings, epsg=epsg_jamaica, spatial_join=False
    )
    matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        write_df = gpd.GeoDataFrame(
            pd.concat(matched_buildings + [buildings_input], axis=0, ignore_index=True),
            geometry="geometry",
            crs=f"EPSG:{epsg_jamaica}",
        )
        write_df.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="commercial",
            driver="GPKG",
        )
        del write_df

    """Step 7: Match the ports and airports
    """
    # buildings_input = gpd.read_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'),
    #                               layer="commercial")

    ports = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "port_polygon.gpkg"),
        layer="areas",
    )
    ports["sector_code"] = "I"
    ports["subsector_code"] = "600"
    ports["assigned_attribute"] = "Seaport"
    ports["infra_type"] = "Port"

    buildings_input = match_buildings_to_polygon_dataframe(
        buildings_input, ports, epsg=epsg_jamaica
    )

    matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    ports = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "airport_polygon.gpkg"
        ),
        layer="areas",
    )
    ports["sector_code"] = "I"
    ports["subsector_code"] = "600"
    ports["assigned_attribute"] = "Airport"
    ports["infra_type"] = "Airport"

    buildings_input = match_buildings_to_polygon_dataframe(
        buildings_input, ports, epsg=epsg_jamaica
    )
    matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    write_output = True
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        write_df = gpd.GeoDataFrame(
            pd.concat(matched_buildings + [buildings_input], axis=0, ignore_index=True),
            geometry="geometry",
            crs=f"EPSG:{epsg_jamaica}",
        )
        write_df.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="ports_and_airports",
            driver="GPKG",
        )
        del write_df

    """Step 8: Match the landuse layers for finding agricultural and residental lands
    """
    # matched_buildings = []
    # buildings_input = gpd.read_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'),
    #                               layer="ports_and_airports"
    #                               )

    # matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    # buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    land_use_types = gpd.read_file(
        os.path.join(
            processed_data_path,
            "land_type_and_use",
            "jamaica_land_use_combined_with_sectors.gpkg",
        ),
        layer="areas",
    )
    land_use_types = land_use_types.to_crs(epsg=epsg_jamaica)

    if "index_left" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_left", axis=1, inplace=True)
    elif "index_right" in land_use_types.columns.values.tolist():
        land_use_types.drop("index_right", axis=1, inplace=True)

    land_use_types["layer_id"] = land_use_types.index.values.tolist()
    # land_use_types["sector_code_res"] = "RES"
    # land_use_types["subsector_code"] = "RES"

    land_use_types = create_sector_subsector_attribute_columns(
        land_use_types,
        "layer_id",
        [
            ("sector_code_forest", "subsector_code_forest"),
            ("sector_code_tnc", "subsector_code_tnc"),
        ],
        [
            "sector_code_forest",
            "subsector_code_forest",
            "sector_code_tnc",
            "subsector_code_tnc",
        ],
        ["Classify", "NAME"],
        ["infra_sector_forest", "infra_sector_tnc"],
    )

    buildings_input = match_buildings_to_polygon_dataframe(
        buildings_input, land_use_types, epsg=epsg_jamaica
    )

    """Step 9: Match the residential building based on the area estimation
    """
    residential_buildings = gpd.read_file(
        os.path.join(incoming_data_path, "JAM_residential", "jam_residential.shp")
    )
    residential_buildings["sector_code"] = "RES"
    residential_buildings["subsector_code"] = "RES"
    residential_buildings["assigned_attribute"] = "Residentail"
    residential_buildings["infra_type"] = "X"

    residential_buildings = residential_buildings[
        ["osm_id", "sector_code", "subsector_code", "assigned_attribute", "infra_type"]
    ]

    for column in ["sector_code", "subsector_code", "assigned_attribute", "infra_type"]:
        residential_buildings.rename(columns={column: f"{column}_match"}, inplace=True)
        buildings_input.rename(columns={column: f"{column}_building"}, inplace=True)

    buildings_input = match_buildings_to_polygon_dataframe(
        buildings_input, residential_buildings, spatial_join=False
    )
    # matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    # buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    # buildings_input['area_sqm'] = buildings_input.progress_apply(lambda x:x.geometry.area,axis=1)

    write_output = True
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        write_df = gpd.GeoDataFrame(
            pd.concat(matched_buildings + [buildings_input], axis=0, ignore_index=True),
            geometry="geometry",
            crs=f"EPSG:{epsg_jamaica}",
        )
        write_df["area_sqm"] = write_df.progress_apply(
            lambda x: x.geometry.area, axis=1
        )
        write_df.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="final",
            driver="GPKG",
        )

    """Step 9: Final matching with manual assignment!
    """
    # matched_buildings = []
    # buildings_input = gpd.read_file(os.path.join(
    #                                 building_data_path,
    #                                 'buildings_assigned_economic_sectors_intermediate.gpkg'),
    #                               layer="final"
    #                               )

    # matched_buildings.append(buildings_input[buildings_input["sector_code"] != "X"])
    # buildings_input = buildings_input[buildings_input["sector_code"] == "X"]

    known_assignment = gpd.read_file(
        os.path.join(building_data_path, "manual_building_assign.shp")
    )
    known_assignment.rename(
        columns={
            "sector_cod": "sector_code",
            "subsector_": "subsector_code",
            "assigned_a": "assigned_attribute",
        },
        inplace=True,
    )
    known_assignment = known_assignment[
        ["osm_id", "sector_code", "subsector_code", "assigned_attribute", "infra_type"]
    ]

    for column in ["sector_code", "subsector_code", "assigned_attribute", "infra_type"]:
        known_assignment.rename(columns={column: f"{column}_match"}, inplace=True)
        buildings_input.rename(columns={column: f"{column}_building"}, inplace=True)

    buildings_input = match_buildings_to_polygon_dataframe(
        buildings_input, known_assignment, spatial_join=False
    )

    write_df = pd.concat(
        matched_buildings + [buildings_input], axis=0, ignore_index=True
    )

    # land_use_mapping = [ "Shrub","Under Construction","<Null>",
    #                     "Beach","Derelict Building","Drain",
    #                     "Drainage","Dump","Earth Drain",
    #                     "Grass Land","Grassland",
    #                     "Inaccessible","Land Cover-Ruinate",
    #                     "Land Use Designation to be Determined",
    #                     "Landfill","No Build Zone","Ruinate",
    #                     "Shrub","Shrub Woodland","Shrub-Woodland",
    #                     "Shrub_Woodland","Shurb",
    #                     "Subdivisional Development Land",
    #                     "Under Construction",
    #                     "Unknown",
    #                     "Vacant Building",
    #                     "Vacant Lot","Vacant Lots",
    #                     "Waterbody","dump","l"]
    land_use_mapping = [
        "Under Construction",
        "Derelict Building",
        "Drain",
        "Drainage",
        "Dump",
        "Earth Drain",
        "Inaccessible",
        "Land Cover-Ruinate",
        "Landfill",
        "No Build Zone",
        "Ruinate",
        "Under Construction",
        "Vacant Building",
        "Waterbody",
        "dump",
    ]
    land_use_mapping += ["abandoned"]

    write_df = gpd.read_file(
        os.path.join(
            building_data_path, "buildings_assigned_economic_sectors_intermediate.gpkg"
        ),
        layer="final",
    )

    write_df["remove_type"] = write_df.progress_apply(
        lambda x: possible_nonbuilding(x, land_use_mapping), axis=1
    )
    write_df[["to_remove", "remove"]] = write_df["remove_type"].apply(pd.Series)
    write_df.drop("remove_type", axis=1, inplace=True)

    write_output = True
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        write_df = gpd.GeoDataFrame(
            write_df, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
        )
        write_df.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="final",
            driver="GPKG",
        )

    """Step 11: Modify the fishery buildings
    """
    buildings_input = gpd.read_file(
        os.path.join(
            building_data_path, "buildings_assigned_economic_sectors_intermediate.gpkg"
        ),
        layer="final",
    )
    buildings_input["find_aqua"] = buildings_input.progress_apply(
        lambda x: 1 if "B" in x.sector_code else 0, axis=1
    )
    assigned_fishing = buildings_input[buildings_input["find_aqua"] == 1]
    all_fishing_assigned_ids = assigned_fishing["osm_id"].values.tolist()

    fishing_locations = gpd.read_file(
        os.path.join(processed_data_path, "land_type_and_use", "aqua_farms.gpkg"),
        layer="areas",
    )
    fishery_matches = gpd.sjoin(
        fishing_locations, buildings_input, how="inner", predicate="intersects"
    ).reset_index()
    fishing_ids = list(set(fishery_matches["osm_id"].values.tolist()))
    buildings_input.loc[buildings_input["osm_id"].isin(fishing_ids), "sector_code"] = (
        "B"
    )
    buildings_input.loc[
        buildings_input["osm_id"].isin(fishing_ids), "subsector_code"
    ] = "50"
    buildings_input.loc[
        buildings_input["osm_id"].isin(fishing_ids), "assigned_attribute"
    ] = "Aqua Farm"

    all_fishing_assigned_ids = [
        a for a in all_fishing_assigned_ids if a not in fishing_ids
    ]

    buildings_input = buildings_input[
        ~buildings_input["osm_id"].isin(all_fishing_assigned_ids)
    ]
    buildings_input.drop("find_aqua", axis=1, inplace=True)
    del all_fishing_assigned_ids

    assigned_fishing = assigned_fishing[assigned_fishing["sector_code"] != "B"]
    assigned_fishing["sector_code"] = assigned_fishing.progress_apply(
        lambda x: ",".join([s for s in str(x.sector_code).split(",") if s != "B"]),
        axis=1,
    )
    assigned_fishing["subsector_code"] = assigned_fishing.progress_apply(
        lambda x: ",".join(
            [s for s in str(x.subsector_code).split(",") if s not in (50, "50")]
        ),
        axis=1,
    )

    assigned_fishing.drop("find_aqua", axis=1, inplace=True)
    buildings_input = gpd.GeoDataFrame(
        pd.concat([buildings_input, assigned_fishing], axis=0, ignore_index=True),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    print(buildings_input)

    write_output = True
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        buildings_input.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="final_mod",
            driver="GPKG",
        )

    """Step 12: Add costs to buildings for damage calculations
    """
    buildings_input = gpd.read_file(
        os.path.join(
            building_data_path, "buildings_assigned_economic_sectors_intermediate.gpkg"
        ),
        layer="final_mod",
    )
    buildings_input = buildings_input.to_crs(epsg=epsg_jamaica)
    buildings_input = buildings_input[buildings_input["to_remove"] == 0]
    buildings_attriibutes = []
    buildings_sector_classes = []
    for row in buildings_input.itertuples():
        buildings_attriibutes += str(row.assigned_attribute).split(",")
        buildings_sector_classes += list(
            zip(str(row.sector_code).split(","), str(row.subsector_code).split(","))
        )

    buildings_attriibutes = pd.DataFrame(
        buildings_attriibutes, columns=["assigned_attribute"]
    )
    buildings_attriibutes = buildings_attriibutes.value_counts().reset_index(
        name="counts"
    )
    buildings_attriibutes.to_csv(
        os.path.join(building_data_path, "building_attributes.csv"), index=False
    )
    buildings_sector_classes = pd.DataFrame(
        buildings_sector_classes, columns=["sector_code", "subsector_code"]
    )
    buildings_sector_classes = buildings_sector_classes.value_counts().reset_index(
        name="counts"
    )
    buildings_sector_classes.to_csv(
        os.path.join(building_data_path, "buidling_sector_subsectors.csv"), index=False
    )

    costs_shapefiles = [
        "Construction_Permit_Mapping_Jan-March_Q4_2019-2020",
        "Construction_Permits_October_-_December_(Q3)_2019-2020",
        "Construction_Permits_Q1_2019_2020",
        "Q2-2019- Construction Permits",
    ]
    cost_columns = [
        ["Nature_of_", "Building_T", "Parish", "B_Area", "Est_Val"],
        ["Nature_of_", "Building_T", "Parish", "B_Area", "Est_Val"],
        ["Nat_Dev", "Bldg_Type", "Parish", "Bldg_Area", "Est_Val"],
        ["Nature_of", "Building_T", "Parish", "Building_A", "Estimated"],
    ]
    all_costs = []
    for c in range(len(costs_shapefiles)):
        shapefile = costs_shapefiles[c]
        columns = cost_columns[c]
        building_costs = gpd.read_file(
            os.path.join(
                incoming_data_path,
                "construction_permits",
                "Construction Permit Mapping 2019-2020 FY",
                f"{shapefile}.shp",
            )
        )[columns]
        building_costs.columns = [
            "development_type",
            "building_type",
            "parish",
            "area_sqm",
            "cost",
        ]
        building_costs = building_costs[
            (building_costs["area_sqm"] > 0) & (building_costs["cost"] > 0)
        ]
        all_costs.append(building_costs)

    all_costs = pd.concat(all_costs, axis=0, ignore_index=True)
    all_costs["cost_persqm"] = all_costs["cost"] / all_costs["area_sqm"]
    cost_ranges = quantiles(
        all_costs, ["development_type", "building_type"], ["cost_persqm"]
    )
    building_types = [
        "Commercial",
        "Industrial",
        "Institutional",
        "Mixed Use",
        "Other",
        "Recreation",
        "Residential",
        "Resort",
    ]
    cost_ranges = cost_ranges[
        (cost_ranges["development_type"] == "Erection")
        & (cost_ranges["building_type"].isin(building_types))
    ]
    print(cost_ranges)
    save_file = True
    if save_file == True:
        cost_ranges.to_csv(
            os.path.join(building_data_path, "building_cost_ranges.csv"), index=False
        )

    cost_ranges["min_damage_cost"] = cost_ranges["cost_persqm_q5"]
    cost_ranges["max_damage_cost"] = cost_ranges["cost_persqm_q95"]
    cost_ranges["mean_damage_cost"] = cost_ranges["cost_persqm_mean"]

    buildings_input["building_type"] = buildings_input.progress_apply(
        lambda x: assign_building_type(x), axis=1
    )
    buildings_input = pd.merge(
        buildings_input,
        cost_ranges[
            ["building_type", "min_damage_cost", "max_damage_cost", "mean_damage_cost"]
        ],
        how="left",
        on=["building_type"],
    )
    buildings_input["cost_unit"] = "JD/m2"
    write_output = False
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        buildings_input.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="final_costs",
            driver="GPKG",
        )

    """Step 13: Add populations assignments to buildings
    """

    buildings_input = gpd.read_file(
        os.path.join(
            building_data_path, "buildings_assigned_economic_sectors_intermediate.gpkg"
        ),
        layer="final_costs",
    )
    buildings_input = buildings_input.to_crs(epsg=epsg_jamaica)
    population_column = "2019"
    population = gpd.read_file(
        os.path.join(processed_data_path, "population", "population_projections.gpkg"),
        layer="mean",
    )
    print(population[population_column].sum())
    population = population.to_crs(epsg=epsg_jamaica)
    buildings_input["residential_type"] = buildings_input.progress_apply(
        lambda x: 1 if "RES" in str(x.sector_code) else 0, axis=1
    )
    residential_buildings = buildings_input[
        (buildings_input["area_sqm"] >= residential_min_area)
        & (buildings_input["residential_type"] == 1)
    ]

    population_buildings = gpd.sjoin(
        population[["ED_ID", population_column, "geometry"]],
        residential_buildings[
            [
                "osm_id",
                "sector_code",
                "subsector_code",
                "building_type",
                "area_sqm",
                "geometry",
            ]
        ],
        how="inner",
        predicate="intersects",
    ).reset_index()
    population_buildings.drop("geometry", axis=1, inplace=True)
    pop_sums = population_buildings.groupby(["ED_ID", "osm_id"]).agg(
        {"area_sqm": "sum"}
    )
    pop_frac = pop_sums.groupby(level=0).apply(lambda x: x / float(x.sum()))
    pop_frac = pop_frac.reset_index(level=["ED_ID", "osm_id"])
    print(pop_frac)
    residential_populations = pd.merge(
        population_buildings[["ED_ID", "osm_id", population_column]],
        pop_frac,
        how="left",
        on=["ED_ID", "osm_id"],
    )
    del pop_sums, pop_frac
    residential_populations["residential_population"] = (
        residential_populations[population_column] * residential_populations["area_sqm"]
    )
    residential_populations = residential_populations.round(
        {"residential_population": 0}
    )
    print(residential_populations["residential_population"].sum())
    residential_populations = (
        residential_populations.groupby("osm_id")["residential_population"]
        .sum()
        .reset_index()
    )
    print(residential_populations)
    buildings_input = pd.merge(
        buildings_input, residential_populations, how="left", on=["osm_id"]
    )
    buildings_input["residential_population"] = buildings_input[
        "residential_population"
    ].fillna(0)
    write_output = True
    if (
        write_output is True
    ):  # We can write the intermediary output for sense checking if we want to
        buildings_input.to_file(
            os.path.join(
                building_data_path,
                "buildings_assigned_economic_sectors_intermediate.gpkg",
            ),
            layer="final_costs_pop",
            driver="GPKG",
        )
        buildings_input.to_file(
            os.path.join(
                processed_data_path,
                "buildings",
                "buildings_assigned_economic_activity.gpkg",
            ),
            layer="areas",
            driver="GPKG",
        )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
