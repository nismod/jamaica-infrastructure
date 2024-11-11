"""Create a road network for Jamaica
    Take road data from the NSDMB database and create topological networks
    Match the geometries and attributes between a NWA national roads network
    With more detailed road network
    Add the final network into a geopackage

"""

import sys
import os

import pandas as pd
import geopandas as gpd
from difflib import SequenceMatcher
import re
from shapely.geometry import Point
from preprocess_utils import *
from tqdm import tqdm

tqdm.pandas()


def modify_road_surface(x):
    if x["constructi"] in [None, np.nan]:
        x["constructi"] = ""
    if x["construc_1"] in [None, np.nan]:
        x["construc_1"] = ""

    x["constructi"] = str(x["constructi"])
    x["construc_1"] = str(x["construc_1"])

    if (x["constructi"]) or (x["construc_1"]):
        if x["construc_1"]:
            if x["construc_1"] in ["SD & AC", "AC & SD"]:
                return "SD & AC"
            else:
                return x["construc_1"]
        else:
            if x["construci"] in ["SD & AC", "AC & SD"]:
                return "SD & AC"
            else:
                return x["constructi"]
    else:
        return "Surface Dressed"


def modify_road_width(x, standard_width=7.3):
    if float(x["average_wi"]) > 0:
        return float(x["average_wi"])
    elif float(x["averagewid"]) > 0:
        return float(x["averagewid"])
    else:
        return standard_width


def nearest_name(x, names, x_column, name_column):
    gdf = names[~names[name_column].isna()]
    highest_match = 0
    for g in gdf.itertuples():
        m = SequenceMatcher(
            None, str(x[x_column]).lower(), str(getattr(g, name_column)).lower()
        ).ratio()
        if m >= highest_match:
            highest_match = m
            nearest_name = getattr(g, name_column)
            nearest_section = g.road_section
            nearest_parish = g.road_parish
        if highest_match == 1:
            break
    return nearest_name, nearest_section, nearest_parish, highest_match


def modify_road_construction(x):
    if x["road_construction"] in ["SD & AC", "AC & SD"]:
        return "SD & AC"
    else:
        return x["road_construction"]


def assign_costs(x, all_costs):
    road_surface = x["road_construction"]
    if road_surface in ["Asphaltic Concrete", "SD & AC", "Surface Dressed"]:
        road_costs = all_costs[~all_costs["cost_code"].isin(["4140"])]
    else:
        road_costs = all_costs[~all_costs["cost_code"].isin(["4130"])]

    reopen_costs = road_costs[road_costs["cost_code"].isin(["4110", "4120", "4310b"])]
    road_costs = road_costs.set_index(["cost_code", "cost_description"])

    asset_costs = (
        road_costs.sum(axis=0, numeric_only=True).values.tolist()
        + reopen_costs.sum(axis=0, numeric_only=True).values.tolist()
    )
    return tuple(asset_costs)


def assign_costs_roads(x, all_costs):
    road_surface = x["road_construction"]
    road_width = x["road_width"]
    if road_surface in ["Asphaltic Concrete", "SD & AC", "Surface Dressed"]:
        road_costs = all_costs[~all_costs["cost_code"].isin(["4140"])]
    else:
        road_costs = all_costs[~all_costs["cost_code"].isin(["4130"])]

    reopen_costs = road_costs[road_costs["cost_code"].isin(["4110", "4120", "4310b"])]
    road_costs = road_costs.set_index(["cost_code", "cost_description"])

    asset_costs = list(
        road_width * road_costs.sum(axis=0, numeric_only=True).values
    ) + list(road_width * reopen_costs.sum(axis=0, numeric_only=True).values)
    return tuple(asset_costs)


def assign_costs_bridges(x, all_costs):
    if x.asset_type == "bridge":
        road_costs = all_costs[~all_costs["cost_code"].isin(["4140"])]
        reopen_costs = road_costs[
            road_costs["cost_code"].isin(["4110", "4120", "4310b"])
        ]
        road_costs = road_costs.set_index(["cost_code", "cost_description"])
        asset_costs = (
            road_costs.sum(axis=0, numeric_only=True).values.tolist()
            + reopen_costs.sum(axis=0, numeric_only=True).values.tolist()
        )
        return tuple(asset_costs)


def find_bridge_in_string(x):
    if "bridge" in str(x.from_to).lower():
        return 1
    else:
        return 0


def create_road_surface_type(x):
    road_surface = x["road_construction"]
    if road_surface in ["Asphaltic Concrete", "SD & AC", "Surface Dressed"]:
        return "asphalt"
    else:
        return "non-asphalt"


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448

    road_cost_file_path = os.path.join(
        incoming_data_path, "Financial Data Post Flood Events"
    )
    road_cost_file = os.path.join(
        incoming_data_path,
        "Financial Data Post Flood Events",
        "TS Laura BILL SUBMISSIONS REPORT.xls",
    )

    data_columns = [
        "section",
        "from",
        "to",
        "description",
        "Pavement - Clearing of slips & blockage of roads",
        "Pavement - Temporary repairs to roads",
        "Pavement - Restoration of Asphalted Road",
        "Pavement - Restoration of Unasphalted Rd.",
        "Drainage - Cleaning of blocked drains etc.",
        "Drainage - Cleaning of critical drains",
        "Drainage - Repair/Kerb & Channel",
        "Drainage - Relay Culverts",
        "Construct & Reconstruction of R.R. wall",
        "total_estimate",
        "cost_to_reopen",
        "balance",
    ]
    xl = pd.ExcelFile(road_cost_file)

    road_sheets = [
        x for x in xl.sheet_names if x not in ["Bills", "SUMMARY", "Pictures"]
    ]

    road_costs = []
    for parish in road_sheets:
        road_cost_data = pd.read_excel(
            road_cost_file, sheet_name=parish, skiprows=6
        ).fillna(0)
        road_cost_data.columns = data_columns
        road_cost_data = road_cost_data[road_cost_data["total_estimate"] > 0]
        road_cost_data = road_cost_data[
            (road_cost_data["from"] != 0) | (road_cost_data["to"] != 0)
        ]
        road_cost_data["parish"] = parish
        road_costs.append(road_cost_data)

    road_costs_all = pd.concat(road_costs, axis=0, ignore_index=True)
    road_costs_all["cost_id"] = road_costs_all.index.values.tolist()
    print(road_costs_all)

    road_costs_all["section"] = road_costs_all.progress_apply(
        lambda x: re.sub("[/-]", "", str(x["section"])), axis=1
    )
    road_costs_all["from_to"] = road_costs_all.progress_apply(
        lambda x: f"{x['from']}-{x['to']}", axis=1
    )
    road_costs_all["contains_bridge"] = road_costs_all.progress_apply(
        lambda x: find_bridge_in_string(x), axis=1
    )

    """Step 0: Find the unit costs in $J/unit for a bridge and get the highest values as estimates
    """
    road_costs_bridges = road_costs_all[road_costs_all["contains_bridge"] == 1]
    data_columns = [
        "Pavement - Clearing of slips & blockage of roads",
        "Pavement - Temporary repairs to roads",
        "Pavement - Restoration of Asphalted Road",
        "Pavement - Restoration of Unasphalted Rd.",
        "Drainage - Cleaning of blocked drains etc.",
        "Drainage - Cleaning of critical drains",
        "Drainage - Repair/Kerb & Channel",
        "Drainage - Relay Culverts",
        "Construct & Reconstruction of R.R. wall",
    ]
    for d in data_columns:
        road_costs_bridges[d] = road_costs_bridges[d].replace(" ", "0.0")
        road_costs_bridges[d] = road_costs_bridges[d].astype(float)

    all_costs = []
    cost_columns_new = []
    for c in data_columns:
        # We only chose the top 5 values of each type because they indicate more severe damage effects
        costs = (
            road_costs_bridges[road_costs_bridges[c] > 0]
            .sort_values(by=c, ascending=False)[c]
            .head(5)
        )
        all_costs.append((c, costs.min(), costs.mean(), costs.max()))

    all_costs = pd.DataFrame(
        all_costs, columns=["cost_description", "min", "mean", "max"]
    )
    all_costs["cost_code"] = [
        "4110",
        "4120",
        "4130",
        "4140",
        "4310",
        "4310b",
        "4320",
        "4330",
        "4520",
    ]

    nodes = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="nodes",
    )
    nodes = gpd.GeoDataFrame(
        nodes, geometry="geometry", crs={"init": f"epsg:{epsg_jamaica}"}
    )
    # print (edges.crs)
    # edges.set_crs(epsg=epsg_jamaica, allow_override=True)
    # print (edges.crs)
    nodes["cost_unit"] = "$J"
    nodes["assigned_costs"] = nodes.progress_apply(
        lambda x: assign_costs_bridges(x, all_costs), axis=1
    )
    nodes[
        [
            "min_damage_cost",
            "mean_damage_cost",
            "max_damage_cost",
            "min_reopen_cost",
            "mean_reopen_cost",
            "max_reopen_cost",
        ]
    ] = nodes["assigned_costs"].apply(pd.Series)
    nodes.drop("assigned_costs", axis=1, inplace=True)
    print(nodes)
    nodes.to_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="nodes",
        driver="GPKG",
    )

    """Step 1: Match the ones with the same section ID values
    """
    road_costs = road_costs_all.copy()
    road_edges = gpd.read_file(
        os.path.join(
            incoming_data_path, "nsdmb", "GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb"
        ),
        layer="roads_main_NWA",
    )
    road_edges.to_crs(epsg=epsg_jamaica)
    # print (road_edges)
    # print (road_edges.columns)
    road_edges.columns = map(str.lower, road_edges.columns)

    road_edges.rename(
        columns={"parish": "road_parish", "first_clas": "road_class"}, inplace=True
    )
    string_columns = ["section", "constructi", "construc_1"]
    road_edges[string_columns] = road_edges[string_columns].replace(
        r"^\s*$", "", regex=True
    )
    # edges[string_columns] = edges[string_columns].fillna('',inplace=True)
    numeric_columns = ["averagewid", "average_wi"]
    for num in numeric_columns:
        road_edges[num] = road_edges[num].replace(" ", "0.0")
        road_edges[num] = road_edges[num].astype(float)

    road_edges.rename(columns={"from_": "from_road", "to": "to_road"}, inplace=True)
    road_edges["from_to_road"] = road_edges.progress_apply(
        lambda x: f"{x.from_road}-{x.to_road}", axis=1
    )
    road_edges["to_from_road"] = road_edges.progress_apply(
        lambda x: f"{x.to_road}-{x.from_road}", axis=1
    )
    road_edges["road_section"] = road_edges.progress_apply(
        lambda x: re.sub("[/-]", "", str(x["section"])), axis=1
    )
    road_edges["road_surface"] = road_edges.progress_apply(
        lambda x: modify_road_surface(x), axis=1
    )
    road_edges["road_width"] = road_edges.progress_apply(
        lambda x: modify_road_width(x), axis=1
    )
    road_edges["length_m"] = road_edges.progress_apply(
        lambda x: x.geometry.length, axis=1
    )
    print(road_edges.columns)
    road_edges = road_edges[
        [
            "road_section",
            "road_class",
            "str_name",
            "from_road",
            "to_road",
            "from_to_road",
            "to_from_road",
            "road_parish",
            "road_surface",
            "length_m",
            "road_width",
        ]
    ]
    road_costs = pd.merge(
        road_costs, road_edges, how="left", left_on="section", right_on="road_section"
    )

    """Step 2: Match the ones with similar names
    """
    for mt in ["from_to", "to_from"]:
        road_costs[f"{mt}_matches"] = road_costs.progress_apply(
            lambda x: nearest_name(x, road_edges, "from_to", f"{mt}_road"), axis=1
        )
        road_costs[[f"{mt}_match", f"{mt}_section", f"{mt}_parish", f"{mt}_value"]] = (
            road_costs[f"{mt}_matches"].apply(pd.Series)
        )
        road_costs.drop(f"{mt}_matches", axis=1, inplace=True)

    for mt in ["from", "to"]:
        road_costs[f"{mt}_matches"] = road_costs.progress_apply(
            lambda x: nearest_name(x, road_edges, mt, "str_name"), axis=1
        )
        road_costs[
            [
                f"{mt}_street_match",
                f"{mt}_street_section",
                f"{mt}_street_parish",
                f"{mt}_street_value",
            ]
        ] = road_costs[f"{mt}_matches"].apply(pd.Series)
        road_costs.drop(f"{mt}_matches", axis=1, inplace=True)

    """Step 3: Filter out the matches
    """
    threshold = 0.73
    section_matches = road_costs[~road_costs["road_section"].isna()][
        ["cost_id", "road_section"]
    ]
    from_to_matches = road_costs[road_costs["from_to_value"] >= threshold][
        ["cost_id", "from_to_section"]
    ]
    to_from_matches = road_costs[road_costs["to_from_value"] >= threshold][
        ["cost_id", "to_from_section"]
    ]
    threshold = 0.70
    from_matches = road_costs[road_costs["from_street_value"] >= threshold][
        ["cost_id", "from_street_section"]
    ]
    to_matches = road_costs[road_costs["to_street_value"] >= threshold][
        ["cost_id", "to_street_section"]
    ]

    section_matches = pd.merge(
        section_matches, road_costs_all, how="left", on="cost_id"
    )
    section_matches = pd.merge(
        section_matches, road_edges, how="left", on="road_section"
    )

    from_to_matches = pd.merge(
        from_to_matches, road_costs_all, how="left", on="cost_id"
    )
    from_to_matches = pd.merge(
        from_to_matches,
        road_edges,
        how="left",
        left_on="from_to_section",
        right_on="road_section",
    )

    to_from_matches = pd.merge(
        to_from_matches, road_costs_all, how="left", on="cost_id"
    )
    to_from_matches = pd.merge(
        to_from_matches,
        road_edges,
        how="left",
        left_on="to_from_section",
        right_on="road_section",
    )

    from_matches = pd.merge(from_matches, road_costs_all, how="left", on="cost_id")
    from_matches = pd.merge(
        from_matches,
        road_edges,
        how="left",
        left_on="from_street_section",
        right_on="road_section",
    )

    to_matches = pd.merge(to_matches, road_costs_all, how="left", on="cost_id")
    to_matches = pd.merge(
        to_matches,
        road_edges,
        how="left",
        left_on="to_street_section",
        right_on="road_section",
    )

    completed_matches = pd.concat(
        [section_matches, from_to_matches, to_from_matches, from_matches, to_matches],
        axis=0,
        ignore_index=True,
    )
    completed_matches.drop_duplicates(
        subset=["cost_id", "road_section"], keep="first", inplace=True
    )

    """Step 4: Find the unit costs in $J/m2 and get the highest values as estimates
    """
    data_columns = [
        "Pavement - Clearing of slips & blockage of roads",
        "Pavement - Temporary repairs to roads",
        "Pavement - Restoration of Asphalted Road",
        "Pavement - Restoration of Unasphalted Rd.",
        "Drainage - Cleaning of blocked drains etc.",
        "Drainage - Cleaning of critical drains",
        "Drainage - Repair/Kerb & Channel",
        "Drainage - Relay Culverts",
        "Construct & Reconstruction of R.R. wall",
    ]
    cost_columns = []
    for d in data_columns:
        completed_matches[d] = completed_matches[d].replace(" ", "0.0")
        completed_matches[d] = completed_matches[d].astype(float)
        completed_matches[f"{d} ($J/m2)"] = completed_matches[d] / (
            completed_matches["length_m"] * completed_matches["road_width"]
        )
        cost_columns.append(f"{d} ($J/m2)")

    all_costs = []
    cost_columns_new = []
    for c in cost_columns:
        # We only chose the top 5 values of each type because they indicate more severe damage effects
        costs = (
            completed_matches[completed_matches[c] > 0]
            .sort_values(by=c, ascending=False)[c]
            .head(5)
        )
        all_costs.append((c, costs.min(), costs.mean(), costs.max()))

    all_costs = pd.DataFrame(
        all_costs, columns=["cost_description", "min", "mean", "max"]
    )
    all_costs["cost_code"] = [
        "4110",
        "4120",
        "4130",
        "4140",
        "4310",
        "4310b",
        "4320",
        "4330",
        "4520",
    ]
    reopen_costs = all_costs[all_costs["cost_code"].isin(["4110", "4120", "4310b"])]
    all_costs = all_costs.set_index(["cost_code", "cost_description"])
    all_costs.loc[("Total", "Estimated total damage cost ($J/m2)"), :] = all_costs.sum(
        axis=0, numeric_only=True
    )
    all_costs.loc[
        ("Reopen", "Estimated cost to reopen (4110 + 4120 + 4310b) ($J/m2)"), :
    ] = reopen_costs.sum(axis=0, numeric_only=True)
    all_costs = all_costs.reset_index()

    # print (all_costs.sum(axis=0,numeric_only=True).values.tolist())
    all_costs[["cost_code", "cost_description", "min", "mean", "max"]].to_csv(
        os.path.join(road_cost_file_path, "roads_damage_cost_estimates.csv"),
        index=False,
    )
    print(all_costs[["cost_code", "cost_description", "min", "mean", "max"]])

    """Step 5: Integrate the costs with the road edges dataset
        We will use the road surface type to assign costs
        We will assign two types of costs:
            - Total damage cost
            - Cost to reopen the roads
    """
    all_costs = pd.read_csv(
        os.path.join(road_cost_file_path, "roads_damage_cost_estimates.csv")
    )
    edges = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="edges",
    )
    edges = gpd.GeoDataFrame(
        edges, geometry="geometry", crs={"init": f"epsg:{epsg_jamaica}"}
    )
    edges["asset_type"] = edges.progress_apply(
        lambda x: create_road_surface_type(x), axis=1
    )
    edges["road_length_m"] = edges.progress_apply(lambda x: x.geometry.length, axis=1)
    edges["road_construction"] = edges.progress_apply(
        lambda x: modify_road_construction(x), axis=1
    )
    edges["cost_unit"] = "$J/m"
    edges["assigned_costs"] = edges.progress_apply(
        lambda x: assign_costs_roads(x, all_costs), axis=1
    )
    edges[
        [
            "min_damage_cost",
            "mean_damage_cost",
            "max_damage_cost",
            "min_reopen_cost",
            "mean_reopen_cost",
            "max_reopen_cost",
        ]
    ] = edges["assigned_costs"].apply(pd.Series)
    edges.drop("assigned_costs", axis=1, inplace=True)
    # print (edges)
    edges.to_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="edges",
        driver="GPKG",
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
