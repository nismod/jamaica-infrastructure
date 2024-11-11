"""Estimate indirect losses from single point failure analysis

"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from analysis_utils import *
from tqdm import tqdm

tqdm.pandas()


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    results_dir = os.path.join(
        output_data_path, "economic_losses", "single_failure_scenarios"
    )
    """Read the buildings data which contains the GDP values attributed to buildings
    """
    buildings = gpd.read_file(
        os.path.join(
            processed_data_path,
            "buildings",
            "buildings_assigned_economic_activity.gpkg",
        ),
        layer="areas",
    )
    buildings["osm_id"] = buildings["osm_id"].astype(int)
    buildings.drop("E_GDP", axis=1, inplace=True)
    buildings_gdp_columns = [
        c for c in buildings.columns.values.tolist() if "_GDP" in c and c != "total_GDP"
    ]
    buildings = buildings[["osm_id"] + buildings_gdp_columns]

    print("* Done with buildings")
    """Read the Potable water data to get the water specific GDP and dependency linkages to buildings
    """
    potable_economic_activity_buildings = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "potable_facilities_buildings_economic_activity_mapping.csv",
        )
    )
    potable_economic_activity_buildings["osm_id"] = potable_economic_activity_buildings[
        "osm_id"
    ].astype(int)
    potable_economic_activity_buildings = (
        potable_economic_activity_buildings.drop_duplicates(
            subset=["node_id", "osm_id"], keep="first"
        )
    )
    potable_economic_activity_buildings = pd.merge(
        potable_economic_activity_buildings[["node_id", "osm_id"]],
        buildings[["osm_id"] + buildings_gdp_columns],
        how="left",
        on=["osm_id"],
    ).fillna(0)
    # del buildings

    """Add the potable water-building GDP values and merge with the water asset GDPs to get failure results     
    """
    potable_economic_activity_dependent = (
        potable_economic_activity_buildings.groupby(["node_id"])[buildings_gdp_columns]
        .sum()
        .reset_index()
    )

    potable_economic_activity_buildings = potable_economic_activity_buildings[
        ["node_id", "osm_id"]
    ]

    potable_economic_activity = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "potable_facilities_dependent_economic_activity.csv",
        )
    )
    potable_economic_activity = (
        potable_economic_activity.groupby(["node_id"])["W_GDP"].sum().reset_index()
    )
    potable_economic_activity = pd.merge(
        potable_economic_activity[["node_id", "W_GDP"]],
        potable_economic_activity_dependent,
        how="left",
        on=["node_id"],
    ).fillna(0)
    del potable_economic_activity_dependent
    potable_economic_activity["economic_loss"] = potable_economic_activity[
        ["W_GDP"] + buildings_gdp_columns
    ].sum(axis=1)
    potable_economic_activity["loss_unit"] = "JD/day"
    potable_economic_activity.to_csv(
        os.path.join(
            results_dir, "single_point_failure_potable_facilities_economic_losses.csv"
        ),
        index=False,
    )
    print("* Done with potable water")
    """Get the linkages between pipelines and potable water assets, and estimate pipeline dependent GDPs
    """
    pipelines_economic_activity = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "potable_pipelines_dependent_economic_activity.csv",
        )
    )
    pipelines_economic_activity = pd.merge(
        pipelines_economic_activity[["edge_id", "node_id"]],
        potable_economic_activity,
        how="left",
        on=["node_id"],
    ).fillna(0)
    pipelines_gdp_columns = [
        c
        for c in pipelines_economic_activity.columns.values.tolist()
        if c not in ["edge_id", "loss_unit"]
    ]
    pipelines_economic_activity = (
        pipelines_economic_activity.groupby(["edge_id", "loss_unit"])[
            pipelines_gdp_columns
        ]
        .sum()
        .reset_index()
    )
    # pipelines_economic_activity.drop("node_id",axis=1,inplace=True)
    pipelines_economic_activity.to_csv(
        os.path.join(
            results_dir, "single_point_failure_potable_pipelines_economic_losses.csv"
        ),
        index=False,
    )
    del pipelines_economic_activity
    print("* Done with pipelines")
    """Get the Irrigation node assets with the Agriculture dependent GDP 
    """

    irrigation_economic_activity = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "irrigation_nodes_dependent_economic_activity.csv",
        )
    )
    irrigation_economic_activity.rename(columns={"GVA (JD/day)": "A_GDP"}, inplace=True)
    irrigation_economic_activity["economic_loss"] = irrigation_economic_activity[
        "A_GDP"
    ]
    irrigation_economic_activity["loss_unit"] = "JD/day"
    irrigation_economic_activity.to_csv(
        os.path.join(
            results_dir, "single_point_failure_irrigation_nodes_economic_losses.csv"
        ),
        index=False,
    )
    print("* Done with irrigation nodes")
    """Get the Irrigation edge assets with the Agriculture dependent GDP 
    """
    irrigation_edges_economic_activity = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "irrigation_edges_dependent_economic_activity.csv",
        )
    )
    irrigation_edges_economic_activity.rename(
        columns={"GVA (JD/day)": "economic_loss"}, inplace=True
    )
    irrigation_edges_economic_activity = (
        irrigation_edges_economic_activity.groupby(["edge_id"])["economic_loss"]
        .sum()
        .reset_index()
    )
    irrigation_edges_economic_activity["loss_unit"] = "JD/day"
    # irrigation_edges_economic_activity.drop("node_id",axis=1,inplace=True)
    irrigation_edges_economic_activity.to_csv(
        os.path.join(
            results_dir, "single_point_failure_irrigation_edges_economic_losses.csv"
        ),
        index=False,
    )

    del irrigation_edges_economic_activity
    print("* Done with irrigation edges")

    """Electricity without water
        Get the electricity sink node assets and their dependent GDPs
        Look-up the results of the single failures of the node assets and link them to the sink nodes 
    """
    electricity_economic_activity = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "electricity_dependent_economic_activity.csv",
        )
    )

    gdp_columns = [
        c
        for c in electricity_economic_activity.columns.values.tolist()
        if c not in ("id", "total_GDP", "GDP_unit")
    ]

    electricity_nodes_failure_results = pd.read_csv(
        os.path.join(
            output_data_path,
            "electricity_failures",
            "single_point_failure_results_nodes.csv",
        )
    )[["attacked_node_id", "affected_node_id"]]
    failure_duplicates = electricity_nodes_failure_results[["attacked_node_id"]]
    failure_duplicates["affected_node_id"] = failure_duplicates["attacked_node_id"]
    failure_duplicates = failure_duplicates.drop_duplicates(
        subset=["attacked_node_id"], keep="first"
    )
    electricity_nodes_failure_results = pd.concat(
        [electricity_nodes_failure_results, failure_duplicates],
        axis=0,
        ignore_index=True,
    )
    del failure_duplicates

    electricity_nodes_failures = pd.merge(
        electricity_nodes_failure_results,
        electricity_economic_activity,
        how="left",
        left_on=["affected_node_id"],
        right_on=["id"],
    ).fillna(0)
    electricity_nodes_failures = (
        electricity_nodes_failures.groupby(["attacked_node_id"])[gdp_columns]
        .sum()
        .reset_index()
    )
    electricity_nodes_failures["economic_loss"] = electricity_nodes_failures[
        gdp_columns
    ].sum(axis=1)
    electricity_nodes_failures["loss_unit"] = "JD/day"
    electricity_nodes_failures.to_csv(
        os.path.join(
            results_dir,
            "single_point_failure_electricity_nodes_economic_losses_no_water.csv",
        ),
        index=False,
    )
    electricity_nodes_failures = electricity_nodes_failures[
        ["attacked_node_id", "E_GDP", "T_GDP"]
    ]
    electricity_nodes_failures.rename(columns={"attacked_node_id": "id"}, inplace=True)
    print("* Done with electricity nodes")
    # del electricity_nodes_failures
    """Electricity without water
        Get the electricity sink node assets and their dependent GDPs
        Look-up the results of the single failures of the edge assets and link them to the sink nodes 
    """
    electricity_edges_failure_results = pd.read_csv(
        os.path.join(
            output_data_path,
            "electricity_failures",
            "single_point_failure_results_edges.csv",
        )
    )[["attacked_edge_id", "affected_node_id"]]

    electricity_edges_failures = pd.merge(
        electricity_edges_failure_results,
        electricity_economic_activity,
        how="left",
        left_on=["affected_node_id"],
        right_on=["id"],
    ).fillna(0)
    electricity_edges_failures = (
        electricity_edges_failures.groupby(["attacked_edge_id"])[gdp_columns]
        .sum()
        .reset_index()
    )
    electricity_edges_failures["economic_loss"] = electricity_edges_failures[
        gdp_columns
    ].sum(axis=1)
    electricity_edges_failures["loss_unit"] = "JD/day"
    electricity_edges_failures.to_csv(
        os.path.join(
            results_dir,
            "single_point_failure_electricity_edges_economic_losses_no_water.csv",
        ),
        index=False,
    )
    # electricity_edges_failures = electricity_edges_failures[["attacked_node_id","E_GDP","T_GDP"]]
    # electricity_edges_failures.rename(columns={"attacked_node_id":"id"},inplace=True)
    del electricity_edges_failures
    print("* Done with electricity edges")

    # electricity_economic_activity = electricity_economic_activity[["id","E_GDP","T_GDP"]]
    """Electricity-water failures combined
        Get the electricity and water mappings and trace the dependencies towards agriculture and buildings
    """
    # Read the mappings
    electricity_water_mapping = pd.read_csv(
        os.path.join(
            processed_data_path, "networks/energy", "mapping_water_to_electricity.csv"
        )
    )
    electricity_water_mapping.rename(columns={"elec_node_id": "id"}, inplace=True)
    # Estimate the irrigation specific dependencies
    electricity_irrigation_mapping = pd.merge(
        electricity_water_mapping,
        irrigation_economic_activity[["node_id", "A_GDP"]],
        how="left",
        left_on=["water_node_id"],
        right_on=["node_id"],
    ).fillna(0)
    electricity_irrigation_mapping.drop("node_id", axis=1, inplace=True)
    electricity_irrigation_mapping = (
        electricity_irrigation_mapping.groupby(["id"])["A_GDP"].sum().reset_index()
    )
    electricity_irrigation_mapping.rename(columns={"A_GDP": "Irr_GDP"}, inplace=True)
    del irrigation_economic_activity
    print("* Done with electricity irrigation mapping")

    # Estimate the potable water specific dependencies
    electricity_potable_mapping = pd.merge(
        electricity_water_mapping,
        potable_economic_activity[["node_id", "W_GDP"]],
        how="left",
        left_on=["water_node_id"],
        right_on=["node_id"],
    ).fillna(0)
    electricity_potable_mapping.drop("node_id", axis=1, inplace=True)
    electricity_potable_mapping = (
        electricity_potable_mapping.groupby(["id"])["W_GDP"].sum().reset_index()
    )
    del potable_economic_activity
    print("* Done with electricity potable mapping")

    # Estimate the common buildings for the electricity and water assets

    electricity_water_buildings = pd.merge(
        electricity_water_mapping,
        potable_economic_activity_buildings,
        how="left",
        left_on=["water_node_id"],
        right_on=["node_id"],
    ).fillna(0)
    electricity_water_buildings = electricity_water_buildings[
        electricity_water_buildings["osm_id"] > 0
    ][["id", "osm_id"]]
    del electricity_water_mapping

    electricity_economic_activity_buildings = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "electricity_buildings_economic_activity_mapping.csv",
        )
    )[["id", "osm_id"]]
    electricity_economic_activity_buildings["osm_id"] = (
        electricity_economic_activity_buildings["osm_id"].astype(int)
    )

    electricity_nodes_failure_results = pd.merge(
        electricity_nodes_failure_results,
        electricity_economic_activity_buildings,
        how="left",
        left_on=["affected_node_id"],
        right_on=["id"],
    ).fillna(0)
    electricity_nodes_failure_results = electricity_nodes_failure_results[
        ["attacked_node_id", "osm_id"]
    ]
    electricity_nodes_failure_results.rename(
        columns={"attacked_node_id": "id"}, inplace=True
    )
    electricity_nodes_failure_results = pd.concat(
        [electricity_nodes_failure_results, electricity_water_buildings],
        axis=0,
        ignore_index=True,
    )
    electricity_nodes_failure_results = (
        electricity_nodes_failure_results.drop_duplicates(
            subset=["id", "osm_id"], keep="first"
        )
    )
    electricity_nodes_failure_results = pd.merge(
        electricity_nodes_failure_results, buildings, how="left", on=["osm_id"]
    ).fillna(0)
    electricity_failures = (
        electricity_nodes_failure_results.groupby(["id"])[buildings_gdp_columns]
        .sum()
        .reset_index()
    )
    print("* Done with electricity-water-buildings mapping")

    # Merge different results together
    electricity_failures = pd.merge(
        electricity_failures, electricity_irrigation_mapping, how="left", on=["id"]
    ).fillna(0)
    electricity_failures = pd.merge(
        electricity_failures, electricity_potable_mapping, how="left", on=["id"]
    ).fillna(0)
    electricity_failures = pd.merge(
        electricity_failures, electricity_nodes_failures, how="left", on=["id"]
    ).fillna(0)

    gdp_columns = buildings_gdp_columns + ["Irr_GDP", "W_GDP", "E_GDP", "T_GDP"]
    # print (electricity_failures.columns.values.tolist())
    electricity_failures["economic_loss"] = electricity_failures[gdp_columns].sum(
        axis=1
    )
    electricity_failures["loss_unit"] = "JD/day"
    electricity_failures.to_csv(
        os.path.join(
            results_dir, "single_point_failure_electricity_nodes_economic_losses.csv"
        ),
        index=False,
    )
    print("* Done with electricity nodes failures with water disruptions")
    electricity_edges_failures = pd.merge(
        electricity_edges_failure_results,
        electricity_failures,
        how="left",
        left_on=["affected_node_id"],
        right_on=["id"],
    ).fillna(0)
    electricity_edges_failures = (
        electricity_edges_failures.groupby(["attacked_edge_id"])[gdp_columns]
        .sum()
        .reset_index()
    )
    electricity_edges_failures.rename(columns={"attacked_edge_id": "id"}, inplace=True)
    electricity_edges_failures["economic_loss"] = electricity_edges_failures[
        gdp_columns
    ].sum(axis=1)
    electricity_edges_failures["loss_unit"] = "JD/day"
    electricity_edges_failures.to_csv(
        os.path.join(
            results_dir, "single_point_failure_electricity_edges_economic_losses.csv"
        ),
        index=False,
    )
    print("* Done with electricity edges failures with water disruptions")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
