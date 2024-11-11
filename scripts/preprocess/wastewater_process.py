"""Create the wastewater asset data for Jamaica
    Add the final asset data into a geopackage
"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from preprocess_utils import *


def estimate_costs_and_units(x):
    if x.curve == "Y (/mgd)":
        # print (str(x['capacity (mgd) from supplement']), x['capacity_inferred'])
        min_damage_cost = float(
            str(x["cost ($J) - lower bound"]).replace(",", "")
        ) * float(x["capacity_inferred"])
        max_damage_cost = float(
            str(x["cost ($J) - upper bound"]).replace(",", "")
        ) * float(x["capacity_inferred"])
        cost_unit = "$J"
    else:
        min_damage_cost = float(str(x["cost ($J) - lower bound"]).replace(",", ""))
        max_damage_cost = float(str(x["cost ($J) - upper bound"]).replace(",", ""))
        cost_unit = "$J"

    return cost_unit, min_damage_cost, max_damage_cost


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448
    wastewater_data_path = os.path.join(incoming_data_path, "water", "wastewater")
    # NOTE assumptions:
    # - note sewers and storage tanks are assumed to not be susceptible to flooding/hurricane hazards and are therefore not provided with a cost or damage function
    # - wastewater treatment plants are assigned a cost based on capacity: limited capacity data was supplemented by data scraped from NWC website
    # - the cost of construction of wastewater treatment plants and pumping stations are assumed to be the same as for water treatment plants
    # - hurricane induced damage for water treatment plants is assumed to be the same as that for wastewater pumping stations and treatment plants
    # - assume cost of construction for sumps is equivalent to potable water intake and flood and hurricane damage curves are equivalent to potable pumping unit and water treatment plant respectively

    # read in raw files
    waste_water_facilities_NWC = gpd.read_file(
        os.path.join(wastewater_data_path, "raw", "waste_water_facilities_NWC.shp")
    ).rename(columns={"Type": "asset_type", "Capacity": "capacity"})
    cost_data = pd.read_csv(
        os.path.join(incoming_data_path, "water", "cost", "water_asset_costs.csv")
    )
    # Currently ignoring sewers
    # sewers = gpd.read_file(
    #     os.path.join(wastewater_data_path, "raw", "wGravityMain.shp")
    # )

    # convert capacity data from gd to mgd
    waste_water_facilities_NWC["capacity"] = pd.to_numeric(
        waste_water_facilities_NWC["capacity"], errors="coerce"
    )
    waste_water_facilities_NWC["capacity (mgd)"] = np.where(
        waste_water_facilities_NWC["asset_type"] == "WW Treatment Plant",
        waste_water_facilities_NWC["capacity"] * 0.000001,
        waste_water_facilities_NWC["capacity"],
    )

    # incorporate supplementary capacity data
    supplement_treatment_plant_capacity = pd.read_csv(
        os.path.join(
            wastewater_data_path,
            "additional_capacity_data",
            "supplement_treatment_plant_capacity.csv",
        )
    ).rename(columns={"Type": "asset_type"})

    # merge
    waste_water_facilities_NWC = pd.merge(
        waste_water_facilities_NWC,
        supplement_treatment_plant_capacity,
        how="left",
        on=["Name", "asset_type"],
    )

    # assign asset name for matching with cost/damage data
    asset_type_name_conversion = pd.DataFrame(
        {
            "asset_type": [
                "WW Relift Station",
                "WW Pump Station",
                "WW Treatment Plant",
                "Sump",
            ],
            "asset_type_cost_data": [
                "pumping unit",
                "pumping unit",
                "treatment plant",
                "intake",
            ],
            "asset_type_flood_damage": [
                "pumping unit",
                "pumping unit",
                "wastewater treatment plant",
                "pumping unit",
            ],
            "asset_type_hurricane_damage": [
                "water treatment plant",
                "water treatment plant",
                "water treatment plant",
                "water treatment plant",
            ],
        }
    )
    waste_water_facilities_NWC = pd.merge(
        waste_water_facilities_NWC,
        asset_type_name_conversion,
        how="left",
        on="asset_type",
    )

    waste_water_facilities_NWC = pd.merge(
        waste_water_facilities_NWC,
        cost_data,
        left_on="asset_type_cost_data",
        right_on="asset",
        how="left",
    )
    waste_water_facilities_NWC["capacity_inferred"] = (
        waste_water_facilities_NWC["capacity (mgd) from supplement"]
        .astype(float)
        .fillna(0.1)
    )
    # print (waste_water_facilities_NWC[['capacity (mgd) from supplement','capacity_inferred']])
    waste_water_facilities_NWC["cost_and_units"] = (
        waste_water_facilities_NWC.progress_apply(
            lambda x: estimate_costs_and_units(x), axis=1
        )
    )
    waste_water_facilities_NWC[["cost_unit", "min_damage_cost", "max_damage_cost"]] = (
        waste_water_facilities_NWC["cost_and_units"].apply(pd.Series)
    )
    waste_water_facilities_NWC.drop(["cost_and_units"], axis=1, inplace=True)

    # provide id
    waste_water_facilities_NWC["node_id"] = waste_water_facilities_NWC.apply(
        lambda node: f"{node.asset_type}_{node.OBJECTID}", axis=1
    )

    # waste_water_facilities_NWC.rename(columns={"curve": "cost_unit"},inplace=True)
    # export as gpkg
    waste_water_facilities_NWC = gpd.GeoDataFrame(
        waste_water_facilities_NWC,
        crs=f"EPSG:{epsg_jamaica}",
        geometry=waste_water_facilities_NWC["geometry"],
    )
    waste_water_facilities_NWC.to_file(
        os.path.join(
            processed_data_path, "networks", "water", "waste_water_facilities_NWC.gpkg"
        ),
        layer="nodes",
        driver="GPKG",
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
