"""Create the potable water asset data for Jamaica
    Add the final asset data into a geopackage
"""

import os

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
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
    water_data_path = os.path.join(incoming_data_path, "water")
    potable_data_path = os.path.join(incoming_data_path, "water", "potable", "raw")

    # read in raw files

    # NOTE assumptions:
    # - pipes and storage tanks are assumed to not be susceptible to flooding/hurricane hazards and are therefore not provided with a cost or damage function
    # - water treatment plants are assigned a cost based on capacity: limited capacity data was supplemented by matching with system attributes from Parish water supply plans where available
    # - hurricane induced damage for pumping stations, intakes and wells are assumed to be the same as that for water treatment plants
    # - flood induced damage for intakes is assumed to be the same as that for pumping stations

    potable_facilities_NWC = gpd.read_file(
        os.path.join(potable_data_path, "potable_facilities.shp")
    )
    asset_population_served = pd.read_csv(
        os.path.join(potable_data_path, "asset_criticality_final.csv"), encoding="latin"
    )
    supplement_treatment_plant_capacity = pd.read_csv(
        os.path.join(
            water_data_path,
            "potable",
            "additional_capacity_data",
            "supplement_treatment_plant_capacity.csv",
        )
    )
    cost_data = pd.read_csv(
        os.path.join(water_data_path, "cost", "water_asset_costs.csv")
    )

    potable_facilities_NWC["lon"] = potable_facilities_NWC["geometry"].centroid.x
    potable_facilities_NWC["lat"] = potable_facilities_NWC["geometry"].centroid.y

    potable_facilities_NWC = pd.merge(
        potable_facilities_NWC,
        asset_population_served[
            [
                "OBJECTID",
                "LOCATION",
                "Type",
                "asset_pop_new_2010",
                "asset_pop_new_2020",
                "asset_pop_new_2030",
            ]
        ],
        on=["OBJECTID", "LOCATION", "Type"],
        how="left",
    )

    # potable_facilities_NWC = pd.merge(
    #     potable_facilities_NWC,
    #     asset_population_served[["OBJECTID", "path"]],
    #     on="OBJECTID",
    # ).rename(columns={"Type": "asset_type"})

    potable_facilities_NWC["Capacity"] = pd.to_numeric(
        potable_facilities_NWC["Capacity"], errors="coerce"
    )
    potable_facilities_NWC["capacity (mgd)"] = np.where(
        (potable_facilities_NWC["Type"].isin(["Filter Plant", "Treatment Plant"]))
        & (potable_facilities_NWC["Capacity"] > 1000000),
        potable_facilities_NWC["Capacity"] * 0.000001,
        potable_facilities_NWC["Capacity"],
    )

    # incorporate supplementary capacity data
    potable_facilities_NWC = pd.merge(
        potable_facilities_NWC,
        supplement_treatment_plant_capacity,
        on=["OBJECTID", "Type"],
        how="left",
    )
    potable_facilities_NWC["capacity (mgd)"] = np.where(
        potable_facilities_NWC["capacity (mgd)"] > 0,
        potable_facilities_NWC["capacity (mgd)"],
        potable_facilities_NWC["capacity (mgd) from supplement"],
    )
    potable_facilities_NWC["capacity (mgd)"] = pd.to_numeric(
        potable_facilities_NWC["capacity (mgd)"], errors="coerce"
    )

    # assign asset name for matching with cost/damage data
    potable_facilities_NWC["asset_type_cost_data"] = np.where(
        potable_facilities_NWC["Type"].isin(["Filter Plant", "Treatment Plant"]),
        "treatment plant",
        np.where(
            potable_facilities_NWC["Type"].isin(["River Source", "Spring", "Sump"]),
            "intake",
            np.where(
                potable_facilities_NWC["Type"].isin(
                    ["Relift Station", "Pump Station", "Booster Station"]
                ),
                "pumping unit",
                np.where(
                    potable_facilities_NWC["Type"].isin(["Production Well"]),
                    "well",
                    "na",
                ),
            ),
        ),
    )
    potable_facilities_NWC["asset_type_flood_damage"] = np.where(
        potable_facilities_NWC["Type"].isin(["Filter Plant", "Treatment Plant"]),
        "water treatment plant",
        np.where(
            potable_facilities_NWC["Type"].isin(["River Source", "Spring", "Sump"]),
            "pumping unit",
            np.where(
                potable_facilities_NWC["Type"].isin(
                    ["Relift Station", "Pump Station", "Booster Station"]
                ),
                "pumping unit",
                np.where(
                    potable_facilities_NWC["Type"].isin(["Production Well"]),
                    "well",
                    "na",
                ),
            ),
        ),
    )
    potable_facilities_NWC["asset_type_hurricane_damage"] = np.where(
        potable_facilities_NWC["Type"].isin(["Filter Plant", "Treatment Plant"]),
        "water treatment plant",
        np.where(
            potable_facilities_NWC["Type"].isin(["River Source", "Spring", "Sump"]),
            "water treatment plant",
            np.where(
                potable_facilities_NWC["Type"].isin(
                    ["Relift Station", "Pump Station", "Booster Station"]
                ),
                "water treatment plant",
                np.where(
                    potable_facilities_NWC["Type"].isin(["Production Well"]),
                    "water treatment plant",
                    "na",
                ),
            ),
        ),
    )
    # merge on cost data
    potable_facilities_NWC = pd.merge(
        potable_facilities_NWC,
        cost_data,
        left_on="asset_type_cost_data",
        right_on="asset",
        how="left",
    )

    potable_facilities_NWC["capacity_inferred"] = (
        potable_facilities_NWC["capacity (mgd) from supplement"]
        .astype(float)
        .fillna(0.1)
    )
    # print (potable_facilities_NWC[['capacity (mgd) from supplement','capacity_inferred']])
    potable_facilities_NWC["cost_and_units"] = potable_facilities_NWC.progress_apply(
        lambda x: estimate_costs_and_units(x), axis=1
    )
    potable_facilities_NWC[["cost_unit", "min_damage_cost", "max_damage_cost"]] = (
        potable_facilities_NWC["cost_and_units"].apply(pd.Series)
    )
    potable_facilities_NWC["min_damage_cost"] = potable_facilities_NWC[
        "min_damage_cost"
    ].fillna(0)
    potable_facilities_NWC["max_damage_cost"] = potable_facilities_NWC[
        "max_damage_cost"
    ].fillna(0)
    potable_facilities_NWC.drop(["cost_and_units"], axis=1, inplace=True)

    potable_facilities_NWC.rename(columns={"Type": "asset_type"}, inplace=True)

    # provide id
    potable_facilities_NWC["node_id"] = potable_facilities_NWC.apply(
        lambda node: f"{node.asset_type}_{node.OBJECTID}", axis=1
    )

    # export as gpkg
    geometry = [
        Point(xy)
        for xy in zip(
            potable_facilities_NWC["lon"].round(4).values,
            potable_facilities_NWC["lat"].round(4).values,
        )
    ]
    potable_facilities_NWC = gpd.GeoDataFrame(
        potable_facilities_NWC, crs=f"EPSG:{epsg_jamaica}", geometry=geometry
    )
    potable_facilities_NWC.to_file(
        os.path.join(
            processed_data_path, "networks", "water", "potable_facilities_NWC.gpkg"
        ),
        layer="nodes",
        driver="GPKG",
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
