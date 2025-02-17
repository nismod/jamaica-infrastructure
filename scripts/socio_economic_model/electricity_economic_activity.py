"""Assign Electricity dependent GDP for Jamaica
"""

import sys
import os
import pandas as pd
import geopandas as gpd

# gpd._compat.USE_PYGEOS = True
# gpd.options.use_pygeos = True
import numpy as np
from utils import *
from tqdm import tqdm

tqdm.pandas()


def match_buildings_to_footprints(buildings, gdf, building_id, gdf_id):
    """Function to assign buildings to output areas
    We do a spatial join to find all buildings intersecting output areas
    We then find the areas of intersection
    For boudary conditions, where a building might intersect more than one output area
        Assign the building to the output area which contains it most, and discard the others

    """
    matches = gpd.sjoin(
        buildings, gdf, how="inner", predicate="intersects"
    ).reset_index()
    matches.rename(columns={"geometry": "building_geometry"}, inplace=True)
    matches = pd.merge(matches, gdf[[gdf_id, "geometry"]], how="left", on=gdf_id)
    matches["area_match"] = matches.progress_apply(
        lambda x: x["building_geometry"].intersection(x["geometry"].buffer(0)).area,
        axis=1,
    )
    matches = matches.sort_values(by=["area_match"], ascending=False)
    matches = matches.drop_duplicates(subset=[building_id], keep="first")
    matches.drop(["area_match", "building_geometry", "geometry"], axis=1, inplace=True)

    return matches


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]

    epsg_jamaica = 3448
    financial_year = 2019  # The year for which we are estimating the GDP values
    population_column = (
        "population"  # The name of the population column for the electricity layer
    )

    # Macroeconomic sector classes assiigned for Jamaica
    sector_codes = [
        "A",
        "B",
        "C",
        "D",
        "F",
        "G",
        "H",
        "I",
        "J",
        "K",
        "L",
        "M",
        "N",
        "O",
    ]
    gdp_columns = [
        f"{scode}_GDP" for scode in sector_codes
    ]  # Names of GDP value columns in the buildings footprint data
    sector_code = "E"  # Electricity sector code
    subsector_code = "401"  # Electricity subsector code

    # Get the Electricity output areas
    electricity_id = "id"
    electricity_output_areas = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "energy", "electricity_network_v3.0.gpkg"
        ),
        layer="voronoi",
    )
    electricity_output_areas = electricity_output_areas.to_crs(epsg=epsg_jamaica)
    print(electricity_output_areas)

    """Step 1: Find the overall GDP generated by the electricity sector itself based on: 
            The Gross Value Added (GVA) output of the electricity sector
            GDP = GVA + TAX, where TAX = TAXES LESS SUBSIDIES ON PRODUCTS

            We only know the TAX as a fraction of the whole GVA
            We assume the same fraction of TAX (tax rate) is added to the GVA of each sector
                to get the sector GDP

            Assign this GDP to electricity output areas based on population weightage of each output area 
            Estimate Daily GDP by converting Annual values in J$ millions to daily values in J$
    """

    economic_output = pd.read_excel(
        os.path.join(
            processed_data_path,
            "macroeconomic_data",
            "detailed_sector_GVA_GDP_current_prices.xlsx",
        ),
        sheet_name="2019",
    )
    economic_output.columns = [
        str(c).strip() for c in economic_output.columns.values.tolist()
    ]
    economic_output["subsector_code"] = economic_output["subsector_code"].apply(str)
    # Get total GVA for the whole economic of Jamaica
    totat_gva = economic_output[economic_output["sector_code"] == "GVA"][
        f"{financial_year}"
    ].sum()
    # Get total TAX for the whole economic of Jamaica
    total_tax = economic_output[economic_output["sector_code"] == "TAX"][
        f"{financial_year}"
    ].sum()
    # Estimate TAX rate
    tax_rate = 1.0 * total_tax / totat_gva

    # Converting Annual values in J$ millions to daily values in J$
    output = (
        (1 + tax_rate)
        * (1.0e6 / 365.0)
        * economic_output[
            (economic_output["sector_code"] == sector_code)
            & (economic_output["subsector_code"] == subsector_code)
        ][f"{financial_year}"].sum()
    )
    # Disaggregate to output areas in proportion to population
    electricity_output_areas[f"{sector_code}_GDP"] = (
        output
        * electricity_output_areas[population_column]
        / electricity_output_areas[population_column].sum()
    )
    print(electricity_output_areas)
    electricity_output_areas["GDP_unit"] = "JD/day"
    # Write the GDP values back to the GPKG file so that we do not need to calcuate it again
    electricity_output_areas.to_file(
        os.path.join(
            processed_data_path, "networks", "energy", "electricity_network_v3.0.gpkg"
        ),
        layer="voronoi",
        driver="GPKG",
    )
    """Step 2: Find the Daily GDP of other sectors dependent upon electricity
            This is done by mapping the assets of other sectors to electricity output areas, having pre-computed their daily GDP
            At present we have two dependeny cases:
            Electricity - Telecoms: 
                Where we find the telecoms masts within electricity output areas to map T_GDP dependent on electricity
            Electricity - Buildings:
                Where we find the buildings within electricity output areas to map other GDP's dependent on electricity

    """
    # Electricity-Telecoms dependencies
    telecoms = gpd.read_file(
        os.path.join(processed_data_path, "networks", "telecoms", "telecoms.gpkg"),
        layer="nodes",
    )
    telecoms = telecoms.to_crs(epsg=epsg_jamaica)
    telecoms.rename(
        columns={"node_id": "telecoms_id"}, inplace=True
    )  # Not to confuse with electricity ID column name

    # Find the telecoms assets within the electricity output areas
    electricity_telecoms = gpd.sjoin(
        telecoms[["telecoms_id", "T_GDP", "geometry"]],
        electricity_output_areas[[electricity_id, "geometry"]],
        how="inner",
        predicate="within",
    ).reset_index()
    electricity_telecoms = electricity_telecoms[
        [electricity_id, "telecoms_id", "T_GDP"]
    ]
    # Store the output for future analysis
    electricity_telecoms.to_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "electricity_telecoms_economic_activity_mapping.csv",
        ),
        index=False,
    )
    print(electricity_telecoms)
    # Electricity-Buildings dependencies
    buildings = gpd.read_file(
        os.path.join(
            processed_data_path,
            "buildings",
            "buildings_assigned_economic_activity.gpkg",
        ),
        layer="areas",
    )
    buildings = buildings.to_crs(epsg=epsg_jamaica)
    building_id = "osm_id"
    electricity_buildings = match_buildings_to_footprints(
        buildings,
        electricity_output_areas[[electricity_id, "geometry"]],
        building_id,
        electricity_id,
    )
    electricity_buildings = electricity_buildings[
        [electricity_id, building_id] + gdp_columns
    ]
    # Store the output for future analysis
    electricity_buildings.to_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "electricity_buildings_economic_activity_mapping.csv",
        ),
        index=False,
    )

    # Get the aggreagted independent and dependent GDP estimate for electricity output areas
    electricity_output_areas = (
        electricity_output_areas.groupby([electricity_id])[f"{sector_code}_GDP"]
        .sum()
        .reset_index()
    )
    electricity_telecoms = (
        electricity_telecoms.groupby([electricity_id])["T_GDP"].sum().reset_index()
    )
    electricity_buildings = (
        electricity_buildings.groupby([electricity_id])[gdp_columns].sum().reset_index()
    )

    electricity_output_areas = pd.merge(
        electricity_output_areas[[electricity_id, f"{sector_code}_GDP"]],
        electricity_telecoms,
        how="left",
        on=[electricity_id],
    ).fillna(0)
    electricity_output_areas = pd.merge(
        electricity_output_areas, electricity_buildings, how="left", on=[electricity_id]
    ).fillna(0)
    electricity_output_areas["total_GDP"] = electricity_output_areas[
        [f"{sector_code}_GDP", "T_GDP"] + gdp_columns
    ].sum(axis=1)
    electricity_output_areas["GDP_unit"] = "JD/day"
    electricity_output_areas.to_csv(
        os.path.join(
            processed_data_path,
            "networks_economic_activity",
            "electricity_dependent_economic_activity.csv",
        ),
        index=False,
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
