"""Assign agriculture GDP to land use layers in Jamaica"""

import logging

import click
import geopandas as gpd
import pandas as pd

from jamaica_infrastructure.geo import (
    LOCAL_PROJ_CRS_EPSG,
    remove_geometry_collections,
)


@click.command()
@click.version_option("1.0")
@click.option(
    "--spam-agriculture-outputs-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, writable=True),
    help="Path to SPAM agriculture outputs GeoPackage.",
)
@click.option(
    "--crop-details-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to crop details CSV.",
)
@click.option(
    "--land-use-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to Jamaica land use combined with sectors GeoPackage.",
)
@click.option(
    "--fishing-locations-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to Jamaica fishing locations GeoPackage.",
)
@click.option(
    "--economic-output-path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to economic output CSV.",
)
@click.option(
    "--output-areas",
    required=True,
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    help="Path to output agriculture GDP GeoPackage.",
)
def main(
    spam_agriculture_outputs_path,
    crop_details_path,
    land_use_path,
    fishing_locations_path,
    economic_output_path,
    output_areas,
):
    """
    Process agricultural production data and allocate GDP to land use areas.

    Example usage:

        python workflow/1_context/agricultural_production.py \
            --spam-agriculture-outputs-path processed_data/agriculture_data/spam_agriculture_outputs.gpkg \
            --crop-details-path processed_data/agriculture_data/crop_details.NIP_2023_codes.csv \
            --land-use-path processed_data/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg \
            --fishing-locations-path processed_data/land_type_and_use/aqua_farms.gpkg \
            --economic-output-path processed_data/macroeconomic_data/NIP_2023.csv \
            --output-areas processed_data/agriculture_data/agriculture_gdp.gpkg

    This function performs the following steps:

    1. Reads SPAM (Spatial Production Allocation Model) crop areas
    2. Maps crop types to economic subsectors and computes production values for
       each crop area / sector
    3. Reads combined land use data and performs a spatial join with crop yield
       areas
    4. Allocates daily GDP (converted from annual JMD millions) to agricultural
       land parcels, weighted by crop production value and area, using economic
       output data filtered by sector code 'A'. Post-harvest output and
       agricultural services allocated by crop tonnage. Forest output allocated
       by forest area. Non-agricultural land use areas are assigned zero GDP.
    5. Computes per-square-meter GDP and crop tonnage densities for each land
       parcel.
    6. Writes areas GeoPackage with GDP and crop tonnage by land use area.
    """

    crop_yields = gpd.read_file(
        spam_agriculture_outputs_path,
        layer=f"tonnage_areas",
    )
    crop_yields = crop_yields.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)

    #
    # Map crop types to subsector codes
    #
    crop_details = pd.read_csv(crop_details_path)
    tech = ["A", "I", "R"]
    poultry_crops = ["maiz", "ocer", "pmil", "smil", "soyb", "sunf", "whea"]

    crop_yields["crop_tons"] = 0

    all_crop_columns = []
    all_sector_columns = []

    idx = crop_yields.index
    zero_series = pd.Series(0.0, index=idx)

    sector_value = {}
    crop_values = {}
    total_crop_volume = pd.Series(0.0, index=idx)

    for crop_detail in crop_details.itertuples():
        crop_columns = [f"{crop_detail.name.upper()}_{t}" for t in tech]
        sector_columns = [
            f"{crop_detail.sector_code}_{crop_detail.subsector_code}_{t}" for t in tech
        ]
        poultry_columns = [f"A_Animal Production_{t}" for t in tech]

        for crop, sector, poultry in zip(crop_columns, sector_columns, poultry_columns):
            all_sector_columns += [sector, poultry]

            if sector not in sector_value:
                sector_value[sector] = zero_series.copy()
            if poultry not in sector_value:
                sector_value[poultry] = zero_series.copy()

            crop_volume = (
                crop_yields[crop].clip(lower=0)
                if crop in crop_yields.columns
                else zero_series
            )
            crop_value = crop_detail.value_usd_per_ton * crop_volume
            logging.debug(crop, sector, crop_value.sum())
            total_crop_volume = total_crop_volume.add(crop_volume, fill_value=0)
            sector_value[sector] = sector_value[sector].add(crop_value, fill_value=0)

            if crop_detail.name in poultry_crops:
                sector_value[poultry] = sector_value[poultry].add(
                    crop_value, fill_value=0
                )

            prod_col = f"{crop}_prod"
            crop_values[prod_col] = crop_value
            all_crop_columns.append(prod_col)

    sector_value_df = pd.DataFrame(sector_value, index=idx)
    crop_value_df = pd.DataFrame(crop_values, index=idx)

    crop_yields = pd.concat(
        [
            crop_yields[["crop_id", "areas", "geometry"]],
            (total_crop_volume / crop_yields["areas"]).rename("crop_tons_persqm"),
            sector_value_df[list(set(all_sector_columns))],
            crop_value_df[all_crop_columns],
        ],
        axis=1,
    )

    crop_yields.to_file(
        spam_agriculture_outputs_path,
        layer=f"values_per_sqm",
        driver="GPKG",
    )

    #
    # Read economic data
    #
    economic_output = pd.read_csv(economic_output_path, comment="#").query(
        "Code == 'A'"
    )
    value_colname = "GVA JMD Millions (2023)"

    #
    # Read landuse
    #
    landuse = gpd.read_file(land_use_path, layer="areas")
    landuse = landuse.to_crs(epsg=LOCAL_PROJ_CRS_EPSG)
    landuse["land_id"] = landuse.index.values.tolist()
    landuse["land_id"] = landuse.apply(lambda x: f"land_{x.land_id}", axis=1)

    # Split agricultural from non-ag
    mask = (landuse["sector_code_forest"] == "A") | (landuse["sector_code_tnc"] == "A")
    # Set non-agricultural value to zero
    non_agri_land_use = landuse[~mask].copy()
    agri_land_use = landuse[mask].copy()

    #
    # Handle subsectors
    #

    # JIC 2005 codes
    # Traditional Export Agriculture
    #    Sugar Cane - 011-1
    #    Other Traditional Exports - 011-2,011-3,011-4,011-5
    # Other Agricultural Crops
    #    Root Crops - 011-6
    #    Other Domestic Crops - 011-7,011-8
    # Animal Farming - 12
    # Post Harvest Crop Activities & Agricultural Services - 14
    # Forestry and logging - 20

    # Names in NIP_2023
    # - Vegetables & Condiments
    # - Root Crops & Other Tubers
    # - Fruits & Other Crops
    # - Animal Production
    # - Agricultural Services, Forestry & Fishing

    # Split forestry / other agricultural landuse
    agri_land_use["known_forest"] = (
        agri_land_use["subsector_code_forest"].astype(str).eq("20")
        & agri_land_use["subsector_code_tnc"].astype(str).eq("20")
    ).astype(int)
    agri_forest = agri_land_use[agri_land_use["known_forest"] == 1]
    agri_land_use = agri_land_use[agri_land_use["known_forest"] == 0]
    agri_land_use.drop(
        columns=["index", "index_left", "index_right"], inplace=True, errors="ignore"
    )

    # Intersect with crop yields
    agri_areas = gpd.sjoin(
        agri_land_use, crop_yields, how="inner", predicate="intersects"
    ).reset_index()
    agri_areas.rename(columns={"geometry": "agri_geometry"}, inplace=True)
    agri_areas = pd.merge(
        agri_areas, crop_yields[["crop_id", "geometry"]], how="left", on=["crop_id"]
    )
    agri_areas["geom"] = agri_areas.apply(
        lambda x: x["agri_geometry"].intersection(x["geometry"].buffer(0)), axis=1
    )
    agri_areas.drop(["agri_geometry", "geometry"], axis=1, inplace=True)
    agri_areas.rename(columns={"geom": "geometry"}, inplace=True)
    values_columns = (
        [
            "land_id",
            "forest_id",
            "Classify",
            "LU_CODE",
            "tnc_id",
            "NAME",
            "TNCCODE",
            "global_id",
            "global_LU_type",
            "sector_code_forest",
            "sector_code_tnc",
            "crop_tons_persqm",
        ]
        + list(set(all_sector_columns))
        + all_crop_columns
        + ["geometry"]
    )
    agri_areas = agri_areas[values_columns]
    agri_areas.set_crs(epsg=LOCAL_PROJ_CRS_EPSG, inplace=True)

    agri_areas["area_m2"] = agri_areas.apply(lambda x: x.geometry.area, axis=1)
    agri_areas["A_GDP"] = 0

    #
    # Join economic output subsector_code against agri_areas columns as codes
    #
    for _, econ in economic_output.iterrows():
        econ_subsector_code = econ["Subsector"]
        econ_code = f"A_{econ_subsector_code}_A"
        output = econ[value_colname]

        if econ_code in agri_areas.columns.values.tolist():
            logging.info("Weighting %s on area and crop production", econ_code)
            subsector_weight = agri_areas[econ_code] * agri_areas["area_m2"]
        else:
            logging.info("Weighting %s on area alone", econ_code)
            subsector_weight = agri_areas["area_m2"]

        agri_areas[f"GDP_{econ_subsector_code}"] = (
            output * subsector_weight / subsector_weight.sum()
        )

    # Estimate split between subsectors of "Agricultural Services, Forestry & Fishing"
    # Source: ./processed_data/macroeconomic_data/detailed_sector_GVA_GDP_current_prices.xlsx
    # 2014-2019 average split:
    # 0.16 post-harvest crop activities and agricultural services
    # 0.09 forestry and logging
    # 0.75 fishing
    service_forestry_fishing_output = economic_output.query(
        "Subsector == 'Agricultural Services, Forestry & Fishing'"
    )[value_colname].sum()
    service_output = service_forestry_fishing_output * 0.16
    forest_output = service_forestry_fishing_output * 0.09
    fishing_output = service_forestry_fishing_output * 0.75

    agri_areas["crop_tons"] = agri_areas["crop_tons_persqm"] * agri_areas["area_m2"]
    agri_areas["GDP_service"] = (
        service_output * agri_areas["crop_tons"] / agri_areas["crop_tons"].sum()
    )

    agri_forest["GDP_forest"] = (
        forest_output * agri_forest["area_m2"] / agri_forest["area_m2"].sum()
    )
    values_columns = [
        "land_id",
        "forest_id",
        "Classify",
        "LU_CODE",
        "tnc_id",
        "NAME",
        "TNCCODE",
        "global_id",
        "global_LU_type",
        "sector_code_forest",
        "sector_code_tnc",
        "area_m2",
        "geometry",
    ]

    agri_areas = pd.concat(
        [
            agri_areas[
                values_columns
                + [
                    "crop_tons",
                    "GDP_service",
                    "GDP_Vegetables & Condiments",
                    "GDP_Root Crops & Other Tubers",
                    "GDP_Fruits & Other Crops",
                    "GDP_Animal Production",
                ]
            ],
            agri_forest[values_columns + ["GDP_forest"]],
            non_agri_land_use[values_columns],
        ],
        axis=0,
        ignore_index=True,
    )
    agri_areas["crop_tons"] = agri_areas["crop_tons"].fillna(0)

    #
    # Fishing
    # - run after all-island areas are joined back, to link fishing locations
    #
    fishing_locations = gpd.read_file(fishing_locations_path, layer="areas")
    agri_areas["GDP_fishing"] = estimate_fishing_gdp(
        fishing_output, fishing_locations, agri_areas
    )

    #
    # Calculate totals and ratios
    #
    gdp_cols = [
        "GDP_forest",
        "GDP_service",
        "GDP_fishing",
        "GDP_Vegetables & Condiments",
        "GDP_Root Crops & Other Tubers",
        "GDP_Fruits & Other Crops",
        "GDP_Animal Production",
    ]
    for gdp_col in gdp_cols:
        agri_areas[gdp_col] = agri_areas[gdp_col].fillna(0)

    agri_areas["A_GDP"] = agri_areas[gdp_cols].sum(axis=1)

    agri_areas["GDP_persqm"] = agri_areas["A_GDP"] / agri_areas["area_m2"]
    agri_areas["crop_tons_persqm"] = agri_areas["crop_tons"] / agri_areas["area_m2"]
    agri_areas["GDP_unit"] = "GVA JMD Millions (2023)"

    agri_areas = gpd.GeoDataFrame(
        agri_areas, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    )
    agri_areas = remove_geometry_collections(agri_areas)

    tot_area = agri_areas["area_m2"].sum()
    logging.info("Given GDP %f", economic_output[value_colname].sum())
    logging.info("Estimated GDP %f", agri_areas["A_GDP"].sum())
    logging.info("Estimated Areas %f", tot_area)

    gpd.GeoDataFrame(
        agri_areas, geometry="geometry", crs=f"EPSG:{LOCAL_PROJ_CRS_EPSG}"
    ).to_file(
        output_areas,
        layer="areas",
        driver="GPKG",
    )


def estimate_fishing_gdp(fishing_output, fishing_locations, agri_areas):
    fishing_locations["farm_wt"] = (
        fishing_locations["Size_Farm"] / fishing_locations["Size_Farm"].sum()
    )
    # Join and maintain "index" as column
    fishing_areas = gpd.sjoin(
        agri_areas.reset_index(), fishing_locations, how="inner", predicate="intersects"
    )

    fishing_areas_total = (
        fishing_areas.groupby("farm_id")["area_m2"].sum().reset_index()
    )
    fishing_areas_total.rename(columns={"area_m2": "area_farms"}, inplace=True)
    fishing_areas = pd.merge(
        fishing_areas, fishing_areas_total, how="left", on=["farm_id"]
    )
    # Proportion GDP
    # - by farm according to normalised "Size_Farm"
    # - by landuse area within farm according to area
    fishing_areas["GDP_fishing"] = (
        fishing_output
        * fishing_areas.farm_wt
        * (fishing_areas["area_m2"] / fishing_areas["area_farms"])
    ).fillna(0)

    # Join back using fishing areas "index" maintained from agri_areas
    sector_df = agri_areas.reset_index().join(
        fishing_areas.groupby("index")[["GDP_fishing"]].sum(),
    )
    return sector_df.GDP_fishing


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    logging.info("Start agricultural production")
    main()
    logging.info("Done.")
