"""Assign GDP values to buildings and aggregate to admin level"""

import logging

import click

import pandas as pd
import geopandas as gpd
from collections import defaultdict

from jamaica_infrastructure.geo import LOCAL_PROJ_CRS_EPSG


def match_buildings_to_areas(buildings, gdf, building_id, gdf_ids):
    matches = gpd.sjoin(
        buildings, gdf, how="inner", predicate="intersects"
    ).reset_index()
    matches.rename(columns={"geometry": "building_geometry"}, inplace=True)
    matches = pd.merge(matches, gdf[gdf_ids + ["geometry"]], how="left", on=gdf_ids)
    matches["area_match"] = matches.progress_apply(
        lambda x: x["building_geometry"].intersection(x["geometry"].buffer(0)).area,
        axis=1,
    )
    matches = matches.sort_values(by=["area_match"], ascending=False)
    matches = matches.drop_duplicates(subset=[building_id], keep="first")
    matches.drop(["area_match", "geometry"], axis=1, inplace=True)
    matches.rename(columns={"building_geometry": "geometry"}, inplace=True)

    return matches


def get_nearest_areas(x, gdf, gdf_column):
    area_index = gdf.distance(x.geometry).sort_values().index[0]
    return gdf.loc[area_index, gdf_column]


def get_sector_gdp(
    economic_output,
    sector_df,
    sector_columns,
    sector_code,
    subsector_code,
    financial_year,
    tax_rate,
):
    output = (
        (1 + tax_rate)
        * (1.0e6 / 365.0)
        * economic_output[
            (economic_output["sector_code"] == sector_code)
            & (economic_output["subsector_code"] == subsector_code)
        ][f"{financial_year}"].sum()
    )
    print("* Given GDP", output)
    if output > 0:
        sector_df[sector_code] = sector_df[sector_columns].sum(axis=1)
        sector_df[sector_code] = sector_df.progress_apply(
            lambda x: 1 if x[sector_code] > 0 else 0, axis=1
        )
        sector_df[f"{sector_code}_t_ij_ext"] = (
            sector_df["t_ij_ext"] * sector_df[sector_code]
        )

        pop_sums = (
            sector_df[sector_df[sector_code] == 1]
            .drop_duplicates(subset=["ED_ID", "ED"], keep="first")[
                f"{sector_code}_t_ij_ext"
            ]
            .sum()
        )
        sector_df[f"{sector_code}_t_ij_ext"] = (
            sector_df[f"{sector_code}_t_ij_ext"] / pop_sums
        )

        area_sums = sector_df[["ED_ID", "ED", "area_sqm", sector_code]]
        area_sums["total_areas"] = area_sums["area_sqm"] * area_sums[sector_code]
        area_sums = area_sums.groupby(["ED_ID", "ED"])["total_areas"].sum()
        sector_df = pd.merge(
            sector_df, area_sums, how="left", on=["ED_ID", "ED"]
        ).fillna(0)
        sector_df["assigned_GDP"] = sector_df.progress_apply(
            lambda x: (
                output
                * x[f"{sector_code}_t_ij_ext"]
                * x["area_sqm"]
                * x[sector_code]
                / x["total_areas"]
                if x["total_areas"] > 0
                else 0
            ),
            axis=1,
        )
        sector_df[f"{sector_code}_GDP"] += sector_df["assigned_GDP"]
        sector_df.drop(["assigned_GDP", "total_areas"], axis=1, inplace=True)
    return sector_df


def get_fishery_gdp(
    economic_output,
    sector_df,
    fishing_locations,
    sector_weight,
    sector_code,
    subsector_code,
    financial_year,
    tax_rate,
):
    output = (
        (1 + tax_rate)
        * (1.0e6 / 365.0)
        * economic_output[
            (economic_output["sector_code"] == sector_code)
            & (economic_output["subsector_code"] == subsector_code)
        ][f"{financial_year}"].sum()
    )
    print("* Given GDP", output)
    if output > 0:
        fishing_areas = gpd.sjoin(
            fishing_locations, sector_df, how="inner", predicate="intersects"
        ).reset_index()
        fishing_areas_total = (
            fishing_areas.groupby("farm_id")["area_sqm"].sum().reset_index()
        )
        fishing_areas_total.rename(columns={"area_sqm": "area_farms"}, inplace=True)
        fishing_areas = pd.merge(
            fishing_areas, fishing_areas_total, how="left", on=["farm_id"]
        )
        fishing_areas["fishing_GDP"] = (
            output
            * fishing_areas[sector_weight]
            * (fishing_areas["area_sqm"] / fishing_areas["area_farms"])
        )
        sector_df = pd.merge(
            sector_df,
            fishing_areas[["osm_id", "fishing_GDP"]],
            how="left",
            on=["osm_id"],
        ).fillna(0)
        sector_df[f"{sector_code}_GDP"] += sector_df["fishing_GDP"]

    return sector_df


def region_attractiveness_by_population(
    population_areas, buffer_distance, population_year, epsg_jamaica
):

    population = population_areas.copy()

    average_pop = (
        population_areas.groupby(["PARISH"])[f"{population_year}"].mean().reset_index()
    )
    average_pop.columns = average_pop.columns.map(str)
    for i, avg in average_pop.iterrows():
        population.loc[
            (population[f"{population_year}"] == 0)
            & (population["PARISH"] == avg["PARISH"]),
            f"{population_year}",
        ] = avg[f"{population_year}"]

    population["geometry"] = population.apply(lambda x: x.geometry.centroid, axis=1)

    """Find the points within a distance buffer
    """
    population_buffer = population.copy()
    population_buffer["geometry"] = population_buffer.geometry.buffer(buffer_distance)
    population_buffer.rename(
        columns={
            "ED_ID": "from_ED_ID",
            "ED": "from_ED",
            f"{population_year}": f"from_{population_year}",
            f"working_{population_year}": f"from_working_{population_year}",
        },
        inplace=True,
    )
    within_distance = gpd.sjoin(
        population,
        population_buffer[
            [
                "from_ED_ID",
                "from_ED",
                f"from_{population_year}",
                f"from_working_{population_year}",
                "geometry",
            ]
        ],
        how="inner",
        predicate="within",
    )
    total_within = (
        within_distance.groupby(["from_ED_ID", "from_ED"])[
            f"{population_year}", f"working_{population_year}"
        ]
        .sum()
        .reset_index()
    )
    total_within.rename(
        columns={
            f"{population_year}": f"total_{population_year}",
            f"working_{population_year}": f"total_working_{population_year}",
        },
        inplace=True,
    )

    within_distance = pd.merge(
        within_distance[
            [
                "ED_ID",
                "ED",
                f"{population_year}",
                f"working_{population_year}",
                "from_ED_ID",
                "from_ED",
                f"from_{population_year}",
                f"from_working_{population_year}",
            ]
        ],
        total_within,
        how="left",
        on=["from_ED_ID", "from_ED"],
    )

    within_distance["t_ij_ext"] = within_distance.progress_apply(
        lambda x: (
            x[f"working_{population_year}"]
            / (x[f"total_working_{population_year}"] - x[f"working_{population_year}"])
        ),
        axis=1,
    )
    within_distance_sums = (
        within_distance.groupby(["from_ED_ID", "from_ED"])["t_ij_ext"]
        .sum()
        .reset_index()
    )
    within_distance_sums.rename(columns={"t_ij_ext": "t_ij_ext_sums"}, inplace=True)
    within_distance = pd.merge(
        within_distance, within_distance_sums, how="left", on=["from_ED_ID", "from_ED"]
    ).fillna(0)
    del within_distance_sums
    within_distance["t_ij_ext"] = within_distance.progress_apply(
        lambda x: x[f"from_{population_year}"] * (x["t_ij_ext"] / x["t_ij_ext_sums"]),
        axis=1,
    )
    within_distance.drop("t_ij_ext_sums", axis=1, inplace=True)
    region_attractiveness = (
        within_distance.groupby(["ED_ID", "ED"])["t_ij_ext"].sum().reset_index()
    )

    region_attractiveness = pd.merge(
        region_attractiveness,
        population_areas[["ED_ID", "ED", "geometry"]],
        how="left",
        on=["ED_ID", "ED"],
    )
    region_attractiveness = gpd.GeoDataFrame(
        region_attractiveness, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
    )
    return region_attractiveness


def join_buildings(buildings, region_attractiveness, population_areas, epsg_jamaica):

    new_columns = ["ED_ID", "ED", "PARISH", "CONST_NAME"]
    for col in new_columns:
        if col in buildings.columns.values.tolist():
            buildings.drop(col, axis=1, inplace=True)

    buildings_regions = match_buildings_to_areas(
        buildings, region_attractiveness, "osm_id", ["ED_ID", "ED"]
    )

    buildings_nomatches = buildings[
        ~(buildings["osm_id"].isin(buildings_regions["osm_id"].values.tolist()))
    ]
    buildings_nomatches["ED_ID"] = buildings_nomatches.progress_apply(
        lambda x: get_nearest_areas(x, region_attractiveness, "ED_ID"), axis=1
    )
    buildings_nomatches["ED"] = buildings_nomatches.progress_apply(
        lambda x: get_nearest_areas(x, region_attractiveness, "ED"), axis=1
    )
    buildings_nomatches = pd.merge(
        buildings_nomatches,
        region_attractiveness[["ED_ID", "ED", "t_ij_ext"]],
        how="left",
        on=["ED_ID", "ED"],
    )

    columns = [
        "osm_id",
        "ED_ID",
        "ED",
        "sector_code",
        "subsector_code",
        "assigned_attribute",
        "building_type",
        "area_sqm",
        "t_ij_ext",
    ]

    buildings_regions = pd.concat(
        [buildings_regions[columns], buildings_nomatches[columns]],
        axis=0,
        ignore_index=True,
    )
    buildings_regions = pd.merge(
        buildings_regions,
        population_areas[["ED_ID", "ED", "PARISH", "CONST_NAME"]],
        how="left",
        on=["ED_ID", "ED"],
    )
    buildings = pd.merge(
        buildings,
        buildings_regions[["osm_id", "ED_ID", "ED", "PARISH", "CONST_NAME"]],
        how="left",
        on=["osm_id"],
    )

    commercial_buildings = buildings_regions[
        buildings_regions["building_type"] != "Residential"
    ]
    sector_df = commercial_buildings.copy()
    sector_dict = defaultdict(list)
    for commerical in commercial_buildings.itertuples():
        sectors = list(
            set(
                zip(
                    str(commerical.sector_code).split(","),
                    str(commerical.subsector_code).split(","),
                )
            )
        )
        sectors = [f"{s[0]}_{s[1]}" for s in sectors]
        sectors = [s for s in sectors if s != "RES_RES"]
        osm_id = commerical.osm_id
        for s in sectors:
            sector_dict[s].append(osm_id)

    for k, v in sector_dict.items():
        df = pd.DataFrame(v, columns=["osm_id"])
        df[k] = 1
        sector_df = pd.merge(sector_df, df, how="left", on=["osm_id"]).fillna(0)

    sector_df = pd.merge(
        sector_df, buildings[["osm_id", "geometry"]], how="left", on=["osm_id"]
    )
    sector_df = gpd.GeoDataFrame(
        sector_df, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
    )
    return buildings, sector_df


def aggregate_to_admin_gdp(buildings, population_areas, epsg_jamaica):
    sector_codes = [
        "A",
        "B",
        "C",
        "D",
        "E",
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
    gdp_columns = [f"{scode}_GDP" for scode in sector_codes] + ["total_GDP"]

    admin_gdp = (
        buildings.groupby(["ED_ID", "ED", "PARISH", "CONST_NAME"])[gdp_columns]
        .sum()
        .reset_index()
    )
    admin_gdp["GDP_unit"] = "JD/day"
    admin_gdp = gpd.GeoDataFrame(
        pd.merge(
            admin_gdp,
            population_areas[["ED_ID", "ED", "geometry"]],
            how="left",
            on=["ED_ID", "ED"],
        ),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    return admin_gdp


def allocate_mining_buildings(mining_areas, sector_df):
    sector_columns = [c for c in sector_df.columns.values.tolist() if "C_" in c[:2]]
    sector_df["C"] = sector_df[sector_columns].sum(axis=1)
    sector_df["C"] = sector_df.progress_apply(lambda x: 1 if x["C"] > 0 else 0, axis=1)
    mining_buildings = gpd.sjoin(
        mining_areas[["mining_id", "GDP_persqm", "geometry"]],
        sector_df[sector_df["C"] == 1][["osm_id", "geometry"]],
        how="inner",
        predicate="intersects",
    ).reset_index()

    mining_buildings.rename(columns={"geometry": "landuse_geometry"}, inplace=True)
    mining_buildings = pd.merge(
        mining_buildings, sector_df[["osm_id", "geometry"]], how="left", on=["osm_id"]
    )
    mining_buildings["area_sqm"] = mining_buildings.progress_apply(
        lambda x: x["landuse_geometry"].intersection(x["geometry"].buffer(0)).area,
        axis=1,
    )
    mining_buildings["C_GDP_building"] = (
        mining_buildings["GDP_persqm"] * mining_buildings["area_sqm"]
    )
    mining_buildings_areas = (
        mining_buildings.groupby(["mining_id"])["C_GDP_building"].sum().reset_index()
    )
    mining_buildings = (
        mining_buildings.groupby(["osm_id"])["C_GDP_building"].sum().reset_index()
    )

    mining_areas_gdp = pd.merge(
        mining_areas, mining_buildings_areas, how="left", on=["mining_id"]
    )
    mining_areas_gdp["C_GDP_building"] = mining_areas_gdp["C_GDP_building"].fillna(0)
    mining_areas_gdp["GDP_building_ratio"] = mining_areas_gdp.progress_apply(
        lambda x: x["C_GDP_building"] / x["C_GDP"] if x["C_GDP"] > 0 else 0, axis=1
    )
    return mining_areas_gdp, mining_buildings


@click.command()
@click.version_option("1.0")
@click.option(
    "--population_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to population GeoPackage",
)
@click.option(
    "--buildings_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to buildings GeoPackage",
)
@click.option(
    "--economic_output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to economic output Excel",
)
@click.option(
    "--fishing_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to fishing/aqua_farms GeoPackage",
)
@click.option(
    "--agriculture_buildings_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to agriculature buildings CSV",
)
@click.option(
    "--mining_areas_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to mining areas GeoPackage",
)
@click.option(
    "--output_buildings",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    help="Path to buildings GeoPackage output",
)
@click.option(
    "--output_admin",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    help="Path to admin GeoPackage output",
)
def main(
    population_path,
    buildings_path,
    economic_output_path,
    fishing_path,
    agriculture_buildings_path,
    mining_areas_path,
    output_buildings,
    output_admin,
):
    """
    --population processed_data/population/population_projections.gpkg \
    --buildings processed_data/buildings/buildings_assigned_economic_activity.gpkg \
    --economic_output_path processed_data/macroeconomic_data/detailed_sector_GVA_GDP_current_prices.xlsx \
    --fishing_path processed_data/land_type_and_use/aqua_farms.gpkg \
    --agriculture_buildings_path processed_data/agriculture_data/building_agricuture_gdp.csv \
    --mining_areas_path processed_data/mining_data/mining_gdp.gpkg \
    --output_buildings processed_data/buildings/buildings_assigned_economic_activity.gpkg \
    --output_admin processed_data/buildings/admin_level_assigned_economic_activity.gpkg
    """
    epsg_jamaica = LOCAL_PROJ_CRS_EPSG
    financial_year = 2019
    buffer_distance = 1e4  # 10 km distance buffer
    population_year = 2019

    population_areas = gpd.read_file(
        population_path,
        layer="mean",
    ).to_crs(epsg=epsg_jamaica)

    region_attractiveness = region_attractiveness_by_population(
        population_areas, buffer_distance, population_year, epsg_jamaica
    )

    buildings = gpd.read_file(
        buildings_path,
        layer="areas",
    )
    sector_df, buildings = join_buildings(
        buildings, region_attractiveness, population_areas, epsg_jamaica
    )

    economic_output = pd.read_excel(
        economic_output_path,
        sheet_name="2019",
    )
    economic_output.columns = [
        str(c).strip() for c in economic_output.columns.values.tolist()
    ]
    economic_output["subsector_code"] = economic_output["subsector_code"].apply(str)
    total_gva = economic_output[economic_output["sector_code"] == "GVA"][
        f"{financial_year}"
    ].sum()
    total_tax = economic_output[economic_output["sector_code"] == "TAX"][
        f"{financial_year}"
    ].sum()
    tax_rate = 1.0 * total_tax / total_gva

    sector_df = sector_df.to_crs(epsg=epsg_jamaica)
    sector_codes = [
        "A",
        "B",
        "C",
        "D",
        "E",
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
    gdp_columns = [f"{scode}_GDP" for scode in sector_codes]
    for scode in sector_codes:
        logging.info("* Start sector", scode)
        sector_columns = [
            c for c in sector_df.columns.values.tolist() if f"{scode}_" in c[:2]
        ]
        sector_df[f"{scode}_GDP"] = 0
        if scode in ("D", "H", "J", "K", "L"):
            logging.debug(sector_columns)
            sector_df = get_sector_gdp(
                economic_output,
                sector_df,
                sector_columns,
                scode,
                "ALL",
                financial_year,
                tax_rate,
            )
        elif scode in ("F", "G", "M", "N", "O"):
            for sc in sector_columns:
                sector_code = sc.split("_")[0]
                subsector_code = sc.split("_")[1]
                logging.debug(sector_code, subsector_code)
                sector_df = get_sector_gdp(
                    economic_output,
                    sector_df,
                    [sc],
                    sector_code,
                    subsector_code,
                    financial_year,
                    tax_rate,
                )
        elif scode == "B":
            fishing_locations = gpd.read_file(
                fishing_path,
                layer="areas",
            ).to_crs(epsg=epsg_jamaica)
            fishing_locations["farm_wt"] = (
                fishing_locations["Size_Farm"] / fishing_locations["Size_Farm"].sum()
            )
            sector_df = get_fishery_gdp(
                economic_output,
                sector_df,
                fishing_locations,
                "farm_wt",
                "B",
                "50",
                financial_year,
                tax_rate,
            )
            del fishing_locations

        elif scode == "A":
            agri_buildings = pd.read_csv(agriculture_buildings_path)
            agri_buildings["osm_id"] = agri_buildings["osm_id"].astype(int)
            sector_df["osm_id"] = sector_df["osm_id"].astype(int)
            sector_df = pd.merge(sector_df, agri_buildings, how="left", on=["osm_id"])
            sector_df[["osm_id", "A_GDP_building"]].to_csv("test.csv")
            logging.debug("GDP to assign", sector_df["A_GDP_building"].sum())
            sector_df["A_GDP_building"] = sector_df["A_GDP_building"].fillna(0)
            sector_df[f"{scode}_GDP"] += sector_df["A_GDP_building"]
            sector_df.drop("A_GDP_building", axis=1, inplace=True)
            del agri_buildings

        elif scode == "C":
            mining_areas = gpd.read_file(
                mining_areas_path,
                layer="areas",
            )
            mining_areas_gdp, mining_buildings = allocate_mining_buildings(
                mining_areas, sector_df.copy()
            )
            # TODO check if we should write mining_areas_gdp to processed_data/mining_data/mining_gdp.gpkg
            # if C_GDP_building and GDP_building_ratio are used later?

            mining_buildings["osm_id"] = mining_buildings["osm_id"].astype(int)
            sector_df["osm_id"] = sector_df["osm_id"].astype(int)
            sector_df = pd.merge(sector_df, mining_buildings, how="left", on=["osm_id"])
            sector_df[["osm_id", "C_GDP_building"]].to_csv("test.csv")
            logging.debug("GDP to assign", sector_df["C_GDP_building"].sum())
            sector_df["C_GDP_building"] = sector_df["C_GDP_building"].fillna(0)
            sector_df[f"{scode}_GDP"] += sector_df["C_GDP_building"]
            sector_df.drop("C_GDP_building", axis=1, inplace=True)
            del mining_buildings

        elif scode == "I":
            post_office_share = 0.157
            sector_df["find_post"] = sector_df.progress_apply(
                lambda x: (
                    1
                    if "post_office" in str(x["assigned_attribute"])
                    or "postal_offices" in str(x["assigned_attribute"])
                    else 0
                ),
                axis=1,
            )
            sector_df["I_640"] = sector_df["I_640"] * sector_df["find_post"]
            sector_df = get_sector_gdp(
                economic_output,
                sector_df,
                ["I_640"],
                scode,
                "640",
                financial_year,
                tax_rate,
            )
            sector_df[f"{scode}_GDP"] = post_office_share * sector_df[f"{scode}_GDP"]
            sector_df.drop("find_post", axis=1, inplace=True)
            sector_df = sector_df[sector_df["sector_code"] != "I"]
        elif scode == "E":
            sector_df = sector_df[~sector_df["sector_code"].isin(["E", "E,I", "I,E"])]

        logging.info("* Estimated GDP", sector_df[f"{scode}_GDP"].sum())
        logging.info("* Done with sector", scode)

    sector_df = sector_df[["osm_id"] + gdp_columns]
    sector_df["total_GDP"] = sector_df[gdp_columns].sum(axis=1)

    for col in gdp_columns + ["total_GDP"]:
        if col in buildings.columns.values.tolist():
            buildings.drop(col, axis=1, inplace=True)
    # Annoying fix to make sure they merge
    sector_df["osm_id"] = sector_df["osm_id"].astype(int)
    buildings["osm_id"] = buildings["osm_id"].astype(int)

    buildings = pd.merge(buildings, sector_df, how="left", on=["osm_id"])
    for col in gdp_columns + ["total_GDP"]:
        buildings[col] = buildings[col].fillna(0)
    buildings["GDP_unit"] = "JD/day"
    buildings = gpd.GeoDataFrame(
        buildings, geometry="geometry", crs=f"EPSG:{epsg_jamaica}"
    )
    buildings.to_file(
        output_buildings,
        layer="areas",
        driver="GPKG",
    )

    admin_gdp = aggregate_to_admin_gdp(buildings, population_areas, epsg_jamaica)
    admin_gdp.to_file(
        output_admin,
        layer="areas",
        driver="GPKG",
    )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    main()
