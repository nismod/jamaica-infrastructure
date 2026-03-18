"""Assign GDP values to buildings and aggregate to admin level"""

import logging
from collections import defaultdict

import click
import pandas as pd
import geopandas as gpd


def match_buildings_to_areas(buildings, gdf, building_id, gdf_ids):
    matches = gpd.sjoin(
        buildings, gdf, how="inner", predicate="intersects"
    ).reset_index()
    matches.rename(columns={"geometry": "building_geometry"}, inplace=True)
    matches = pd.merge(matches, gdf[gdf_ids + ["geometry"]], how="left", on=gdf_ids)
    matches["area_match"] = matches.apply(
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


def get_sector_gdp(output: float, buildings: gpd.GeoDataFrame, sector_code: str):
    if output <= 0:
        return 0

    # Filter for buildings in this sector
    sector_buildings = buildings[buildings[sector_code] == 1]

    # Aggregate to regions
    # sum over building footprint area
    # pick first (any) t_ij_ext (attractiveness measure)
    sector_regions = (
        sector_buildings.groupby(["ED_ID", "ED"])
        .agg({"t_ij_ext": "first", "area_sqm": "sum"})
        .rename(columns={"area_sqm": "region_total_area_sqm"})
    )

    # t_ij_ext is originally per-region
    # Sum t_ij_ext over all regions which contain this sector
    sector_tij_total: float = sector_regions.t_ij_ext.sum()

    # Join region_total_area_sqm back on to *all* buildings
    # columns: "ED_ID", "ED", "region_total_area_sqm"
    sector_region_area_sqm: pd.DataFrame = sector_regions[
        ["region_total_area_sqm"]
    ].reset_index()
    all_buildings = pd.merge(
        buildings, sector_region_area_sqm, how="left", on=["ED_ID", "ED"]
    ).fillna(0)

    # Calculate output for sector, in regions according to attractiveness, in buildings according to area
    sector_building_gdp = (
        output
        * (all_buildings[sector_code])
        * (all_buildings.t_ij_ext / sector_tij_total)
        * (all_buildings.area_sqm / all_buildings.region_total_area_sqm)
    ).fillna(0)

    assert (
        abs(sector_building_gdp.sum() - output) < 1
    ), f"Sense check total assigned sector GDP {sector_building_gdp.sum()} ~= amount to assign {output}"

    return sector_building_gdp


def filter_nonres_and_join_attractiveness(buildings, region_attractiveness):
    columns = [
        "osm_id",
        "ED_ID",
        "ED",
        "jic2016_sector",
        "building_type",
        "area_sqm",
        "geometry",
    ]
    nonres_buildings = (
        buildings[columns]
        .query("building_type != 'Residential'")
        .dropna(subset="jic2016_sector")
        .copy()
        .merge(region_attractiveness[["ED_ID", "ED", "t_ij_ext"]], on=["ED_ID", "ED"])
    )
    return nonres_buildings


def add_sector_indicator_columns(nonres_buildings, sector_codes):
    """Add a column per sector code, value 0 or 1 indicating if the building is
    tagged with that sector"""
    for sector_code in sector_codes:
        nonres_buildings[sector_code] = nonres_buildings.jic2016_sector.str.contains(
            sector_code
        ).astype(int)

    return nonres_buildings


def allocate_mining_buildings(mining_areas, sector_df):
    sector_columns = [c for c in sector_df.columns.values.tolist() if "C_" in c[:2]]
    sector_df["C"] = sector_df[sector_columns].sum(axis=1)
    sector_df["C"] = sector_df.apply(lambda x: 1 if x["C"] > 0 else 0, axis=1)
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
    mining_buildings["area_sqm"] = mining_buildings.apply(
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
    mining_areas_gdp["GDP_building_ratio"] = mining_areas_gdp.apply(
        lambda x: x["C_GDP_building"] / x["C_GDP"] if x["C_GDP"] > 0 else 0, axis=1
    )
    return mining_areas_gdp, mining_buildings


@click.command()
@click.version_option("1.0")
@click.option(
    "--region_attractiveness_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to region_attractiveness Geoparquet",
)
@click.option(
    "--buildings_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to buildings Geoparquet",
)
@click.option(
    "--economic_output_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to economic output Excel",
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
    "--output",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    help="Path to buildings GeoPackage output",
)
def main(
    region_attractiveness_path,
    buildings_path,
    economic_output_path,
    agriculture_buildings_path,
    mining_areas_path,
    output,
):
    """
    Allocate national, sectoral GDP to buildings.

    Example usage:

        python workflow/1_context/spatial_economic_allocation.py \
            --region_attractiveness_path processed_data/population/region_attractiveness.geoparquet \
            --buildings_path processed_data/buildings//buildings_nic2016.geoparquet \
            --economic_output_path processed_data/macroeconomic_data/NIP_2023.csv \
            --agriculture_buildings_path processed_data/agriculture_data/building_agricuture_gdp.csv \
            --mining_areas_path processed_data/mining_data/mining_gdp.gpkg \
            --output processed_data/buildings/buildings_assigned_economic_activity.geoparquet
    """
    ANNUAL_TO_DAILY = 1 / 365
    MILLIONS_TO_JMD = 1e6

    region_attractiveness = gpd.read_parquet(region_attractiveness_path)

    buildings = gpd.read_parquet(buildings_path)
    logging.info(
        "Read buildings: all %s",
        buildings.shape,
    )

    sector_codes = sorted(
        list(
            pd.Series(buildings.jic2016_sector.dropna().unique())
            .apply(lambda s: s.split(","))
            .explode()
            .unique()
        )
    )

    nonres_buildings = add_sector_indicator_columns(
        filter_nonres_and_join_attractiveness(buildings, region_attractiveness),
        sector_codes,
    )
    logging.info(
        "Joined buildings: non-residential %s all %s",
        nonres_buildings.shape,
        buildings.shape,
    )

    economic_output = pd.read_csv(economic_output_path, comment="#")
    assert (
        economic_output.columns
        == [
            "Code",
            "Sector",
            "Subsector",
            "GVA JMD Millions (2023)",
        ]
    ).all(), 'Expected columns "Code","Sector","Subsector","GVA JMD Millions (2023)" in NIP CSV'

    for sector_code in sector_codes:
        sector = economic_output.query(f"Code == '{sector_code}'")
        sector_output_per_day = (
            sector["GVA JMD Millions (2023)"].sum() * ANNUAL_TO_DAILY * MILLIONS_TO_JMD
        )
        sector_description = sector.Sector.iloc[0]
        subsector_description = "; ".join(list(sector.Subsector))
        logging.info(
            "Assigning GDP/day (JMD) for sector %s : %f (%s: %s)",
            sector_code,
            sector_output_per_day,
            sector_description,
            subsector_description,
        )
        if sector_code == "A":
            nonres_buildings[f"{sector_code}_GDP"] = (
                nonres_buildings[sector_code]
                * sector_output_per_day
                / nonres_buildings[sector_code].sum()
            )
            # # TODO figure out generation of agriculture buildings
            # agri_buildings = pd.read_csv(agriculture_buildings_path)
            # agri_buildings["osm_id"] = agri_buildings["osm_id"].astype(int)
            # nonres_buildings["osm_id"] = nonres_buildings["osm_id"].astype(int)
            # nonres_buildings = pd.merge(
            #     nonres_buildings, agri_buildings, how="left", on=["osm_id"]
            # )
            # nonres_buildings[["osm_id", "A_GDP_building"]].to_csv("test.csv")
            # logging.debug("GDP to assign %f", nonres_buildings["A_GDP_building"].sum())
            # nonres_buildings["A_GDP_building"] = nonres_buildings[
            #     "A_GDP_building"
            # ].fillna(0)
            # nonres_buildings[f"{sector_code}_GDP"] += nonres_buildings["A_GDP_building"]
            # nonres_buildings.drop("A_GDP_building", axis=1, inplace=True)
            # del agri_buildings

        elif sector_code == "B":
            nonres_buildings[f"{sector_code}_GDP"] = (
                nonres_buildings[sector_code]
                * sector_output_per_day
                / nonres_buildings[sector_code].sum()
            )
            # mining_areas = gpd.read_file(
            #     mining_areas_path,
            #     layer="areas",
            # )
            # _, mining_buildings = allocate_mining_buildings(
            #     mining_areas, nonres_buildings.copy()
            # )
            # # TODO check if we should write mining_areas_gdp to processed_data/mining_data/mining_gdp.gpkg
            # # if C_GDP_building and GDP_building_ratio are used later?

            # mining_buildings["osm_id"] = mining_buildings["osm_id"].astype(int)
            # nonres_buildings["osm_id"] = nonres_buildings["osm_id"].astype(int)
            # nonres_buildings = pd.merge(
            #     nonres_buildings, mining_buildings, how="left", on=["osm_id"]
            # )
            # nonres_buildings[["osm_id", "C_GDP_building"]].to_csv("test.csv")
            # logging.debug("GDP to assign %f", nonres_buildings["C_GDP_building"].sum())
            # nonres_buildings["C_GDP_building"] = nonres_buildings[
            #     "C_GDP_building"
            # ].fillna(0)
            # nonres_buildings[f"{sector_code}_GDP"] += nonres_buildings["C_GDP_building"]
            # nonres_buildings.drop("C_GDP_building", axis=1, inplace=True)

        else:
            nonres_buildings[f"{sector_code}_GDP"] = get_sector_gdp(
                sector_output_per_day, nonres_buildings, sector_code
            )

    gdp_columns = [f"{sector_code}_GDP" for sector_code in sector_codes]

    buildings = pd.merge(
        buildings, nonres_buildings[gdp_columns + ["osm_id"]], how="left", on=["osm_id"]
    )
    buildings["total_GDP"] = buildings[gdp_columns].sum(axis=1)
    buildings["GDP_unit"] = "JD/day"

    buildings.to_parquet(output)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    main()
