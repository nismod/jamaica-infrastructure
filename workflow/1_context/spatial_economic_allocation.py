"""Assign GDP values to buildings and aggregate to admin level"""

import logging
from collections import defaultdict

import click
import pandas as pd
import geopandas as gpd


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


def allocate_agriculture_to_buildings(buildings, agri_areas):

    sector_columns = [c for c in buildings.columns.values.tolist() if "A_" in c[:2]]
    buildings["A"] = buildings[sector_columns].sum(axis=1)
    buildings["A"] = buildings.progress_apply(lambda x: 1 if x["A"] > 0 else 0, axis=1)
    agri_buildings = gpd.sjoin(
        agri_areas[["land_id", "GDP_persqm", "geometry"]],
        buildings[buildings["A"] == 1][["osm_id", "geometry"]],
        how="inner",
        predicate="intersects",
    ).reset_index()

    agri_buildings.rename(columns={"geometry": "landuse_geometry"}, inplace=True)
    agri_buildings = pd.merge(
        agri_buildings, buildings[["osm_id", "geometry"]], how="left", on=["osm_id"]
    )
    agri_buildings["area_sqm"] = agri_buildings.progress_apply(
        lambda x: x["landuse_geometry"].intersection(x["geometry"].buffer(0)).area,
        axis=1,
    )
    agri_buildings["A_GDP_building"] = (
        agri_buildings["GDP_persqm"] * agri_buildings["area_sqm"]
    )
    agri_buildings = (
        agri_buildings.groupby(["osm_id"])["A_GDP_building"].sum().reset_index()
    )
    return agri_buildings


def allocate_mining_to_buildings(
    sector_output_per_day: float,
    mining_areas: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame,
    sector_code: str,
):
    total_mining_gdp = mining_areas["mining_gdp"].sum() * 1e6 / 365

    assert (
        abs(total_mining_gdp - sector_output_per_day) < 1
    ), f"Sense check total sector GDP {sector_output_per_day} ~= mining GDP assigned to areas {total_mining_gdp}"

    # join buildings intersecting with mining areas
    mining_buildings = (
        gpd.sjoin(
            mining_areas[["mining_id", "mining_gdp", "GDP_persqm", "geometry"]],
            buildings[buildings[sector_code] == 1][["osm_id", "geometry"]],
            how="inner",
            predicate="intersects",
        )
        .reset_index()
        .rename(columns={"geometry": "landuse_geometry"})
    )
    # join building osm_id and building geometry
    mining_buildings = pd.merge(
        mining_buildings, buildings[["osm_id", "geometry"]], how="left", on=["osm_id"]
    )
    # calculate mining footprint area
    mining_buildings["area_sqm"] = mining_buildings.apply(
        lambda x: x["landuse_geometry"].intersection(x["geometry"].buffer(0)).area,
        axis=1,
    )
    # allocate **part** of mining GDP from mining area to buildings
    # using building area and productivity per m2 of the mining area
    # NB this assigns some value to buildings, but the mining areas themselves remain
    # the source of economic activity for this sector.
    mining_buildings[f"{sector_code}_GDP"] = (
        mining_buildings["GDP_persqm"] * mining_buildings["area_sqm"]
    )

    mining_buildings = (
        mining_buildings.groupby(["osm_id"])[f"{sector_code}_GDP"].sum().reset_index()
    )

    # merge back to all buildings
    all_buildings = pd.merge(
        buildings,
        mining_buildings[["osm_id", f"{sector_code}_GDP"]],
        how="left",
        on=["osm_id"],
    ).fillna(0)

    assigned_to_buildings = mining_buildings[f"{sector_code}_GDP"].sum()
    logging.info(
        "Assigned %f fraction of sector GDP to buildings (%f of %f)",
        assigned_to_buildings / total_mining_gdp,
        assigned_to_buildings,
        total_mining_gdp,
    )

    return all_buildings[f"{sector_code}_GDP"]


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
    "--agriculture_areas_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to agriculture areas GeoPackage",
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
    agriculture_areas_path,
    mining_areas_path,
    output,
):
    """
    Allocate national, sectoral GDP to buildings.

    Example usage:

        python workflow/1_context/spatial_economic_allocation.py \
            --region_attractiveness_path processed_data/population/region_attractiveness.geoparquet \
            --buildings_path processed_data/buildings/buildings_nic2016.geoparquet \
            --economic_output_path processed_data/macroeconomic_data/NIP_2023.csv \
            --agriculture_buildings_path processed_data/agriculture_data/building_agricuture_gdp.csv \
            --mining_areas_path processed_data/mining_data/mining_gdp.gpkg \
            --output processed_data/buildings/buildings_assigned_economic_activity.geoparquet
    """
    ANNUAL_TO_DAILY = 1 / 365
    MILLIONS_TO_UNIT = 1e6

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

    # Read national industrial product (sectoral GVA)
    economic_output = pd.read_csv(economic_output_path, comment="#")

    # Check we read the expected columns
    assert (
        economic_output.columns
        == [
            "Code",
            "Sector",
            "Subsector",
            "GVA JMD Millions (2023)",
        ]
    ).all(), 'Expected columns "Code","Sector","Subsector","GVA JMD Millions (2023)" in NIP CSV'

    # Check we can assign all sectors to some buildings
    economic_output_codes = economic_output.Code
    msg = f"Expected economic output CSV to contain same set of codes as buildings, {set(economic_output_codes)} != {set(sector_codes)}"
    assert set(economic_output_codes) == set(sector_codes), msg

    # TODO REMOVE
    sector_codes = ["A"]
    # TODO REMOVE

    for sector_code in sector_codes:
        sector = economic_output.query(f"Code == '{sector_code}'")
        sector_output_per_day = (
            sector["GVA JMD Millions (2023)"].sum() * ANNUAL_TO_DAILY * MILLIONS_TO_UNIT
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

        #
        # Agriculture
        #
        if sector_code == "A":
            agriculture_areas = gpd.read_file(
                agriculture_areas_path,
                layer="areas",
            )

            nonres_buildings[f"{sector_code}_GDP"] = allocate_agriculture_to_buildings(
                sector_output_per_day, agriculture_areas, nonres_buildings, sector_code
            )

        #
        # Mining
        # - areas remain the representation of all mining GDP, this just assigns
        #   some to buildings in those areas
        #
        elif sector_code == "B":
            mining_areas = gpd.read_file(
                mining_areas_path,
                layer="areas",
            )
            nonres_buildings[f"{sector_code}_GDP"] = allocate_mining_to_buildings(
                sector_output_per_day, mining_areas, nonres_buildings, sector_code
            )

        #
        # All other sectors
        #
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
