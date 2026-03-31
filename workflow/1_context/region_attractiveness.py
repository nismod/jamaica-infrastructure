"""Assign 't_ij_ext' attractiveness parameter based on proximity of populations
to each region"""

import logging

import click
import geopandas as gpd
import pandas as pd
from jamaica_infrastructure.geo import LOCAL_PROJ_CRS_EPSG


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
            [f"{population_year}", f"working_{population_year}"]
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

    within_distance["t_ij_ext"] = within_distance.apply(
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
    within_distance["t_ij_ext"] = within_distance.apply(
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


@click.command()
@click.version_option("1.0")
@click.option(
    "--population_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to population GeoPackage",
)
@click.option(
    "--output_region_attractiveness_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, writable=True),
    help="Path to population GeoPackage",
)
def main(population_path, output_region_attractiveness_path):
    """Assign 't_ij_ext' attractiveness parameter based on proximity of populations to a given region

    Example usage:

        python workflow/1_context/region_attractiveness.py \
            --population_path processed_data/population/population_projections.gpkg \
            --output_region_attractiveness_path processed_data/population/region_attractiveness.geoparquet
    """
    buffer_distance = 1e4  # 10 km distance buffer
    population_year = 2019

    population_areas = gpd.read_file(
        population_path,
        layer="mean",
    ).to_crs(epsg=LOCAL_PROJ_CRS_EPSG)

    region_attractiveness = region_attractiveness_by_population(
        population_areas, buffer_distance, population_year, LOCAL_PROJ_CRS_EPSG
    )
    region_attractiveness.to_parquet(output_region_attractiveness_path)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    main()
