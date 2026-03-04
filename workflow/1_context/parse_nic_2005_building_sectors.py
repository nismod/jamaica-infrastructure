import itertools
import logging

import click
import geopandas
import pandas

from jamaica_infrastructure.utils import parse_jic2005


@click.command()
@click.version_option("1.0")
@click.option(
    "--buildings",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to buildings GeoPackage",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(exists=False, dir_okay=True, file_okay=True, writable=True),
    help="Path to buildings GeoParquet output",
)
def main(buildings, output):
    # List columns to keep - comment out those to drop which are in the base file (as of Mar-2026)
    keep_columns = [
        "osm_id",
        "osm_way_id",
        "type",
        "name",
        "addr_housenumber",
        "addr_street",
        "addr_city",
        "usage",
        "construction",
        "other_tags",
        "sector_code",
        "subsector_code",
        "assigned_attribute",
        "infra_type",
        "area_sqm",
        # "to_remove",
        # "remove",
        "building_type",
        "min_damage_cost",
        "max_damage_cost",
        "mean_damage_cost",
        "cost_unit",
        "residential_type",
        "residential_population",
        # "GDP_unit",
        "ED_ID",
        "ED",
        "PARISH",
        "CONST_NAME",
        # "A_GDP",
        # "B_GDP",
        # "C_GDP",
        # "D_GDP",
        # "E_GDP",
        # "F_GDP",
        # "G_GDP",
        # "H_GDP",
        # "I_GDP",
        # "J_GDP",
        # "K_GDP",
        # "L_GDP",
        # "M_GDP",
        # "N_GDP",
        # "O_GDP",
        # "total_GDP",
        "geometry",
    ]

    buildings_gdf = geopandas.read_file(
        buildings, columns=keep_columns, use_arrow=True, engine="pyogrio"
    )
    buildings_gdf[["jic2005_two_digit", "jic2005_three_digit"]] = (
        buildings_gdf.subsector_code.apply(subsector_to_jic2005)
    )
    buildings_gdf.to_parquet(output)


def set_encode(cs):
    return ",".join(sorted(list(set(cs))))


def subsector_to_jic2005(subsector_codes_string):
    codes = set(subsector_codes_string.split(","))
    two_or_three_digit = list(
        itertools.chain.from_iterable(parse_jic2005(c) for c in codes)
    )
    two_digit = (c[:2] for c in two_or_three_digit)
    three_digit = (c for c in two_or_three_digit if len(c) == 3)
    parsed = pandas.Series((set_encode(two_digit), set_encode(three_digit)))

    return parsed


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    main()
