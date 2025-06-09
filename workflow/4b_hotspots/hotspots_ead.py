"""Estimate expected annual (direct) damages from return period damages, output on hotspots grid."""

import logging

import click
import numpy as np
import pandas as pd
from scipy.integrate import simpson
from tqdm import tqdm

from jamaica_infrastructure.raster import (
    write_indexed_splits_to_tiff,
    write_constant_tiff,
)

tqdm.pandas()


def list_hazard_columns(df: pd.DataFrame) -> list[str]:
    """Pick columns named like hazard keys (containing return-period damage estimates)"""
    # hard-codes assumption in column naming
    return [colname for colname in df.columns if "__rp_" in colname]


def get_hazard_meta(df: pd.DataFrame) -> pd.DataFrame:
    """Construct dataframe with metadata about hazard for return-period damage estimates"""
    hazard_records = []
    for colname in list_hazard_columns(df):
        split_column_name = colname.split("__")
        hazard_type = split_column_name[0]
        return_period = int(split_column_name[1].split("_")[1])
        hazard_records.append(
            {"hazard": hazard_type, "rp": return_period, "key": colname}
        )

    hazard_meta = pd.DataFrame(hazard_records)
    hazard_meta["probability"] = 1 / hazard_meta.rp
    return hazard_meta


def ead(rp_damages: pd.DataFrame, meta: pd.DataFrame, hazard: str) -> np.ndarray:
    meta = meta.query(f'hazard == "{hazard}"').sort_values(by="probability")

    rp_cols = meta.key
    rps = meta.rp
    probabilities = 1 / rps

    return simpson(rp_damages[rp_cols], x=probabilities, axis=1)


@click.command()
@click.version_option("1.0")
@click.option(
    "--splits-path",
    "-s",
    required=True,
    help="Input return-period damages splits file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--grid-path",
    "-gr",
    required=True,
    help="Hotspots raster grid",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--hazard", "-h", required=True, help="hazard value in the hotspots_layers CSV"
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to save EAD hotspots",
)
def ead_to_grid(splits_path, grid_path, hazard, output_path):
    splits = pd.read_parquet(splits_path)
    ead_colname = f"{hazard}__ead"
    meta = get_hazard_meta(splits)

    if hazard not in set(meta.hazard):
        write_constant_tiff(np.nan, grid_path, output_path)
        return

    splits[ead_colname] = ead(splits, meta, hazard)
    write_indexed_splits_to_tiff(
        splits,
        grid_path=grid_path,
        output_path=output_path,
        grid_name=hazard,
        value_colname=ead_colname,
        index_x="cell_index_0_x",
        index_y="cell_index_0_y",
        no_data_value=np.nan,
    )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO
    )
    ead_to_grid()
