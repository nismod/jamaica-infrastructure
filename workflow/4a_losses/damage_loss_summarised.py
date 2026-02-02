"""
Find summary statistics of direct damages and losses over various realisations
of cost and uncertainty input data.
"""

import logging
import os
from pathlib import Path
from typing import Optional

import click
import geopandas as gpd
import pandas as pd
import numpy as np
from tqdm import tqdm

from jamaica_infrastructure.utils import get_asset, numeric_only_dataframe

tqdm.pandas()


def prepend_a_to_min_max(x: str) -> str:
    if x.endswith("_min"):
        return x[:-4] + "_amin"
    elif x.endswith("_max"):
        return x[:-4] + "_amax"
    else:
        return x


def quantiles(dataframe, grouping_by_columns, grouped_columns):
    assert numeric_only_dataframe(dataframe[grouped_columns])
    grouped = dataframe.groupby(grouping_by_columns, dropna=False)[grouped_columns].agg(["min", "mean", "max"]).reset_index()
    grouped.columns = grouping_by_columns + [
        f"{prefix}_{agg_name}" for prefix, agg_name in grouped.columns if prefix not in grouping_by_columns
    ]
    # downstream processes (including in irv-jamaica) assume 'amin' and 'amax' naming scheme
    grouped = grouped.rename(columns=prepend_a_to_min_max)
    return grouped


def read_single_failure_scenarios(path: Optional[str], asset: pd.Series) -> Optional[pd.DataFrame]:
    if not path:
        return None

    _, ext = os.path.splitext(path)
    if asset.sector == "buildings":
        if ext.lower() != ".gpkg":
            raise ValueError(f"Expect buildings single_failure_scenarios files to be GPKG format, received: {ext}")
        single_failure_df = gpd.read_file(path, layer="areas")
        single_failure_df = single_failure_df.rename(columns={"total_GDP": "economic_loss"})
    else:
        if ext.lower() != ".csv":
            raise ValueError(f"Expect most single_failure_scenarios files to be CSV format, received: {ext}")
        single_failure_df = pd.read_csv(path)
        if asset.asset_gpkg == "potable_facilities_NWC":
            single_failure_df[asset.asset_id_column] = single_failure_df.progress_apply(
                lambda x: str(x[asset.asset_id_column]).lower().replace(" ", "_").replace(".0", ""),
                axis=1,
            )
    return single_failure_df


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the asset definition file",
)
@click.option(
    "--damages",
    "-d",
    "damage_paths",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Input damage files (one per parameter set)",
)
@click.option(
    "--ead-eael",
    "-ee",
    "ead_eael_paths",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="EAD and EAEL files (one per parameter set)",
)
@click.option(
    "--single-failure-scenarios",
    "-s",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the single failure scenarios file",
)
@click.option(
    "--asset-gpkg",
    "-g",
    required=True,
    help="asset_gpkg value in the asset definition CSV",
)
@click.option(
    "--asset-layer",
    "-l",
    required=True,
    help="asset_layer value in the asset definition CSV",
)
@click.option(
    "--output-exposures",
    "-oe",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to write summarised exposures to",
)
@click.option(
    "--output-damages",
    "-od",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to write summarised damages to",
)
@click.option(
    "--output-losses",
    "-ol",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to write summarised losses to",
)
@click.option(
    "--output-ead-eael",
    "-oee",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to write summarised EAD and EAEL to",
)
def loss_summary(
    network_csv,
    damage_paths,
    ead_eael_paths,
    single_failure_scenarios,
    asset_gpkg,
    asset_layer,
    output_exposures,
    output_damages,
    output_losses,
    output_ead_eael,
):
    """
    Collate direct damages and losses to an asset across all hazards under all parameter sets.
    """

    logging.info(f"{asset_gpkg=} {asset_layer=}")
    asset: pd.Series = get_asset(network_csv, asset_gpkg, asset_layer)

    logging.info("Calculating exposures")
    exposures = pd.read_parquet(damage_paths[0])  # exposure the same for all damage files
    non_hazard_columns = [
        asset.asset_id_column,
        "exposure_unit",
        "damage_cost_unit",
        "damage_uncertainty_parameter",
        "cost_uncertainty_parameter",
        "exposure",
    ]
    hazard_columns = [c for c in exposures.columns.values.tolist() if c not in non_hazard_columns]
    exposures[hazard_columns] = exposures["exposure"].to_numpy()[:, None] * np.where(exposures[hazard_columns] > 0, 1, 0)

    sum_dict = dict([(hk, "sum") for hk in hazard_columns])
    exposures = exposures.groupby([asset.asset_id_column, "exposure_unit"], dropna=False).agg(sum_dict).reset_index()

    logging.info("Writing exposures to disk")
    exposures.to_parquet(output_exposures, index=False)
    del exposures

    logging.info("Reading single failure scenarios")
    single_failure_df: pd.DataFrame = read_single_failure_scenarios(single_failure_scenarios, asset)

    logging.info("Collating damages and losses")
    damages = []
    losses = []
    for damage_path in damage_paths:
        logging.info(damage_path)
        df = pd.read_parquet(damage_path)
        df = (
            df.groupby(
                [
                    asset.asset_id_column,
                    "damage_cost_unit",
                ],
                dropna=False,
            )
            .agg(sum_dict)
            .reset_index()
        )
        damages.append(df)
        if single_failure_df is not None:
            df = pd.merge(
                df,
                single_failure_df[[asset.asset_id_column, "economic_loss"]],
                how="left",
                on=[asset.asset_id_column],
            ).fillna(0)
            df["economic_loss_unit"] = "J$/day"
            loss = df.copy()
            loss[hazard_columns] = loss["economic_loss"].to_numpy()[:, None] * np.where(loss[hazard_columns] > 0, 1, 0)
            losses.append(loss[[asset.asset_id_column, "economic_loss_unit"] + hazard_columns])

    logging.info("Writing damages to disk")
    damages = pd.concat(damages, axis=0, ignore_index=True)
    if len(damages.index) > 0:
        damages = quantiles(damages, [asset.asset_id_column, "damage_cost_unit"], hazard_columns)
        damages.to_parquet(output_damages, index=False)
    else:
        Path(output_damages).touch()
    del damages

    logging.info("Writing losses to disk")
    if len(losses) > 0:
        losses = pd.concat(losses, axis=0, ignore_index=True)
        if len(losses.index) > 0:
            losses = quantiles(losses, [asset.asset_id_column, "economic_loss_unit"], hazard_columns)
            losses.to_parquet(output_losses, index=False)
        else:
            Path(output_losses).touch()
    else:
        Path(output_losses).touch()
    del losses

    logging.info("Reading EAD and EAEL")
    EAD_EAEL_damages = [pd.read_csv(file, dtype={"rcp": str, "epoch": str}) for file in ead_eael_paths]

    logging.info("Processing EAD and EAEL")
    haz_rcp_epochs = list(set(EAD_EAEL_damages[0].set_index(["hazard", "rcp", "epoch"]).index.values.tolist()))
    summarised_damages = []
    for haz, rcp, epoch in haz_rcp_epochs:
        logging.info(f"{haz} {rcp} {epoch}")
        damages = [df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch)] for df in EAD_EAEL_damages]
        damages = pd.concat(damages, axis=0, ignore_index=True)
        damages.drop("confidence", axis=1, inplace=True)

        index_columns = [c for c in damages.columns.values.tolist() if ("EAD_" not in c) and ("EAEL_" not in c)]
        index_columns = [i for i in index_columns if i not in ["cost_uncertainty_parameter", "damage_uncertainty_parameter"]]
        damage_columns = [c for c in damages.columns.values.tolist() if ("EAD_" in c) or ("EAEL_" in c)]

        if len(damages.index) > 0:
            summarised_damages.append(quantiles(damages, index_columns, damage_columns))
    summarised_damages = pd.concat(summarised_damages, axis=0, ignore_index=True)

    logging.info("Writing EAD and EAEL to disk")
    summarised_damages.to_csv(output_ead_eael, index=False)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    loss_summary()
