"""Estimate direct damages to physical assets exposed to hazards

"""

import logging
import os
from pathlib import Path

import click
import geopandas as gpd
import pandas as pd
import numpy as np
from tqdm import tqdm

from jamaica_infrastructure.utils import get_asset, numeric_only_dataframe

tqdm.pandas()


def quantiles(dataframe, grouping_by_columns, grouped_columns):
    assert numeric_only_dataframe(dataframe[grouped_columns])
    grouped = dataframe.groupby(grouping_by_columns, dropna=False)[grouped_columns].agg(["min", "mean", "max"]).reset_index()
    grouped.columns = grouping_by_columns + [
        f"{prefix}_{agg_name}" for prefix, agg_name in grouped.columns if prefix not in grouping_by_columns
    ]
    return grouped


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
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Input damage files (one per parameter set)",
)
@click.option(
    "--ead-eael",
    "-ee",
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
    damages,
    ead_eael,
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
    asset = get_asset(network_csv, asset_gpkg, asset_layer)

    logging.info("Reading exposure and direct damages")
    direct_damages = [pd.read_parquet(file) for file in damages]

    logging.info("Reading EAD and EAEL")
    EAD_EAEL_damages = [pd.read_csv(file) for file in ead_eael]

    logging.info("Reading single failure scenarios")
    if single_failure_scenarios:
        _, ext = os.path.splitext(single_failure_scenarios)
        if asset.sector == "buildings":
            if ext.lower() != ".gpkg":
                raise ValueError(f"Expect buildings single_failure_scenarios files to be GPKG format, received: {ext}")
            single_failure_df = gpd.read_file(single_failure_scenarios, layer="areas")
            single_failure_df = single_failure_df.rename(columns={"total_GDP": "economic_loss"}, inplace=True)
        else:
            if ext.lower() != ".csv":
                raise ValueError(f"Expect most single_failure_scenarios files to be CSV format, received: {ext}")
            single_failure_df = pd.read_csv(single_failure_scenarios)
            if asset_gpkg == "potable_facilities_NWC":
                single_failure_df[asset.asset_id_column] = single_failure_df.progress_apply(
                    lambda x: str(x[asset.asset_id_column]).lower().replace(" ", "_").replace(".0", ""),
                    axis=1,
                )
    else:
        single_failure_df = None

    logging.info("Calculating exposures")
    exposures = direct_damages[0].copy()
    hazard_columns = [
        c
        for c in exposures.columns.values.tolist()
        if c
        not in [
            asset.asset_id_column,
            "exposure_unit",
            "damage_cost_unit",
            "damage_uncertainty_parameter",
            "cost_uncertainty_parameter",
            "exposure",
        ]
    ]
    exposures[hazard_columns] = exposures["exposure"].to_numpy()[:, None] * np.where(exposures[hazard_columns] > 0, 1, 0)

    sum_dict = dict([(hk, "sum") for hk in hazard_columns])
    exposures = exposures.groupby([asset.asset_id_column, "exposure_unit"], dropna=False).agg(sum_dict).reset_index()
    exposures.to_parquet(output_exposures, index=False)

    logging.info("Collating damages and losses")
    damages = []
    losses = []
    for df in direct_damages:
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

    logging.info("Writing outputs to disk")
    damages = pd.concat(damages, axis=0, ignore_index=True)
    if len(damages.index) > 0:
        damages = quantiles(damages, [asset.asset_id_column, "damage_cost_unit"], hazard_columns)
        damages.to_parquet(output_damages, index=False)
    else:
        Path(output_damages).touch()

    if len(losses) > 0:
        losses = pd.concat(losses, axis=0, ignore_index=True)
        if len(losses.index) > 0:
            losses = quantiles(losses, [asset.asset_id_column, "economic_loss_unit"], hazard_columns)
            losses.to_parquet(output_losses, index=False)
        else:
            Path(output_losses).touch()
    else:
        Path(output_losses).touch()

    # Process the EAD and EAEL results
    for df in EAD_EAEL_damages:
        df["rcp"] = df["rcp"].astype(str)
        df["epoch"] = df["epoch"].astype(str)

    haz_rcp_epochs = list(set(EAD_EAEL_damages[0].set_index(["hazard", "rcp", "epoch"]).index.values.tolist()))
    summarised_damages = []
    for i, (haz, rcp, epoch) in enumerate(haz_rcp_epochs):
        damages = [df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch)] for df in EAD_EAEL_damages]
        damages = pd.concat(damages, axis=0, ignore_index=True)
        damages.drop("confidence", axis=1, inplace=True)

        index_columns = [c for c in damages.columns.values.tolist() if ("EAD_" not in c) and ("EAEL_" not in c)]
        index_columns = [i for i in index_columns if i not in ["cost_uncertainty_parameter", "damage_uncertainty_parameter"]]
        damage_columns = [c for c in damages.columns.values.tolist() if ("EAD_" in c) or ("EAEL_" in c)]

        if len(damages.index) > 0:
            summarised_damages.append(quantiles(damages, index_columns, damage_columns))
    summarised_damages = pd.concat(summarised_damages, axis=0, ignore_index=True)
    summarised_damages.to_csv(output_ead_eael, index=False)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    loss_summary()
