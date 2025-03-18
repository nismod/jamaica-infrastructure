"""Estimate direct damages to physical assets exposed to hazards

"""

import os
from pathlib import Path

import click
import pandas as pd
import geopandas as gpd
import numpy as np

from analysis_utils import get_asset
from tqdm import tqdm

tqdm.pandas()

def quantiles(dataframe, grouping_by_columns, grouped_columns):
    grouped = (
        dataframe.groupby(grouping_by_columns, dropna=False)[grouped_columns]
        .agg([np.min, np.mean, np.max])
        .reset_index()
    )
    grouped.columns = grouping_by_columns + [
        f"{prefix}_{agg_name}"
        for prefix, agg_name in grouped.columns
        if prefix not in grouping_by_columns
    ]

    return grouped


@click.command()
@click.version_option("1.0")
@click.option(
    "--gpkg",
    "-g",
    required=True,
    help="asset_gpkg value in the asset definition CSV",
)
@click.option(
    "--layer",
    "-l",
    required=True,
    help="asset_layer value in the asset definition CSV",
)
@click.option(
    "--out-dir",
    "-o",
    required=True,
    help="Directory to save the results.",
)
@click.option(
    "--damage-dir",
    "-d",
    required=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=False, readable=True),
    help="Directory containing direct damage and EAD-EAEL damage results",
)
@click.option(
    "--single-failure-scenarios",
    "-s",
    required=False,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the single failure scenarios file",
)
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the asset definition file",
)
@click.option(
    "--parameter-file",
    "-p",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the parameter set definition file",
)
def loss_summary(
        gpkg,
        layer,
        out_dir,
        damage_dir,
        single_failure_scenarios,
        network_csv,
        parameter_file,
):
    """
    Collate direct damages and losses to an asset across all hazards under all parameter sets.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Read input files
    asset = get_asset(gpkg, layer, network_csv)
    parameter_sets = pd.read_csv(
        parameter_file,
        header=None,
        names=["parameter_set", "cost_uncertainty_parameter", "damage_uncertainty_parameter"]
    )
    direct_damages = [
        pd.read_parquet(
            os.path.join(
                damage_dir,
                f"{gpkg}_{layer}_direct_damages_parameter_set_{param.parameter_set}.parquet",
            )
        )for param in parameter_sets
    ]
    EAD_EAEL_damages = [
        pd.read_csv(
            os.path.join(
                damage_dir,
                f"{gpkg}_{layer}_EAD_EAEL_parameter_set_{param.parameter_set}.csv",
            )
        ) for param in parameter_sets
    ]
    if single_failure_scenarios:
        if asset.sector == "buildings":
            single_failure_df = gpd.read_file(single_failure_scenarios, layer="areas")
            single_failure_df = single_failure_df.rename(columns={"total_GDP": "economic_loss"}, inplace=True)
        else:
            single_failure_df = gpd.read_file(single_failure_scenarios)
            if gpkg == "potable_facilities_NWC":
                single_failure_df[asset.asset_id_column] = single_failure_df.progress_apply(
                    lambda x: str(x[asset.asset_id_column])
                    .lower()
                    .replace(" ", "_")
                    .replace(".0", ""),
                    axis=1,
                )
    else:
        single_failure_df = None

    # Calculate exposures
    exposures = direct_damages[0].copy()
    hazard_columns = [
        c
        for c in exposures.columns.values.tolist()
        if c
           not in [
               asset.asset.asset_id_column_column,
               "exposure_unit",
               "damage_cost_unit",
               "damage_uncertainty_parameter",
               "cost_uncertainty_parameter",
               "exposure",
           ]
    ]
    exposures[hazard_columns] = exposures["exposure"].to_numpy()[
                                :, None
                                ] * np.where(exposures[hazard_columns] > 0, 1, 0)

    sum_dict = dict([(hk, "sum") for hk in hazard_columns])
    exposures = (
        exposures.groupby([asset.asset.asset_id_column_column, "exposure_unit"], dropna=False)
        .agg(sum_dict)
        .reset_index()
    )
    exposures.to_parquet(
        os.path.join(
            out_dir,
            f"{gpkg}_{layer}_exposures.parquet",
        ),
        index=False,
    )
    
    # Collate direct damages
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
                single_failure_df[[asset.asset.asset_id_column_column, "economic_loss"]],
                how="left",
                on=[asset.asset.asset_id_column_column],
            ).fillna(0)
            df["economic_loss_unit"] = "J$/day"
            loss = df.copy()
            loss[hazard_columns] = loss["economic_loss"].to_numpy()[
                                   :, None
                                   ] * np.where(loss[hazard_columns] > 0, 1, 0)
            losses.append(
                loss[[asset.asset_id_column, "economic_loss_unit"] + hazard_columns]
            )

    damages = pd.concat(damages, axis=0, ignore_index=True)
    if len(damages.index) > 0:
        damages = quantiles(
            damages, [asset.asset_id_column, "damage_cost_unit"], hazard_columns
        )
        damages.to_parquet(
            os.path.join(out_dir, f"{gpkg}_{layer}_damages.parquet"),
            index=False,
        )
    else:
        Path(os.path.join(out_dir, f"{gpkg}_{layer}_damages.parquet")).touch()

    if len(losses) > 0:
        losses = pd.concat(losses, axis=0, ignore_index=True)
        if len(losses.index) > 0:
            losses = quantiles(
                losses, [asset.asset_id_column, "economic_loss_unit"], hazard_columns
            )
            losses.to_parquet(
                os.path.join(out_dir, f"{gpkg}_{layer}_losses.parquet"),
                index=False,
            )
        else:
            Path(os.path.join(out_dir, f"{gpkg}_{layer}_damages.parquet")).touch()
    else:
        Path(os.path.join(out_dir, f"{gpkg}_{layer}_damages.parquet")).touch()

    # Process the EAD and EAEL results
    for df in EAD_EAEL_damages:
        df["rcp"] = df["rcp"].astype(str)
        df["epoch"] = df["epoch"].astype(str)

    haz_rcp_epochs = list(
        set(
            EAD_EAEL_damages[0]
            .set_index(["hazard", "rcp", "epoch"])
            .index.values.tolist()
        )
    )
    summarised_damages = []
    for i, (haz, rcp, epoch) in enumerate(haz_rcp_epochs):
        damages = [
            df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch)]
            for df in EAD_EAEL_damages
        ]
        damages = pd.concat(damages, axis=0, ignore_index=True)
        damages.drop("confidence", axis=1, inplace=True)

        index_columns = [
            c
            for c in damages.columns.values.tolist()
            if ("EAD_" not in c) and ("EAEL_" not in c)
        ]
        index_columns = [
            i for i in index_columns if i not in ["cost_uncertainty_parameter", "damage_uncertainty_parameter"]
        ]
        damage_columns = [
            c
            for c in damages.columns.values.tolist()
            if ("EAD_" in c) or ("EAEL_" in c)
        ]

        if len(damages.index) > 0:
            summarised_damages.append(
                quantiles(damages, index_columns, damage_columns)
            )
    summarised_damages = pd.concat(
        summarised_damages, axis=0, ignore_index=True
    )
    summarised_damages.to_csv(
        os.path.join(out_dir, f"{gpkg}_{layer}_EAD_EAEL.csv"),
        index=False,
    )


if __name__ == "__main__":
    loss_summary()
