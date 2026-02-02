"""
Calculate Expected Annual Damages from direct damages at a set of return periods
Calculate Expected Annual Economic Losses from single point of failure results and return period maps
"""

import logging
import os
import warnings

import click
import geopandas as gpd
import numpy as np
import pandas as pd
from tqdm import tqdm

from jamaica_infrastructure.analysis.utils import risks

warnings.simplefilter(action="ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
tqdm.pandas()


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the asset layer table",
)
@click.option(
    "--hazard-csv",
    "-h",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the hazard table",
)
@click.option(
    "--asset-gpkg-label",
    "-g",
    required=True,
    help="asset_gpkg value in the network CSV",
)
@click.option(
    "--asset-layer",
    "-l",
    required=True,
    help="asset_layer value in the network CSV",
)
@click.option(
    "--sensitivity-csv",
    "-s",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to sensitivity parameters table",
)
@click.option(
    "--sensitivity-id",
    "-id",
    required=True,
    type=int,
    help="Sensitivity parameter set ID",
)
@click.option(
    "--damage-file",
    "-d",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to previously calculated direct damages",
)
@click.option(
    "--single-failure-scenarios",
    "-sfs",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to single failure scenarios",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to save direct damages to",
)
def EAD_EAEL(
    network_csv,
    hazard_csv,
    asset_gpkg_label,
    asset_layer,
    sensitivity_csv,
    sensitivity_id,
    damage_file,
    single_failure_scenarios,
    output_path,
):

    bridge_flood_protection = (
        50  # Bridges in Jamaica are designed to withstand 1 in 50 year floods
    )

    logging.info("Read sensitivity data")
    sensitivity = pd.read_csv(sensitivity_csv)
    cost_uncertainty_parameter = sensitivity.loc[sensitivity_id, "cost_uncertainty_parameter"].squeeze()
    damage_uncertainty_parameter = sensitivity.loc[sensitivity_id, "damage_uncertainty_parameter"].squeeze()

    logging.info("Read asset data")
    asset_data_details = pd.read_csv(network_csv)
    asset_info = asset_data_details.loc[
        (asset_data_details.asset_gpkg == asset_gpkg_label) \
        & (asset_data_details.asset_layer == asset_layer)
    ].squeeze()
    asset_id = asset_info.asset_id_column

    if not os.path.isfile(damage_file):
        return  # short circuit if the file does not exist -- is this necessary?

    logging.info("Read damage data")
    expected_damages = []
    df = pd.read_parquet(damage_file)
    hazard_columns = [
        c
        for c in df.columns.values.tolist()
        if c not in [asset_id, "damage_cost_unit", "exposure"]
    ]

    df = (
        df.groupby([asset_id, "damage_cost_unit"])[hazard_columns]
        .sum()
        .reset_index()
    )
    loss_column = []
    if asset_info.single_failure_scenarios != "none":
        logging.info("Read single failure scenarios")
        loss_column = ["economic_loss"]
        if asset_info.sector != "buildings":
            loss_df = pd.read_csv(single_failure_scenarios)
            if asset_info.asset_gpkg == "potable_facilities_NWC":
                loss_df[asset_id] = loss_df.progress_apply(
                    lambda x: str(x[asset_id])
                    .lower()
                    .replace(" ", "_")
                    .replace(".0", ""),
                    axis=1,
                )
        else:
            loss_df = gpd.read_file(single_failure_scenarios, layer="areas")
            loss_df.rename(columns={"total_GDP": "economic_loss"}, inplace=True)
        df = pd.merge(
            df,
            loss_df[[asset_info.asset_id_column, "economic_loss"]],
            how="left",
            on=[asset_info.asset_id_column],
        ).fillna(0)
        df["economic_loss_unit"] = "J$/day"
    else:
        df["economic_loss_unit"] = "None"

    hazard_data_details = pd.read_csv(hazard_csv, encoding="latin1").fillna(0)
    hazard_data_details = hazard_data_details[
        hazard_data_details.key.isin(hazard_columns)
    ]
    haz_rcp_epoch_confidence = list(
        set(
            hazard_data_details.set_index(
                ["hazard", "rcp", "epoch", "confidence"]
            ).index.values.tolist()
        )
    )
    for i, (haz, rcp, epoch, confidence) in enumerate(haz_rcp_epoch_confidence):
        logging.info(f"Processing {haz} {rcp} {epoch} {confidence}")
        index_columns = [asset_id, "damage_cost_unit", "economic_loss_unit"]
        haz_df = hazard_data_details[
            (hazard_data_details.hazard == haz)
            & (hazard_data_details.rcp == rcp)
            & (hazard_data_details.epoch == epoch)
            & (hazard_data_details.confidence == confidence)
        ]
        haz_cols, haz_rps = map(
            list,
            list(
                zip(
                    *sorted(
                        list(
                            zip(
                                haz_df.key.values.tolist(),
                                haz_df.rp.values.tolist(),
                            )
                        ),
                        key=lambda x: x[-1],
                        reverse=True,
                    )
                )
            ),
        )

        haz_prob = [1.0 / rp for rp in haz_rps]
        damages = df[index_columns + loss_column + haz_cols]
        damages["hazard"] = haz
        damages["rcp"] = rcp
        damages["epoch"] = epoch
        damages["confidence"] = confidence
        damages.columns = (
            index_columns
            + loss_column
            + haz_prob
            + ["hazard", "rcp", "epoch", "confidence"]
        )
        index_columns += ["hazard", "rcp", "epoch", "confidence"]
        damages = damages[damages[haz_prob].sum(axis=1) > 0]
        losses = damages.copy()
        losses.columns = losses.columns.map(str)

        if (
            (asset_info.asset_gpkg == "roads")
            and (asset_info.asset_layer == "nodes")
            and (haz in ["coastal", "fluvial", "surface"])
        ):
            damages["protection_standard"] = bridge_flood_protection
            expected_damage_df = risks(
                damages,
                index_columns + ["protection_standard"],
                haz_prob,
                "EAD",
                flood_protection_period=bridge_flood_protection,
                flood_protection_name=flood_protection_name,
            )
            expected_damage_df["protection_standard"] = bridge_flood_protection
            if "economic_loss" in damages.columns.values.tolist():
                losses[[str(h) for h in haz_prob]] = losses[
                    "economic_loss"
                ].to_numpy()[:, None] * np.where(
                    losses[[str(h) for h in haz_prob]] > 0, 1, 0
                )
                losses["protection_standard"] = bridge_flood_protection
                protected_loss_df = risks(
                    losses,
                    index_columns + ["protection_standard"],
                    haz_prob,
                    "EAEL",
                    flood_protection_period=bridge_flood_protection,
                    flood_protection_name=flood_protection_name,
                )
                expected_damage_df = pd.merge(
                    expected_damage_df,
                    protected_loss_df,
                    how="left",
                    on=index_columns + ["protection_standard"],
                ).fillna(0)
                del protected_loss_df
        else:
            expected_damage_df = risks(
                damages,
                index_columns,
                haz_prob,
                "EAD",
                flood_protection_name=flood_protection_name,
            )
            if "economic_loss" in damages.columns.values.tolist():
                losses[[str(h) for h in haz_prob]] = losses[
                    "economic_loss"
                ].to_numpy()[:, None] * np.where(
                    losses[[str(h) for h in haz_prob]] > 0, 1, 0
                )
                economic_loss_df = risks(
                    losses,
                    index_columns,
                    haz_prob,
                    "EAEL",
                    flood_protection_name=flood_protection_name,
                )
                expected_damage_df = pd.merge(
                    expected_damage_df,
                    economic_loss_df,
                    how="left",
                    on=index_columns,
                ).fillna(0)
                del economic_loss_df

        expected_damages.append(expected_damage_df)
        del expected_damage_df

    expected_damages = pd.concat(expected_damages, axis=0, ignore_index=True)
    expected_loss_columns = [
        c
        for c in expected_damages.columns.values.tolist()
        if "EAD_" in c or "EAEL_" in c
    ]
    expected_damages = expected_damages[
        expected_damages[expected_loss_columns].sum(axis=1) > 0
    ]
    expected_damages["cost_uncertainty_parameter"] = cost_uncertainty_parameter
    expected_damages["damage_uncertainty_parameter"] = (
        damage_uncertainty_parameter
    )

    logging.info("Writing expected damages and losses to disk")
    expected_damages.to_csv(output_path, index=False)


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    EAD_EAEL()
