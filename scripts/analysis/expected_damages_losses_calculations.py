"""Estimate direct damages to physical assets exposed to hazards

"""

import sys
import os

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

import geopandas as gpd
import numpy as np

from analysis_utils import *
from tqdm import tqdm

tqdm.pandas()


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    direct_damages_results = os.path.join(output_data_path, "direct_damages")

    asset_data_details = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks",
            "network_layers_hazard_intersections_details.csv",
        )
    )
    bridge_flood_protection = (
        50  # Bridges in Jamaica are designed to withstand 1 in 50 year floods
    )

    hazard_data_file = os.path.join(
        processed_data_path, "hazards", "hazard_descriptions", "hazard_layers.csv"
    )
    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column
        asset_damages_results = os.path.join(
            direct_damages_results, f"{asset_info.asset_gpkg}_{asset_info.asset_layer}"
        )
        param_values = pd.read_csv(f"{processed_data_path}/sensitivity_parameters.csv")
        for param in param_values.itertuples():
            parameter_set = param.set_id
            damage_file = os.path.join(
                asset_damages_results,
                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_direct_damages_parameter_set_{parameter_set}.parquet",
            )
            if os.path.isfile(damage_file) is True:
                expected_damages = []
                df = pd.read_parquet(damage_file)
                hazard_columns = [
                    c
                    for c in df.columns.values.tolist()
                    if c not in [asset_id, "damage_cost_unit" "exposure"]
                ]
                df = (
                    df.groupby([asset_id, "damage_cost_unit"])[hazard_columns]
                    .sum()
                    .reset_index()
                )
                loss_column = []
                if asset_info.single_failure_scenarios != "none":
                    loss_column = ["economic_loss"]
                    if asset_info.sector != "buildings":
                        loss_df = pd.read_csv(
                            os.path.join(
                                output_data_path, asset_info.single_failure_scenarios
                            )
                        )
                        if asset_info.asset_gpkg == "potable_facilities_NWC":
                            loss_df[asset_id] = loss_df.progress_apply(
                                lambda x: str(x[asset_id])
                                .lower()
                                .replace(" ", "_")
                                .replace(".0", ""),
                                axis=1,
                            )
                    else:
                        loss_df = gpd.read_file(
                            os.path.join(
                                processed_data_path, asset_info.single_failure_scenarios
                            ),
                            layer="areas",
                        )
                        loss_df.rename(
                            columns={"total_GDP": "economic_loss"}, inplace=True
                        )
                        # loss_df["loss_unit"] = "JD/day"
                    df = pd.merge(
                        df,
                        loss_df[[asset_info.asset_id_column, "economic_loss"]],
                        how="left",
                        on=[asset_info.asset_id_column],
                    ).fillna(0)
                    df["economic_loss_unit"] = "J$/day"
                else:
                    df["economic_loss_unit"] = "None"
                # haz_rcp_epoch_confidence = list(set(df.set_index(["hazard","rcp","epoch","confidence"]).index.values.tolist()))
                hazard_data_details = pd.read_csv(
                    os.path.join(hazard_data_file), encoding="latin1"
                ).fillna(0)
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
                # print (haz_rcp_epoch_confidence)
                for i, (haz, rcp, epoch, confidence) in enumerate(
                    haz_rcp_epoch_confidence
                ):
                    index_columns = [asset_id, "damage_cost_unit", "economic_loss_unit"]
                    haz_df = hazard_data_details[
                        (hazard_data_details.hazard == haz)
                        & (hazard_data_details.rcp == rcp)
                        & (hazard_data_details.epoch == epoch)
                        & (hazard_data_details.confidence == confidence)
                    ]
                    # haz_cols = haz_df.key.values.tolist()
                    # haz_rps = haz_df.rp.values.tolist()
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
                    # expected_damage_df = risks(damages,index_columns,haz_prob,
                    #                             None,'EAD',
                    #                             flood_protection=None)
                    expected_damage_df = risks(damages, index_columns, haz_prob, "EAD")
                    # print (expected_damage_df)
                    if "economic_loss" in damages.columns.values.tolist():
                        losses = damages.copy()
                        # for hz in haz_prob:
                        #     losses[str(hz)] = losses["economic_loss"]*np.where(losses[str(hz)]>0,1,0)
                        losses[haz_prob] = losses["economic_loss"].to_numpy()[
                            :, None
                        ] * np.where(losses[haz_prob] > 0, 1, 0)
                        # economic_loss_df = risks(losses,index_columns,haz_prob,
                        #                         None,'EAEL',
                        #                         flood_protection=None)
                        economic_loss_df = risks(
                            losses, index_columns, haz_prob, "EAEL"
                        )
                        expected_damage_df = pd.merge(
                            expected_damage_df,
                            economic_loss_df,
                            how="left",
                            on=index_columns,
                        ).fillna(0)
                        del economic_loss_df
                    if (
                        (asset_info.asset_gpkg == "roads")
                        and (asset_info.asset_layer == "nodes")
                        and (haz in ["coastal", "fluvial", "surface"])
                    ):
                        damages["protection_standard"] = bridge_flood_protection
                        expected_damage_df["protection_standard"] = (
                            bridge_flood_protection
                        )
                        # protected_damage_df = risks(damages,index_columns + ["protection_standard"],haz_prob,
                        #                         "protection_standard",'EAD',
                        #                         flood_protection="yes",flood_protection_name="designed_protection")
                        protected_damage_df = risks(
                            damages,
                            index_columns + ["protection_standard"],
                            haz_prob,
                            "EAD",
                            flood_protection_period=bridge_flood_protection,
                            flood_protection_name="designed_protection",
                        )
                        expected_damage_df = pd.merge(
                            expected_damage_df,
                            protected_damage_df,
                            how="left",
                            on=index_columns + ["protection_standard"],
                        ).fillna(0)
                        del protected_damage_df
                        if "economic_loss" in damages.columns.values.tolist():
                            losses["protection_standard"] = bridge_flood_protection
                            # protected_loss_df = risks(losses,index_columns + ["protection_standard"],haz_prob,
                            #                         "protection_standard",'EAEL',
                            #                         flood_protection="yes",flood_protection_name="designed_protection")
                            protected_loss_df = risks(
                                losses,
                                index_columns + ["protection_standard"],
                                haz_prob,
                                "EAEL",
                                flood_protection_period=bridge_flood_protection,
                                flood_protection_name="designed_protection",
                            )
                            expected_damage_df["protection_standard"] = (
                                bridge_flood_protection
                            )
                            expected_damage_df = pd.merge(
                                expected_damage_df,
                                protected_loss_df,
                                how="left",
                                on=index_columns + ["protection_standard"],
                            ).fillna(0)
                            del protected_loss_df
                    expected_damages.append(expected_damage_df)
                    del expected_damage_df

                expected_damages = pd.concat(
                    expected_damages, axis=0, ignore_index=True
                )
                expected_loss_columns = [
                    c
                    for c in expected_damages.columns.values.tolist()
                    if "EAD_" in c or "EAEL_" in c
                ]
                expected_damages = expected_damages[
                    expected_damages[expected_loss_columns].sum(axis=1) > 0
                ]
                asset_damages_results = os.path.join(
                    direct_damages_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}",
                )
                if os.path.exists(asset_damages_results) == False:
                    os.mkdir(asset_damages_results)
                expected_damages.to_csv(
                    os.path.join(
                        asset_damages_results,
                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_parameter_set_{parameter_set}_0.csv",
                    ),
                    index=False,
                )
        print(f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer}")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
