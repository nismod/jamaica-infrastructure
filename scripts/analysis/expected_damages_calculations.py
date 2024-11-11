"""Estimate direct damages to physical assets exposed to hazards

"""

import sys
import os

import pandas as pd
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

    param_values = pd.read_csv("parameter_combinations.txt", sep=" ")
    for asset_info in asset_data_details.itertuples():
        asset_damages_results = os.path.join(
            direct_damages_results, f"{asset_info.asset_gpkg}_{asset_info.asset_layer}"
        )
        for param in param_values.itertuples():
            damage_file = os.path.join(
                asset_damages_results,
                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_direct_damages_parameter_set_{param.parameter_set}.csv",
            )
            if os.path.isfile(damage_file) is True:
                expected_damages = []
                df = pd.read_csv(damage_file).fillna(0)
                haz_rcp_epoch_confidence = list(
                    set(
                        df.set_index(
                            ["hazard", "rcp", "epoch", "confidence"]
                        ).index.values.tolist()
                    )
                )
                # print (haz_rcp_epoch_confidence)
                for i, (haz, rcp, epoch, confidence) in enumerate(
                    haz_rcp_epoch_confidence
                ):
                    damages = df[
                        (df.hazard == haz)
                        & (df.rcp == rcp)
                        & (df.epoch == epoch)
                        & (df.confidence == confidence)
                    ]
                    # print (damages)
                    damages["probability"] = 1.0 / damages["rp"]
                    index_columns = [
                        c
                        for c in damages.columns.values.tolist()
                        if c
                        not in [
                            "rp",
                            "probability",
                            "direct_damage_cost",
                            "direct_reopen_cost",
                            "exposure",
                        ]
                    ]
                    expected_damage_df = risks_pivot(
                        damages,
                        index_columns,
                        "probability",
                        "direct_damage_cost",
                        None,
                        "EAD",
                        flood_protection=None,
                    )
                    # print (expected_damage_df)
                    if "direct_reopen_cost" in damages.columns.values.tolist():
                        expected_reopen_df = risks_pivot(
                            damages,
                            index_columns,
                            "probability",
                            "direct_reopen_cost",
                            None,
                            "EAR",
                            flood_protection=None,
                        )
                        expected_damage_df = pd.merge(
                            expected_damage_df,
                            expected_reopen_df,
                            how="left",
                            on=index_columns,
                        ).fillna(0)
                        del expected_reopen_df
                    expected_damages.append(expected_damage_df)
                    del expected_damage_df

                expected_damages = pd.concat(
                    expected_damages, axis=0, ignore_index=True
                )
                asset_damages_results = os.path.join(
                    direct_damages_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}",
                )
                if os.path.exists(asset_damages_results) == False:
                    os.mkdir(asset_damages_results)
                expected_damages.to_csv(
                    os.path.join(
                        asset_damages_results,
                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_parameter_set_{param.parameter_set}.csv",
                    ),
                    index=False,
                )
        print(f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer}")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
