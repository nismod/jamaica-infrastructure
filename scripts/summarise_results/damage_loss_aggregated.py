"""Estimate direct damages to physical assets exposed to hazards

"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from summarise_utils import *
from tqdm import tqdm

tqdm.pandas()


def main(config, damages_losses_folder, summary_results_folder, network_csv):

    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]

    damages_losses_results = os.path.join(output_data_path, damages_losses_folder)

    summary_results = os.path.join(output_data_path, summary_results_folder)
    if os.path.exists(summary_results) == False:
        os.mkdir(summary_results)

    damages_totals = []
    loss_totals = []
    asset_data_details = pd.read_csv(network_csv)
    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column

        if (
            os.path.isfile(
                os.path.join(
                    damages_losses_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_damages.parquet",
                )
            )
            is True
        ):
            damages = pd.read_parquet(
                os.path.join(
                    damages_losses_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_damages.parquet",
                )
            )
            haz_columns = [
                c
                for c in damages.columns.values.tolist()
                if c not in [asset_id, "damage_cost_unit"]
            ]
            df = pd.DataFrame()
            df["sector"] = [asset_info.sector]
            df["asset_gpkg"] = [asset_info.asset_gpkg]
            df["asset_layer"] = [asset_info.asset_layer]
            df["damage_cost_unit"] = [damages["damage_cost_unit"].values[0]]
            df[haz_columns] = damages[haz_columns].sum(axis=0)
            damages_totals.append(df)
            del df

        if (
            os.path.isfile(
                os.path.join(
                    damages_losses_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_losses.parquet",
                )
            )
            is True
        ):
            losses = pd.read_parquet(
                os.path.join(
                    damages_losses_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_losses.parquet",
                )
            )
            haz_columns = [
                c
                for c in losses.columns.values.tolist()
                if c not in [asset_id, "economic_loss_unit"]
            ]
            df = pd.DataFrame()
            df["sector"] = [asset_info.sector]
            df["asset_gpkg"] = [asset_info.asset_gpkg]
            df["asset_layer"] = [asset_info.asset_layer]
            df["economic_loss_unit"] = [losses["economic_loss_unit"].values[0]]
            df[haz_columns] = losses[haz_columns].sum(axis=0)
            loss_totals.append(df)
            del df

    if len(damages_totals) > 0:
        damages = pd.concat(damages_totals, axis=0, ignore_index=True).fillna(0)
        damages.to_csv(
            os.path.join(summary_results, f"sector_total_damages.csv"), index=False
        )
        del damages

    if len(loss_totals) > 0:
        losses = pd.concat(loss_totals, axis=0, ignore_index=True).fillna(0)
        losses.to_csv(
            os.path.join(summary_results, f"sector_level_losses.csv"), index=False
        )
        del losses


if __name__ == "__main__":
    CONFIG = load_config()
    try:
        direct_damages_folder = str(sys.argv[1])
        summary_results_folder = str(sys.argv[2])
        network_csv = str(sys.argv[3])
    except IndexError:
        print("Got arguments", sys.argv)
        exit()

    main(CONFIG, direct_damages_folder, summary_results_folder, network_csv)
