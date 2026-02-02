"""
Find summary statistics of direct damages and losses over various realisations
of cost and uncertainty input data.
"""

import logging
import os

import click
import pandas as pd
import numpy as np
from tqdm import tqdm


def quantiles(dataframe, grouping_by_columns, grouped_columns):
    # quantiles_list = ['mean','min','max','median','q5','q95']
    grouped = dataframe.groupby(grouping_by_columns, dropna=False)[grouped_columns].agg(["min", "mean", "max"]).reset_index()
    grouped.columns = grouping_by_columns + [f"{prefix}_{agg_name}" for prefix, agg_name in grouped.columns if prefix not in grouping_by_columns]

    return grouped


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
    "--sensitivity-csv",
    "-s",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to sensitivity parameters table",
)
@click.option(
    "--asset-gpkg",
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
    "--damages",
    "-d",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Computed direct damages by parameter set",
)
@click.option(
    "--output-damages-path",
    "-od",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to save direct damages to",
)
@click.option(
    "--output-exposure-path",
    "-oe",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to save exposure to",
)
def collapse_sensitivity(
    network_csv,
    sensitivity_csv,
    asset_gpkg,
    asset_layer,
    damages,
    output_damages_path,
    output_exposure_path,
):

    logging.info("Reading metadata")
    asset_data_details = pd.read_csv(network_csv)
    param_values = pd.read_csv(sensitivity_csv)
    uncertainty_columns = ("cost_uncertainty_parameter", "damage_uncertainty_parameter")

    asset_info = asset_data_details.loc[(asset_data_details.asset_gpkg == asset_gpkg) & (asset_data_details.asset_layer == asset_layer)].squeeze()
    asset_id = asset_info.asset_id_column

    logging.info("Reading damage data")
    damage_results = [pd.read_parquet(file) for file in damages if os.path.isfile(file) is True]

    if damage_results:
        logging.info("Processing exposure")
        exposures = damage_results[0].copy()
        hazard_columns = [
            c
            for c in exposures.columns.values.tolist()
            if c
            not in [
                asset_id,
                "exposure_unit",
                "damage_cost_unit",
                "damage_uncertainty_parameter",
                "cost_uncertainty_parameter",
                "exposure",
            ]
        ]
        # for hz in hazard_columns:
        #     exposures[hz] = exposures["exposure"]*np.where(exposures[hz] > 0,1,0)
        exposures[hazard_columns] = exposures["exposure"].to_numpy()[:, None] * np.where(exposures[hazard_columns] > 0, 1, 0)

        sum_dict = dict([(hk, "sum") for hk in hazard_columns])
        exposures = exposures.groupby([asset_id, "exposure_unit"], dropna=False).agg(sum_dict).reset_index()
        exposures.to_csv(output_exposure_path, index=False)
        del exposures

        logging.info("Processing damages")
        damages = [
            df.groupby(
                [
                    asset_id,
                    "damage_cost_unit",
                ],
                dropna=False,
            )
            .agg(sum_dict)
            .reset_index()
            for df in damage_results
        ]

        damages = pd.concat(damages, axis=0, ignore_index=True)
        if len(damages.index) > 0:
            damages = quantiles(damages, [asset_id, "damage_cost_unit"], hazard_columns)
            # damages.to_parquet(os.path.join(summary_results,
            #             f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_damages.parquet"),index=False)
            damages.to_csv(output_damages_path, index=False)
        del damages

    logging.info("Processing EAD & EAEL")
    damage_files = [
        os.path.join(
            asset_damages_results,
            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_parameter_set_{param.parameter_set}.csv",
        )
        for param in param_values.itertuples()
    ]
    damage_results = [pd.read_csv(file) for file in damage_files if os.path.isfile(file) is True]

    if damage_results:
        print([len(df.index) for df in damage_results])
        for df in damage_results:
            df["rcp"] = df["rcp"].astype(str)
            df["epoch"] = df["epoch"].astype(str)
        haz_rcp_epochs = list(set(damage_results[0].set_index(["hazard", "rcp", "epoch"]).index.values.tolist()))
        print(haz_rcp_epochs)
        summarised_damages = []
        for i, (haz, rcp, epoch) in enumerate(haz_rcp_epochs):
            damages = [df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch)] for df in damage_results]
            damages = pd.concat(damages, axis=0, ignore_index=True)
            print("* Done with concatinating all dataframes")
            damages.drop("confidence", axis=1, inplace=True)

            index_columns = [c for c in damages.columns.values.tolist() if ("EAD_" not in c) and ("EAEL_" not in c)]
            index_columns = [i for i in index_columns if i not in uncertainty_columns]
            damage_columns = [c for c in damages.columns.values.tolist() if ("EAD_" in c) or ("EAEL_" in c)]

            if len(damages.index) > 0:
                summarised_damages.append(quantiles(damages, index_columns, damage_columns))
        summarised_damages = pd.concat(summarised_damages, axis=0, ignore_index=True)
        summarised_damages.to_csv(
            os.path.join(
                summary_results,
                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL.csv",
            ),
            index=False,
        )
        del summarised_damages


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    tqdm.pandas()
    collapse_sensitivity()
