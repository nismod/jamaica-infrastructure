"""Estimate direct damages to physical assets exposed to hazards

"""

import sys
import os

import pandas as pd
import geopandas as gpd
from shapely import wkb
import numpy as np
from analysis_utils import *
from tqdm import tqdm

tqdm.pandas()


def quantiles(dataframe, grouping_by_columns, grouped_columns):
    # quantiles_list = ['mean','min','max','median','q5','q95']
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


# def quantiles(dataframe,grouping_by_columns,grouped_columns):
#     # quantiles_list = ['mean','min','max','median','q5','q95']
#     quantiles_list = ['mean','min','max']
#     groouped = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].agg([np.mean, np.min, np.max])
#     grouped.columns = [f"{prefix}_{agg_name}" for prefix, agg_name in grouped.columns]
#     df_list = []
#     for quant in quantiles_list:
#         if quant == 'mean':
#             # print (dataframe)
#             df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].mean()
#         elif quant == 'min':
#             df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].min()
#         elif quant == 'max':
#             df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].max()
#         # elif quant == 'median':
#         #     df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].quantile(0.5)
#         # elif quant == 'q5':
#         #     df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].quantile(0.05)
#         # elif quant == 'q95':
#         #     df = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].quantile(0.95)

#         df.rename(columns=dict((g,'{}_{}'.format(g,quant)) for g in grouped_columns),inplace=True)
#         df_list.append(df)
#     return pd.concat(df_list,axis=1).reset_index()


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    direct_damages_results = os.path.join(output_data_path, "direct_damages")

    summary_results = os.path.join(output_data_path, "direct_damages_summary")
    if os.path.exists(summary_results) == False:
        os.mkdir(summary_results)

    asset_data_details = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks",
            "network_layers_hazard_intersections_details_1.csv",
        )
    )
    # asset_data_details = asset_data_details[(asset_data_details["sector"] == "energy") & (asset_data_details["asset_layer"] == "edges")]
    # print (asset_data_details)

    param_values = open("parameter_combinations.txt")
    param_values = [tuple(line.split(",")) for line in param_values.readlines()]
    param_values = pd.DataFrame(
        param_values,
        columns=[
            "parameter_set",
            "cost_uncertainty_parameter",
            "damage_uncertainty_parameter",
        ],
    )

    uncertainty_columns = ["cost_uncertainty_parameter", "damage_uncertainty_parameter"]
    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column
        asset_damages_results = os.path.join(
            direct_damages_results, f"{asset_info.asset_gpkg}_{asset_info.asset_layer}"
        )

        # Process the exposure and damage results
        damage_files = [
            os.path.join(
                asset_damages_results,
                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_direct_damages_parameter_set_{param.parameter_set}.parquet",
            )
            for param in param_values.itertuples()
        ]
        damage_results = [
            pd.read_parquet(file)
            for file in damage_files
            if os.path.isfile(file) is True
        ]
        print("* Done with creating list of all dataframes")

        if damage_results:
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
            exposures[hazard_columns] = exposures["exposure"].to_numpy()[
                :, None
            ] * np.where(exposures[hazard_columns] > 0, 1, 0)

            sum_dict = dict([(hk, "sum") for hk in hazard_columns])
            exposures = (
                exposures.groupby([asset_id, "exposure_unit"], dropna=False)
                .agg(sum_dict)
                .reset_index()
            )
            # exposures.to_parquet(os.path.join(summary_results,
            #                 f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_exposures.parquet"),index=False)
            exposures.to_csv(
                os.path.join(
                    summary_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_exposures.csv",
                ),
                index=False,
            )
            del exposures
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
            print("* Done with concatinating all dataframes")
            if len(damages.index) > 0:
                damages = quantiles(
                    damages, [asset_id, "damage_cost_unit"], hazard_columns
                )
                # damages.to_parquet(os.path.join(summary_results,
                #             f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_damages.parquet"),index=False)
                damages.to_csv(
                    os.path.join(
                        summary_results,
                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_damages.csv",
                    ),
                    index=False,
                )
            del damages
        # Process the EAD and EAEL results
        damage_files = [
            os.path.join(
                asset_damages_results,
                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_parameter_set_{param.parameter_set}.csv",
            )
            for param in param_values.itertuples()
        ]
        damage_results = [
            pd.read_csv(file) for file in damage_files if os.path.isfile(file) is True
        ]

        if damage_results:
            print([len(df.index) for df in damage_results])
            for df in damage_results:
                df["rcp"] = df["rcp"].astype(str)
                df["epoch"] = df["epoch"].astype(str)
            haz_rcp_epochs = list(
                set(
                    damage_results[0]
                    .set_index(["hazard", "rcp", "epoch"])
                    .index.values.tolist()
                )
            )
            print(haz_rcp_epochs)
            summarised_damages = []
            for i, (haz, rcp, epoch) in enumerate(haz_rcp_epochs):
                damages = [
                    df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch)]
                    for df in damage_results
                ]
                damages = pd.concat(damages, axis=0, ignore_index=True)
                print("* Done with concatinating all dataframes")
                damages.drop("confidence", axis=1, inplace=True)

                index_columns = [
                    c
                    for c in damages.columns.values.tolist()
                    if ("EAD_" not in c) and ("EAEL_" not in c)
                ]
                index_columns = [
                    i for i in index_columns if i not in uncertainty_columns
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
                os.path.join(
                    summary_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL.csv",
                ),
                index=False,
            )
            print(len(summarised_damages.index))
            del summarised_damages
        print(f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer}")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
