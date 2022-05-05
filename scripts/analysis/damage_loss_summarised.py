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

def quantiles(dataframe,grouping_by_columns,grouped_columns):
    grouped = dataframe.groupby(grouping_by_columns,dropna=False)[grouped_columns].agg([np.min, np.mean, np.max]).reset_index()
    grouped.columns = grouping_by_columns + [f"{prefix}_{agg_name}" for prefix, agg_name in grouped.columns if prefix not in grouping_by_columns]

    return grouped

def main(config,direct_damages_folder,
        summary_results_folder,
        network_csv,
        parameter_txt_file):

    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']

    direct_damages_results = os.path.join(output_data_path,direct_damages_folder)

    summary_results = os.path.join(output_data_path,summary_results_folder)
    if os.path.exists(summary_results) == False:
        os.mkdir(summary_results)

    asset_data_details = pd.read_csv(network_csv)
    param_values = open(parameter_txt_file)
    param_values = [tuple(line.split(',')) for line in param_values.readlines()]
    param_values = pd.DataFrame(param_values,columns=["parameter_set","cost_uncertainty_parameter","damage_uncertainty_parameter"])

    uncertainty_columns = ["cost_uncertainty_parameter","damage_uncertainty_parameter"]
    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column
        asset_damages_results = os.path.join(direct_damages_results,f"{asset_info.asset_gpkg}_{asset_info.asset_layer}")

        # Process the exposure and damage results
        damage_files = [os.path.join(
                                asset_damages_results,
                                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_direct_damages_parameter_set_{param.parameter_set}.parquet"
                                ) for param in param_values.itertuples()]
        damage_results = [pd.read_parquet(file) for file in damage_files if os.path.isfile(file) is True]
        # print ("* Done with creating list of all dataframes")

        if damage_results:
            exposures = damage_results[0].copy()
            hazard_columns = [c for c in exposures.columns.values.tolist() if c not in [asset_id,
                                                                                    'exposure_unit',
                                                                                    'damage_cost_unit',
                                                                                    'damage_uncertainty_parameter',
                                                                                    'cost_uncertainty_parameter',
                                                                                    'exposure']]
            exposures[hazard_columns] = exposures["exposure"].to_numpy()[:,None]*np.where(exposures[hazard_columns]>0,1,0)

            sum_dict = dict([(hk,"sum") for hk in hazard_columns])
            exposures = exposures.groupby([asset_id,
                                        'exposure_unit'
                                        ],
                                        dropna=False).agg(sum_dict).reset_index()
            exposures.to_parquet(os.path.join(summary_results,
                            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_exposures.parquet"),index=False)
            del exposures

            damages = []
            losses = []
            for df in damage_results:
                df = df.groupby([asset_id,
                                'damage_cost_unit',
                                ],
                                dropna=False).agg(sum_dict).reset_index()
                damages.append(df)
                if asset_info.single_failure_scenarios != "none":
                    loss_column = ["economic_loss"]
                    if asset_info.sector != "buildings":
                        loss_df = pd.read_csv(os.path.join(output_data_path,asset_info.single_failure_scenarios))
                        if asset_info.asset_gpkg == "potable_facilities_NWC":
                            loss_df[asset_id] = loss_df.progress_apply(
                                lambda x: str(x[asset_id]).lower().replace(" ","_").replace(".0",""),
                                axis=1)
                    else:
                        loss_df = gpd.read_file(os.path.join(processed_data_path,asset_info.single_failure_scenarios),layer="areas")
                        loss_df.rename(columns={"total_GDP":"economic_loss"},inplace=True)
                    df = pd.merge(df,loss_df[[asset_info.asset_id_column,"economic_loss"]],
                                    how="left",on=[asset_info.asset_id_column]).fillna(0)
                    df["economic_loss_unit"] = "J$/day"
                    del loss_df
                    loss = df.copy()
                    loss[hazard_columns] = loss["economic_loss"].to_numpy()[:,None]*np.where(loss[hazard_columns]>0,1,0)
                    losses.append(loss[[asset_id,"economic_loss_unit"] + hazard_columns])

            damages = pd.concat(damages,axis=0,ignore_index=True)
            print ("* Done with concatinating all dataframes")
            if len(damages.index) > 0:
                damages = quantiles(damages,[asset_id,'damage_cost_unit'],hazard_columns)
                damages.to_parquet(os.path.join(summary_results,
                            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_damages.parquet"),index=False)
            del damages

            if len(losses) > 0:
                losses = pd.concat(losses,axis=0,ignore_index=True)
                print ("* Done with concatinating all dataframes")
                if len(losses.index) > 0:
                    losses = quantiles(losses,[asset_id,'economic_loss_unit'],hazard_columns)
                    losses.to_parquet(os.path.join(summary_results,
                                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_losses.parquet"),index=False)
            del losses

            # damages = [df.groupby([asset_id,
            #                     'damage_cost_unit',
            #                     ],
            #                     dropna=False).agg(sum_dict).reset_index() for df in damage_results]

            # damages = pd.concat(damages,axis=0,ignore_index=True)
            # # print ("* Done with concatinating all dataframes")
            # if len(damages.index) > 0:
            #     damages = quantiles(damages,[asset_id,'damage_cost_unit'],hazard_columns)
            #     damages.to_parquet(os.path.join(summary_results,
            #                 f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_damages.parquet"),index=False)
            # del damages
        # Process the EAD and EAEL results
        damage_files = [os.path.join(
                                asset_damages_results,
                                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_parameter_set_{param.parameter_set}.csv"
                                ) for param in param_values.itertuples()]
        damage_results = [pd.read_csv(file) for file in damage_files if os.path.isfile(file) is True]

        if damage_results:
            # print ([len(df.index) for df in damage_results])
            for df in damage_results:
                df["rcp"] = df["rcp"].astype(str)
                df["epoch"] = df["epoch"].astype(str)
            haz_rcp_epochs = list(set(damage_results[0].set_index(["hazard","rcp","epoch"]).index.values.tolist()))
            # print (haz_rcp_epochs)
            summarised_damages = []
            for i,(haz,rcp,epoch) in enumerate(haz_rcp_epochs):
                damages = [df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch)] for df in damage_results]
                damages = pd.concat(damages,axis=0,ignore_index=True)
                # print ("* Done with concatinating all dataframes")
                damages.drop("confidence",axis=1,inplace=True)

                index_columns = [c for c in damages.columns.values.tolist() if ("EAD_" not in c) and ("EAEL_" not in c)]
                index_columns = [i for i in index_columns if i not in uncertainty_columns]
                damage_columns = [c for c in damages.columns.values.tolist() if ("EAD_" in c) or ("EAEL_" in c)]

                if len(damages.index) > 0:
                    summarised_damages.append(quantiles(damages,index_columns,damage_columns))
            summarised_damages = pd.concat(summarised_damages,axis=0,ignore_index=True)
            summarised_damages.to_csv(os.path.join(summary_results,
                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL.csv"),index=False)
            # print (len(summarised_damages.index))
            del summarised_damages
        print (f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer}")

if __name__ == "__main__":
    CONFIG = load_config()
    try:
        direct_damages_folder = str(sys.argv[1])
        summary_results_folder = str(sys.argv[2])
        network_csv = str(sys.argv[3])
        parameter_txt_file = str(sys.argv[4])
    except IndexError:
        print("Got arguments", sys.argv)
        exit()

    main(CONFIG,direct_damages_folder,
        summary_results_folder,
        network_csv,
        parameter_txt_file)
