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
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
    epsg_jamaica = 3448

    direct_damages_results = os.path.join(output_data_path,"direct_damages")

    asset_data_details = pd.read_csv(os.path.join(processed_data_path,
                        "networks",
                        "network_layers_hazard_intersections_details.csv"))
    bridge_flood_protection = 50 # Bridges in Jamaica are designed to withstand 1 in 50 year floods

    for asset_info in asset_data_details.itertuples():
        asset_damages_results = os.path.join(direct_damages_results,f"{asset_info.asset_gpkg}_{asset_info.asset_layer}")
        param_values = open('parameter_combinations.txt', "r")
        for param in param_values:
            parameter_set = param.split(",")[0]
            damage_file = os.path.join(
                                asset_damages_results,
                                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_direct_damages_parameter_set_{parameter_set}.csv"
                                )
            if os.path.isfile(damage_file) is True:
                expected_damages = []
                df = pd.read_csv(damage_file).fillna(0)
                if asset_info.single_failure_scenarios != "none":
                    if asset_info.sector != "buildings":
                        loss_df = pd.read_csv(os.path.join(output_data_path,asset_info.single_failure_scenarios))
                    else:
                        loss_df = gpd.read_file(os.path.join(processed_data_path,asset_info.single_failure_scenarios),layer="areas")
                        loss_df.rename(columns={"total_GDP":"economic_loss"},inplace=True)
                        loss_df["loss_unit"] = "JD/day"
                    df = pd.merge(df,loss_df[[asset_info.asset_id_column,"economic_loss","loss_unit"]],
                                    how="left",on=[asset_info.asset_id_column]).fillna(0)
                haz_rcp_epoch_confidence = list(set(df.set_index(["hazard","rcp","epoch","confidence"]).index.values.tolist()))
                # print (haz_rcp_epoch_confidence)
                for i,(haz,rcp,epoch,confidence) in enumerate(haz_rcp_epoch_confidence):
                    damages = df[(df.hazard == haz) & (df.rcp == rcp) & (df.epoch == epoch) & (df.confidence == confidence)]
                    # print (damages)
                    damages['probability'] = 1.0/damages['rp']
                    index_columns = [c for c in damages.columns.values.tolist() if c not in [
                                                                        'rp',
                                                                        'probability',
                                                                        'direct_damage_cost',
                                                                        'economic_loss',
                                                                        'exposure']
                                    ]
                    expected_damage_df = risks_pivot(damages,index_columns,'probability',
                                                'direct_damage_cost',None,'EAD',
                                                flood_protection=None)
                    # print (expected_damage_df)
                    if 'economic_loss' in damages.columns.values.tolist():
                        economic_loss_df = risks_pivot(damages,index_columns,'probability',
                                                'economic_loss',None,'EAEL',
                                                flood_protection=None)
                        expected_damage_df = pd.merge(expected_damage_df,economic_loss_df,how='left',on=index_columns).fillna(0)
                        del economic_loss_df
                    if (asset_info.asset_gpkg == "roads") and (asset_info.asset_layer == "nodes") and (haz in ["coastal","fluvial","surface"]):
                        damages["protection_standard"] = bridge_flood_protection
                        protected_damage_df = risks_pivot(damages,index_columns + ["protection_standard"],'probability',
                                                'direct_damage_cost',"protection_standard",'EAD',
                                                flood_protection="yes",flood_protection_name="designed_protection")
                        protected_loss_df = risks_pivot(damages,index_columns + ["protection_standard"],'probability',
                                                'economic_loss',"protection_standard",'EAEL',
                                                flood_protection="yes",flood_protection_name="designed_protection")

                        expected_damage_df = pd.merge(expected_damage_df,protected_damage_df,how='left',on=index_columns).fillna(0)
                        expected_damage_df = pd.merge(expected_damage_df,protected_loss_df,how='left',on=index_columns).fillna(0)
                        del protected_damage_df, protected_loss_df
                    expected_damages.append(expected_damage_df)
                    del expected_damage_df

                expected_damages = pd.concat(expected_damages,axis=0,ignore_index=True)
                asset_damages_results = os.path.join(direct_damages_results,f"{asset_info.asset_gpkg}_{asset_info.asset_layer}")
                if os.path.exists(asset_damages_results) == False:
                    os.mkdir(asset_damages_results)
                expected_damages.to_csv(os.path.join(asset_damages_results,
                            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_parameter_set_{parameter_set}.csv"),
                            index=False)
        param_values.close()
        print (f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer}")
                



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)