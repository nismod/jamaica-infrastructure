"""Estimate adaptation options costs and benefits

"""
import sys
import os

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
import geopandas as gpd
import numpy as np
from analysis_utils import *
from tqdm import tqdm
tqdm.pandas()
     
def get_benefits(id_column,df1,df2,df1_columns,df2_columns):
    modified_columns = dict([(c,c.replace("risk","protection_risk")) for c in df2_columns])
    df2.replace(columns=modified_columns,inplace=True)
    benefit_columns = [c.replace("risk","avoided_risk") for c in df1_columns]

    df1 = pd.merge(df1,df2,how="left",on=[id_column]).fillna(0)
    for i, (b,c) in enumerate(list(zip(benefit_columns,df1_columns))):
        if c in df2_columns:
            df1[b] = df1[c] - df1[c.replace("risk","protection_risk")]
        else:
            df1[b] = df1[c]

    return df1, benefit_columns


def get_all_column_combinations(hzd,rcps,risk_type,val_type):
    all_c = [] 
    rcp_c = []
    for rcp in rcps:
        r_c = []
        for h in hzd:
            for rt in risk_type:
                for vt in val_type:
                    all_c.append(f"{h}__rcp_{rcp}__{rt}_{vt}")
                    r_c.append(f"{h}__rcp_{rcp}__{rt}_{vt}")
        rcp_c.append(r_c)
    return all_c,rcp_c

def get_risks(df,asset_id,hazard,hazard_types,rcps,risk_type,val_type,days=10):
    all_columns, rcp_columns = get_all_column_combinations(hazard_types,rcps,risk_type,val_type)
    all_columns = [c for c in df.columns.values.tolist() if c in all_columns]
    eael_columns = [c for c in all_columns if "EAEL_" in c]
    if len(eael_columns) > 0:
        df[eael_columns] = days*df[eael_columns]

    risk_rcp_columns = []
    for ri,(rcp,rcp_c) in enumerate(list(zip(rcps,rcp_columns))):
        rcp_c = [c for c in rcp_c if c in df.columns.values.tolist()]
        if len(rcp_c) > 0:
            for vt in val_type:
                vt_cols = [c for c in rcp_c if f"_{vt}" in c]
                risk_rcp_columns.append(f"{hazard}__rcp_{rcp}__risk_{vt}")
                df[f"{hazard}__rcp_{rcp}__risk_{vt}"] = df[vt_cols].sum(axis=1)

    return df[[asset_id] + risk_rcp_columns], risk_rcp_columns

def bcr_estimates(asset_id,option_cost_df,risk_df,hazard_thresholds_column_name,adapt_benefit_columns):
    option_cost_df["cost_units"] = "J$"
    risk_df = pd.merge(option_cost_df[[asset_id,"adaptation_option",
                                    hazard_thresholds_column_name,"cost_units","adapt_cost_npv"]],
                                risk_df,how="left",on=[asset_id]).fillna(0)
    risk_df = risk_df[risk_df["adapt_cost_npv"] > 0]
    bcr_columns = [c.replace("avoided_risk","BCR") for c in adapt_benefit_columns]
    risk_df[bcr_columns] = risk_df[adapt_benefit_columns].div(risk_df["adapt_cost_npv"],axis=0)

    return risk_df

def get_bcr_values(results_path,asset_id,bcr_results,asset_info,
                    hazard,rcps,risk_type,val_type,
                    no_adapt_risk_df,risk_columns,
                    option_df,
                    hazard_thresholds,cost_multiplication_factors,hazard_thresholds_column_name,
                    protection_type_name,days=10):
    print (option_df)
    for idx, (ft,cmf) in enumerate(list(zip(hazard_thresholds,cost_multiplication_factors))):
        folder_name = f"{protection_type_name}_{str(ft).replace('.','p')}"
        results_folder = os.path.join(results_path,folder_name)
        risk_file = os.path.join(results_folder,"loss_damage_npvs",
                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_npvs.csv")
        if os.path.isfile(risk_file) is True:
            adapt_risk_df = pd.read_csv(risk_file)
            print ("Adapt")
            print (adapt_risk_df)
            adapt_risk_df, adapt_risk_columns = get_risks(adapt_risk_df,asset_id,
                                                hazard['hazard'],hazard['hazard_type'],
                                                rcps,risk_type,val_type,days=days)
            risk_df, adapt_benefit_columns = get_benefits(asset_id,no_adapt_risk_df.copy(),adapt_risk_df,
                                                                risk_columns,adapt_risk_columns)
            
            option_cost_df = option_df.copy()
            option_cost_df[hazard_thresholds_column_name] = ft
            option_cost_df["adapt_cost_npv"] = cmf*option_cost_df["adapt_cost_npv"]
            risk_df = bcr_estimates(asset_id,option_cost_df,risk_df,hazard_thresholds_column_name,adapt_benefit_columns)
            risk_df = risk_df[[asset_id,
                            "adaptation_option",
                            hazard_thresholds_column_name,
                            "cost_units",
                            "adapt_cost_npv"] + adapt_benefit_columns,bcr_columns]

            bcr_results.append(risk_df)

    return bcr_results

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    output_data_path = config['paths']['output']
     
    days = 15
    adaptation_results = os.path.join(output_data_path,"adaptation_costs")
    adaptation_bcr_results = os.path.join(output_data_path,"adaptation_benefits_costs_bcr")
    if os.path.exists(adaptation_bcr_results) == False:
            os.mkdir(adaptation_bcr_results)

    non_adapt_risk_results = os.path.join(output_data_path,"loss_damage_npvs")

    asset_data_details = pd.read_csv(os.path.join(processed_data_path,
                        "networks",
                        "network_layers_hazard_intersections_details.csv"))
    adapt_hazards = [
                        {
                            "hazard":"flooding",
                            "hazard_type":["coastal","fluvial","surface"]
                        },
                        {
                            "hazard":"TC",
                            "hazard_type":["cyclone"]
                        }
                    ]
    rcps = [2.6,4.5,8.5]
    minimum_threshold = 0.5
    risk_type = ["EAD","EAEL"]
    val_type = ["amin","mean","amax"]
    for hazard in adapt_hazards:
        hazard_adapt_costs = os.path.join(adaptation_results,f"{hazard['hazard']}_costs")
        for asset_info in asset_data_details.itertuples():
            asset_id = asset_info.asset_id_column
            cost_file = os.path.join(hazard_adapt_costs,
                                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_adaptation_timeseries_and_npvs.csv")
            risk_file = os.path.join(non_adapt_risk_results,
                                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_npvs.csv")
            if (os.path.isfile(cost_file) is True) and (os.path.isfile(risk_file) is True):
                cost_df = pd.read_csv(cost_file)
                adaptation_options = list(set(cost_df["adaptation_option"].values.tolist()))
                no_adapt_risk_df = pd.read_csv(risk_file)
                no_adapt_risk_df, risk_columns = get_risks(no_adapt_risk_df,asset_id,hazard['hazard'],hazard['hazard_type'],rcps,risk_type,val_type,days=days)
                # print (no_adapt_risk_df)
                # print (risk_columns)
                bcr_results = []
                for option in adaptation_options:
                    option_df = cost_df[cost_df["adaptation_option"] == option]
                    asset_adaptation_cost = option_df["asset_adaptation_cost"].values[0]
                    # print (option_df)
                    if hazard['hazard'] == "flooding" and asset_adaptation_cost == "J$/m":
                        flood_thresholds = [1.0,1.5,2.0,2.5]
                        bcr_results = get_bcr_values(output_data_path,asset_id,bcr_results,asset_info,
                                                    hazard,rcps,risk_type,val_type,
                                                    no_adapt_risk_df,risk_columns,
                                                    option_df,
                                                    flood_thresholds,flood_thresholds,"flood_depth_protection_level",
                                                    "flood_threshold",days=days)
                    elif hazard['hazard'] == "TC" and asset_info.sector == "energy":
                        cyclone_damage_curve_change = [0.76]
                        bcr_results = get_bcr_values(output_data_path,asset_id,bcr_results,asset_info,
                                                    hazard,rcps,risk_type,val_type,
                                                    no_adapt_risk_df,risk_columns,
                                                    option_df,
                                                    cyclone_damage_curve_change,[1],"cyclone_damage_curve_reduction",
                                                    "cyclone_damage_curve_change",days=days)
                    else:
                        print (risk_columns)
                        option_df["flood_protection_level"] = 'All'
                        adapt_benefit_columns = [c.replace("risk","avoided_risk") for c in risk_columns]
                        bcr_results = bcr_estimates(asset_id,option_df,
                                                    no_adapt_risk_df,
                                                    "flood_protection_level",
                                                    adapt_benefit_columns)
                if len(bcr_results) > 0:
                    bcr_results = pd.concat(bcr_results,axis=0,ignore_index=False)

                    bcr_results.to_csv(adaptation_bcr_results,
                        f"{hazard}_{asset_info.asset_gpkg}_{asset_info.asset_layer}_adaptation_benefits_costs_bcr.csv")

                    print (f"* Done with {hazard} {asset_info.asset_gpkg} {asset_info.asset_layer}")



                



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)