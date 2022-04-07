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

    discounted_results_nbs = os.path.join(output_data_path,"loss_damage_npvs_with_mangroves")
    discounted_results = os.path.join(output_data_path,"loss_damage_npvs")

    nbs_benefits_results = os.path.join(output_data_path,"coastal_changes_with_mangroves")
    if os.path.exists(nbs_benefits_results) == False:
        os.mkdir(nbs_benefits_results)

    
    asset_data_details = pd.read_csv(os.path.join(processed_data_path,
                        "networks",
                        "network_layers_hazard_intersections_details.csv"))
    

    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column
        file = os.path.join(discounted_results_nbs,
                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_npvs.csv")
        index_columns = [asset_id,"damage_cost_unit","economic_loss_unit"]
        if os.path.isfile(file) is True:
            nbs_damages = pd.read_csv(file)

            damages = pd.read_csv(os.path.join(discounted_results,
                        f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_npvs.csv"))
            damages = damages[damages[asset_id].isin(list(set(nbs_damages[asset_id].values.tolist())))]
            damages = damages[nbs_damages.columns.values.tolist()] 
            diff = damages.set_index(index_columns).subtract(nbs_damages.set_index(index_columns),fill_value=0)
            
            diff.reset_index().to_csv(os.path.join(nbs_benefits_results,
                                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_coastal_risk_change_with_mangroves.csv"),index=False)

            gdf = gpd.read_file(os.path.join(processed_data_path,asset_info.path),layer=asset_info.asset_layer)[[asset_id,"geometry"]]
            gpd.GeoDataFrame(pd.merge(diff,gdf,how="left",on=[asset_id]),
                geometry="geometry",crs=f"EPSG:{epsg_jamaica}").to_file(os.path.join(nbs_benefits_results,
                                f"{asset_info.asset_gpkg}_coastal_risk_change_with_mangroves.gpkg"),layer=asset_info.asset_layer,driver="GPKG")
            print (f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer} values")
if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)