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


def get_all_column_combinations(hzd, rcps, risk_type, val_type):
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
    return all_c, rcp_c


def get_risks(df, asset_id, hazard, hazard_types, rcps, risk_type, val_type, days=10):
    all_columns, rcp_columns = get_all_column_combinations(
        hazard_types, rcps, risk_type, val_type
    )
    all_columns = [c for c in df.columns.values.tolist() if c in all_columns]
    eael_columns = [c for c in all_columns if "EAEL_" in c]
    if len(eael_columns) > 0:
        df[eael_columns] = days * df[eael_columns]

    risk_rcp_columns = []
    for ri, (rcp, rcp_c) in enumerate(list(zip(rcps, rcp_columns))):
        rcp_c = [c for c in rcp_c if c in df.columns.values.tolist()]
        if len(rcp_c) > 0:
            for vt in val_type:
                vt_cols = [c for c in rcp_c if f"_{vt}" in c]
                risk_rcp_columns.append(f"{hazard}__rcp_{rcp}__risk_{vt}")
                df[f"{hazard}__rcp_{rcp}__risk_{vt}"] = df[vt_cols].sum(axis=1)

    return df[[asset_id] + risk_rcp_columns], risk_rcp_columns


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    discounted_results_nbs = os.path.join(
        output_data_path, "loss_damage_npvs_with_mangroves"
    )
    discounted_results = os.path.join(output_data_path, "loss_damage_npvs")

    nbs_benefits_results = os.path.join(
        output_data_path, "coastal_changes_with_mangroves"
    )
    if os.path.exists(nbs_benefits_results) == False:
        os.mkdir(nbs_benefits_results)

    asset_data_details = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks",
            "network_layers_hazard_intersections_details.csv",
        )
    )

    rcps = [2.6, 4.5, 8.5]
    risk_type = ["EAD", "EAEL"]
    val_type = ["amin", "mean", "amax"]
    days = 15
    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column
        file = os.path.join(
            discounted_results_nbs,
            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_npvs.csv",
        )
        index_columns = [asset_id, "damage_cost_unit", "economic_loss_unit"]
        if os.path.isfile(file) is True:
            nbs_damages = pd.read_csv(file)

            nbs_risks, nbs_risk_columns = get_risks(
                nbs_damages.copy(),
                asset_id,
                "coastal",
                ["coastal"],
                rcps,
                risk_type,
                val_type,
                days=days,
            )
            damages = pd.read_csv(
                os.path.join(
                    discounted_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_npvs.csv",
                )
            )
            damages = damages[
                damages[asset_id].isin(list(set(nbs_damages[asset_id].values.tolist())))
            ]
            damages = damages[nbs_damages.columns.values.tolist()]

            risks, risk_columns = get_risks(
                damages.copy(),
                asset_id,
                "coastal",
                ["coastal"],
                rcps,
                risk_type,
                val_type,
                days=days,
            )

            diff = damages.set_index(index_columns).subtract(
                nbs_damages.set_index(index_columns), fill_value=0
            )
            diff.reset_index().to_csv(
                os.path.join(
                    nbs_benefits_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_coastal_EAD_EAEL_change_with_mangroves.csv",
                ),
                index=False,
            )

            benefit_columns = [
                c.replace("risk", "avoided_risk") for c in nbs_risk_columns
            ]
            diff = risks.set_index([asset_id]).subtract(
                nbs_risks.set_index([asset_id]), fill_value=0
            )
            diff = diff.reset_index()

            diff.rename(
                columns=dict(zip(nbs_risk_columns, benefit_columns)), inplace=True
            )
            diff.to_csv(
                os.path.join(
                    nbs_benefits_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_coastal_risks_change_with_mangroves.csv",
                ),
                index=False,
            )

            gdf = gpd.read_file(
                os.path.join(processed_data_path, asset_info.path),
                layer=asset_info.asset_layer,
            )[[asset_id, "geometry"]]
            gpd.GeoDataFrame(
                pd.merge(diff, gdf, how="left", on=[asset_id]),
                geometry="geometry",
                crs=f"EPSG:{epsg_jamaica}",
            ).to_file(
                os.path.join(
                    nbs_benefits_results,
                    f"{asset_info.asset_gpkg}_coastal_risk_change_with_mangroves.gpkg",
                ),
                layer=asset_info.asset_layer,
                driver="GPKG",
            )
            print(
                f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer} values"
            )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
