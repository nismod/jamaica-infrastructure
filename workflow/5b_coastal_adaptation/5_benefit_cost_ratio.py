"""
Calculate Benefit Cost Ratio (BCR) of coastal defence adaptation options for a given asset class.
"""

import logging
import math
import os
import sys
import warnings

import click
import pandas as pd
from tqdm import tqdm

from jamaica_infrastructure.adaptation import (
    get_risks,
    get_benefits,
    get_ead_eael,
    get_ead_eael_benefits,
    bcr_estimates,
    write_empty,
)

warnings.simplefilter(action="ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

tqdm.pandas()


def scale_coastal_adaptation_costs(
    protection_asset_breakdown,
    option_cost_df,
    hazard_thresholds_column_name,
    asset_id,
    protection_asset_dict,
    rcp,
    rp,
    epoch,
):
    protection_asset_breakdown = pd.read_csv(protection_asset_breakdown)
    protect_dict = pd.read_parquet(protection_asset_dict)
    for index, row in option_cost_df.iterrows():
        f_id = row.get("flood_id")
        if pd.isna(f_id):
            option_cost_df.loc[index, "adapt_cost_npv"] = 0
            continue

        match = protection_asset_breakdown[
            (protection_asset_breakdown["polygon_id"] == f_id)
            & (protection_asset_breakdown["rcp"] == rcp)
            & (protection_asset_breakdown["epoch"] == epoch)
            & (protection_asset_breakdown["rp"] == rp)
        ]

        if match.empty:
            continue

        f_info = match.iloc[0]
        tot_cost = f_info["total_cost"]
        scaled_cost = row["adapt_cost_npv"]
        k = 0

        protect_asset = protect_dict[protect_dict[asset_id] == row[asset_id]]
        if protect_asset.empty or tot_cost == 0:
            final_cost = scaled_cost
        else:
            k = protect_asset["mean_cost"].values[0] / tot_cost
            final_cost = scaled_cost * k

        option_cost_df.loc[index, "adapt_cost_npv"] = final_cost * row[hazard_thresholds_column_name]

    return option_cost_df


def assign_coastal_protection_ft(
    option_cost_df, protection_asset_dict, hazard_thresholds_column_name, asset_id, rcp, rp, epoch
):
    df = pd.read_parquet(protection_asset_dict)

    # Create column names based on parameters
    rcp_str = str(int(rcp * 10))
    flood_height_col = f"flood_height_rcp_{rcp_str}{epoch}_rp_{rp}"
    flood_id_col = f"flood_id_rcp_{rcp_str}{epoch}_rp_{rp}"

    # Select and rename columns
    columns_to_keep = [asset_id, "mean_cost", flood_height_col, flood_id_col]
    df = df[columns_to_keep].rename(columns={flood_height_col: hazard_thresholds_column_name, flood_id_col: "flood_id"})

    # Fill and round flood heights
    df[hazard_thresholds_column_name] = df[hazard_thresholds_column_name].fillna(0)
    df[hazard_thresholds_column_name] = df[hazard_thresholds_column_name].map(lambda x: math.ceil(x * 2) / 2)

    # Merge with option_cost_df
    merged_df = option_cost_df.merge(df, on=asset_id, how="left")

    # Fill missing flood height with 0 and convert flood_id to nullable integer
    merged_df[hazard_thresholds_column_name] = merged_df[hazard_thresholds_column_name].fillna(0)
    merged_df["flood_id"] = pd.to_numeric(merged_df["flood_id"], errors="coerce").astype("Int64")

    return merged_df


def get_bcr_values(
    results_path,
    asset_id,
    bcr_results,
    risk_filepath,
    hazard,
    rcps,
    risk_type,
    val_type,
    no_adapt_risk_df,
    risk_columns,
    option_df,
    hazard_thresholds_column_name,
    protection_type_name,
    protection_asset_dict,
    flood_params,
    protection_asset_breakdown,
    days=10,
):
    folder_name = protection_type_name

    risk_file = os.path.join(
        results_path,
        folder_name,
        risk_filepath,
    )

    if os.path.isfile(risk_file) is True:
        adapt_risk_df = pd.read_csv(risk_file)
        adapt_risk_df, adapt_risk_columns = get_risks(
            adapt_risk_df,
            asset_id,
            hazard["hazard"],
            hazard["hazard_type"],
            rcps,
            risk_type,
            val_type,
            days=days,
        )
        risk_df, adapt_benefit_columns = get_benefits(
            asset_id,
            no_adapt_risk_df.copy(),
            adapt_risk_df,
            risk_columns,
            adapt_risk_columns,
        )
        option_cost_df = option_df.copy()
        rcp = flood_params[0]
        rp = flood_params[1]
        epoch = flood_params[2]

        option_cost_df = assign_coastal_protection_ft(
            option_cost_df, protection_asset_dict, hazard_thresholds_column_name, asset_id, rcp, rp, epoch
        )
        option_cost_df = scale_coastal_adaptation_costs(
            protection_asset_breakdown,
            option_cost_df,
            hazard_thresholds_column_name,
            asset_id,
            protection_asset_dict,
            rcp,
            rp,
            epoch,
        )
        risk_df, bcr_columns = bcr_estimates(
            asset_id,
            option_cost_df,
            risk_df,
            hazard_thresholds_column_name,
            adapt_benefit_columns,
        )
        risk_df = risk_df[
            [
                asset_id,
                "adaptation_option",
                hazard_thresholds_column_name,
                "cost_units",
                "adapt_cost_npv",
            ]
            + adapt_benefit_columns
            + bcr_columns
        ]

        bcr_results.append(risk_df)

    return bcr_results


def get_ead_eael_costs(
    results_path,
    asset_id,
    ead_eael_results,
    risk_filepath,
    hazard,
    rcps,
    risk_type,
    val_type,
    no_adapt_ead_eael_df,
    ead_eael_columns,
    option_df,
    hazard_thresholds_column_name,
    protection_type_name,
    protection_asset_dict,
    flood_params,
    protection_asset_breakdown,
):
    folder_name = protection_type_name
    risk_file = os.path.join(
        results_path,
        folder_name,
        risk_filepath,
    )
    if os.path.isfile(risk_file) is True:
        adapt_risk = pd.read_csv(risk_file)
        adapt_ead_eael_df, adapt_ead_eael_columns = get_ead_eael(
            adapt_risk,
            asset_id,
            hazard["hazard"],
            hazard["hazard_type"],
            rcps,
            risk_type,
            val_type,
        )
        # print (no_adapt_ead_eael_df.columns.values.tolist())
        ead_eael_df, ead_eael_benefit_columns = get_ead_eael_benefits(
            asset_id,
            no_adapt_ead_eael_df.copy(),
            adapt_ead_eael_df,
            ead_eael_columns,
            adapt_ead_eael_columns,
        )
        option_cost_df = option_df.copy()
        rcp, rp, epoch = flood_params

        option_cost_df = assign_coastal_protection_ft(
            option_cost_df, protection_asset_dict, hazard_thresholds_column_name, asset_id, rcp, rp, epoch
        )
        option_cost_df = scale_coastal_adaptation_costs(
            protection_asset_breakdown,
            option_cost_df,
            hazard_thresholds_column_name,
            asset_id,
            protection_asset_dict,
            rcp,
            rp,
            epoch,
        )

        option_cost_df["adapt_cost_unit"] = "J$"
        option_cost_df["ead_cost_unit"] = "J$"
        option_cost_df["eael_cost_unit"] = "J$/day"
        ead_eael_df = pd.merge(
            option_cost_df[
                [
                    asset_id,
                    "adaptation_option",
                    hazard_thresholds_column_name,
                    "adapt_cost_unit",
                    "ead_cost_unit",
                    "eael_cost_unit",
                    "adapt_cost_npv",
                ]
            ],
            ead_eael_df,
            how="left",
            on=[asset_id],
        ).fillna(0)
        ead_eael_df = ead_eael_df[ead_eael_df["adapt_cost_npv"] > 0]
        ead_eael_df = ead_eael_df[
            [
                asset_id,
                "adaptation_option",
                hazard_thresholds_column_name,
                "adapt_cost_unit",
                "ead_cost_unit",
                "eael_cost_unit",
                "adapt_cost_npv",
            ]
            + ead_eael_benefit_columns
        ]

        ead_eael_results.append(ead_eael_df)
    return ead_eael_results


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the asset layer data",
)
@click.option(
    "--cost-file",
    "-c",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to adaptation cost data",
)
@click.option(
    "--no-adapt-risk-file",
    "-r",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to NPV of EAD & EAEL with no adaptation",
)
@click.option(
    "--protection-asset-breakdown",
    "-pb",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to directroy with breakdown of assets for each flood polygon",
)
@click.option(
    "--protection-asset-dict",
    "-pa",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to directroy with network assets to flood portectoin area relatoinal dictionary",
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
    "--proj-end-year",
    "-py",
    required=True,
    type=int,
    help="Projection End Year",
)
@click.option(
    "--rcp",
    "-rcp",
    required=True,
    type=float,
    help="RPS value",
)
@click.option(
    "--rp",
    "-rp",
    required=True,
    type=int,
    help="RP value",
)
@click.option(
    "--disruption-duration-days",
    "-d",
    "days",
    required=True,
    type=float,
    help="Assumed duration in days of any wider economic loss.",
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=False, readable=True),
    help="Path to the output directory",
)
def benefit_cost_ratio(
    network_csv,
    cost_file,
    no_adapt_risk_file,
    protection_asset_breakdown,
    protection_asset_dict,
    asset_gpkg,
    asset_layer,
    proj_end_year,
    rcp,
    rp,
    days,
    output_dir,
):
    if not os.path.isfile(cost_file):
        raise FileNotFoundError(f"Cost file {cost_file} does not exist.")

    hazard_label = "coastal"
    hazard = {"hazard": f"{hazard_label}", "hazard_type": [f"{hazard_label}"]}
    rcps = [2.6, 4.5, 8.5]
    risk_type = ["EAD", "EAEL"]
    val_type = ["amin", "mean", "amax"]

    risk_filepath = os.path.join("loss_damage_npvs", f"{asset_gpkg}_{asset_layer}_EAD_EAEL_npvs.csv")

    if not os.path.isfile(no_adapt_risk_file):
        raise FileNotFoundError(f"Risk file {no_adapt_risk_file} does not exist.")

    file_prefix = f"{hazard_label}_{asset_gpkg}_{asset_layer}"
    output_csv_path = os.path.join(
        output_dir,
        "adaptation_benefits_costs_bcr",
    )
    output_bcr = os.path.join(
        output_csv_path,
        f"{file_prefix}_adaptation_benefits_costs_bcr.csv",
    )
    output_ead = os.path.join(
        output_csv_path,
        f"{file_prefix}_adaptation_costs_avoided_EAD_EAEL.csv",
    )

    asset_df = pd.read_csv(network_csv)
    asset_data_details = asset_df[(asset_df["asset_gpkg"] == asset_gpkg) & (asset_df["asset_layer"] == asset_layer)]
    if len(asset_data_details) > 1:
        raise ValueError((f"Multiple assets found for gpkg={asset_gpkg} " f"and layer={asset_layer}"))
    elif len(asset_data_details) == 0:
        raise ValueError(f"No asset found for gpkg={asset_gpkg} and layer={asset_layer}")
    asset_info = asset_data_details.squeeze()

    asset_id = asset_info.asset_id_column

    logging.info(f"{hazard_label} {asset_gpkg} {asset_layer}")

    logging.info("Reading costs")
    cost_df = pd.read_csv(cost_file)

    if cost_df.empty:
        logging.info("No adaptation options (costs) available, skipping...")
        write_empty(asset_id, output_bcr, output_ead)
        return

    logging.info("Reading risks with no adaptation")
    no_adapt_risk = pd.read_csv(no_adapt_risk_file)

    logging.info("Calculate BCR")
    adaptation_options = list(set(cost_df["adaptation_option"].values.tolist()))
    no_adapt_risk_df, risk_columns = get_risks(
        no_adapt_risk,
        asset_id,
        hazard_label,
        hazard["hazard_type"],
        rcps,
        risk_type,
        val_type,
        days=days,
    )
    no_adapt_ead_eael_df, ead_eael_columns = get_ead_eael(
        no_adapt_risk,
        asset_id,
        hazard_label,
        hazard["hazard_type"],
        rcps,
        risk_type,
        val_type,
    )
    bcr_results = []
    ead_eael_results = []
    flood_params = [rcp, rp, proj_end_year]

    for option in adaptation_options:
        logging.info(f"{option=}")
        option_df = cost_df[cost_df["adaptation_option"] == option]
        asset_adaptation_cost = option_df["asset_adaptation_cost"].values[0]
        if asset_adaptation_cost == "J$/m":
            bcr_results = get_bcr_values(
                output_dir,
                asset_id,
                bcr_results,
                risk_filepath,
                hazard,
                rcps,
                risk_type,
                val_type,
                no_adapt_risk_df,
                risk_columns,
                option_df,
                "flood_depth_protection_level",
                "coastal_adaptation",
                protection_asset_dict,
                flood_params,
                protection_asset_breakdown,
                days=days,
            )
            ead_eael_results = get_ead_eael_costs(
                output_dir,
                asset_id,
                ead_eael_results,
                risk_filepath,
                hazard,
                rcps,
                risk_type,
                val_type,
                no_adapt_ead_eael_df,
                ead_eael_columns,
                option_df,
                "flood_depth_protection_level",
                "coastal_adaptation",
                protection_asset_dict,
                flood_params,
                protection_asset_breakdown,
            )
        else:
            logging.error("This script may have been called on the wrong data. Quitting.")
            sys.exit(1)

    if len(bcr_results) > 0:
        bcr_results = pd.concat(bcr_results, axis=0, ignore_index=False)
    pd.DataFrame(bcr_results).to_csv(output_bcr, index=False)
    logging.info(f"Writing BCR results to disk: {output_bcr}")

    if len(ead_eael_results) > 0:
        ead_eael_results = pd.concat(ead_eael_results, axis=0, ignore_index=False)
    pd.DataFrame(ead_eael_results).to_csv(output_ead, index=False)
    logging.info(f"Writing EAD & EAEL results to disk: {output_ead}")


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )
    benefit_cost_ratio()
