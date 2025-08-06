"""
Calculate Benefit Cost Ratio (BCR) of adaptation options for given asset class.
"""

import logging
import os
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


def ead_eael_estimates(
    asset_id,
    option_cost_df,
    risk_df,
    hazard_thresholds_column_name,
):
    option_cost_df["adapt_cost_unit"] = "J$"
    option_cost_df["ead_cost_unit"] = "J$"
    option_cost_df["eael_cost_unit"] = "J$/day"
    risk_df = pd.merge(
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
        risk_df,
        how="left",
        on=[asset_id],
    ).fillna(0)
    risk_df = risk_df[risk_df["adapt_cost_npv"] > 0]

    return risk_df


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
    hazard_thresholds,
    cost_multiplication_factors,
    hazard_thresholds_column_name,
    protection_type_name,
    days=10,
):
    for idx, (ft, cmf) in enumerate(list(zip(hazard_thresholds, cost_multiplication_factors))):
        folder_name = f"{protection_type_name}_{str(ft).replace('.','p')}"
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
            option_cost_df[hazard_thresholds_column_name] = ft
            option_cost_df["adapt_cost_npv"] = cmf * option_cost_df["adapt_cost_npv"]
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
    hazard_thresholds,
    cost_multiplication_factors,
    hazard_thresholds_column_name,
    protection_type_name,
):
    for idx, (ft, cmf) in enumerate(list(zip(hazard_thresholds, cost_multiplication_factors))):
        folder_name = f"{protection_type_name}_{str(ft).replace('.','p')}"
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
            option_cost_df[hazard_thresholds_column_name] = ft
            option_cost_df["adapt_cost_npv"] = cmf * option_cost_df["adapt_cost_npv"]
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
    "--hazard-label",
    "-h",
    required=True,
    help="hazard label",
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
    "--disruption-duration-days",
    "-d",
    "days",
    required=True,
    type=float,
    help="Assumed duration in days of any wider economic loss.",
)
@click.option(
    "--flood-defence-threshold",
    "-f",
    "flood_thresholds",
    required=True,
    type=float,
    multiple=True,
    help="Flood defence heights to consider in meters.",
)
@click.option(
    "--tc-damage-curve-factor",
    "-t",
    "cyclone_damage_curve_change",
    required=True,
    type=float,
    multiple=True,
    help="TC wind damage curve modifiers. See config.yaml for more information.",
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
    hazard_label,
    asset_gpkg,
    asset_layer,
    days,
    flood_thresholds,
    cyclone_damage_curve_change,
    output_dir,
):
    if not os.path.isfile(cost_file):
        raise FileNotFoundError(f"Cost file {cost_file} does not exist.")

    adapt_hazards = [
        {"hazard": "flooding", "hazard_type": ["coastal", "fluvial", "surface"]},
        {"hazard": "TC", "hazard_type": ["cyclone"]},
    ]
    rcps = [2.6, 4.5, 8.5]
    risk_type = ["EAD", "EAEL"]
    val_type = ["amin", "mean", "amax"]

    risk_filepath = os.path.join("loss_damage_npvs", f"{asset_gpkg}_{asset_layer}_EAD_EAEL_npvs.csv")

    no_adapt_risk_file = os.path.join(output_dir, risk_filepath)
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

    hazard = next((h for h in adapt_hazards if h["hazard"] == hazard_label), None)

    asset_df = pd.read_csv(network_csv)
    asset_data_details = asset_df[(asset_df["asset_gpkg"] == asset_gpkg) & (asset_df["asset_layer"] == asset_layer)]
    if len(asset_data_details) > 1:
        raise ValueError((f"Multiple assets found for gpkg={asset_gpkg} " f"and layer={asset_layer}"))
    elif len(asset_data_details) == 0:
        raise ValueError(f"No asset found for gpkg={asset_gpkg} and layer={asset_layer}")
    asset_info = asset_data_details.squeeze()

    asset_id = asset_info.asset_id_column

    logging.info(f"{hazard['hazard']} {asset_gpkg} {asset_layer}")

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
        hazard["hazard"],
        hazard["hazard_type"],
        rcps,
        risk_type,
        val_type,
        days=days,
    )
    no_adapt_ead_eael_df, ead_eael_columns = get_ead_eael(
        no_adapt_risk,
        asset_id,
        hazard["hazard"],
        hazard["hazard_type"],
        rcps,
        risk_type,
        val_type,
    )
    bcr_results = []
    ead_eael_results = []

    for option in adaptation_options:
        option_df = cost_df[cost_df["adaptation_option"] == option]
        asset_adaptation_cost = option_df["asset_adaptation_cost"].values[0]
        if hazard_label == "flooding" and asset_adaptation_cost == "J$/m":
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
                flood_thresholds,
                flood_thresholds,
                "flood_depth_protection_level",
                "flood_threshold",
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
                flood_thresholds,
                flood_thresholds,
                "flood_depth_protection_level",
                "flood_threshold",
            )
        elif hazard_label == "TC" and asset_info.sector == "energy":
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
                cyclone_damage_curve_change,
                [1],
                "cyclone_damage_curve_reduction",
                "cyclone_damage_curve_change",
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
                cyclone_damage_curve_change,
                [1],
                "cyclone_damage_curve_reduction",
                "cyclone_damage_curve_change",
            )
        else:
            option_df["flood_protection_level"] = "All"
            adapt_benefit_columns = [c.replace("risk", "avoided_risk") for c in risk_columns]
            risk_df = no_adapt_risk_df.copy()
            risk_df.rename(
                columns=dict(list(zip(risk_columns, adapt_benefit_columns))),
                inplace=True,
            )
            bcr_r, _ = bcr_estimates(
                asset_id,
                option_df,
                risk_df,
                "flood_protection_level",
                adapt_benefit_columns,
            )
            bcr_results.append(bcr_r)

            adapt_ead_eael_columns = [
                c.replace("EAD", "avoided_EAD").replace("EAEL", "avoided_EAEL") for c in ead_eael_columns
            ]
            ead_eael_df = no_adapt_ead_eael_df.copy()
            ead_eael_df.rename(
                columns=dict(list(zip(ead_eael_columns, adapt_ead_eael_columns))),
                inplace=True,
            )
            ead_eael_r = ead_eael_estimates(
                asset_id,
                option_df,
                ead_eael_df,
                "flood_protection_level",
            )
            ead_eael_results.append(ead_eael_r)

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
