"""
Calculate Benefit Cost Ratio (BCR) of adaptation options for given asset class.
"""

import logging
import os
import warnings
import math
import click
import pandas as pd
from tqdm import tqdm

warnings.simplefilter(action="ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

tqdm.pandas()


def get_benefits(id_column, df1, df2, df1_columns, df2_columns):
    modified_columns = dict([(c, c.replace("risk", "protection_risk")) for c in df2_columns])
    df2.rename(columns=modified_columns, inplace=True)
    benefit_columns = [c.replace("risk", "avoided_risk") for c in df1_columns]

    df1 = pd.merge(df1, df2, how="left", on=[id_column]).fillna(0)
    for i, (b, c) in enumerate(list(zip(benefit_columns, df1_columns))):
        if c in df2_columns:
            df1[b] = df1[c] - df1[c.replace("risk", "protection_risk")]
        else:
            df1[b] = df1[c]

    return df1, benefit_columns


def get_ead_eael_benefits(id_column, df1, df2, df1_columns, df2_columns):
    modified_columns = dict(
        [(c, c.replace("EAD", "protection_EAD").replace("EAEL", "protection_EAEL")) for c in df2_columns]
    )
    df2.rename(columns=modified_columns, inplace=True)
    benefit_columns = [c.replace("EAD", "avoided_EAD").replace("EAEL", "avoided_EAEL") for c in df1_columns]

    df1 = pd.merge(df1, df2, how="left", on=[id_column]).fillna(0)
    for i, (b, c) in enumerate(list(zip(benefit_columns, df1_columns))):
        if c in df2_columns:
            df1[b] = df1[c] - df1[c.replace("EAD", "protection_EAD").replace("EAEL", "protection_EAEL")]
        else:
            df1[b] = df1[c]

    return df1, benefit_columns


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
    all_columns, rcp_columns = get_all_column_combinations(hazard_types, rcps, risk_type, val_type)
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


def get_ead_eael(df, asset_id, hazard, hazard_types, rcps, risk_type, val_type):
    all_columns, rcp_columns = get_all_column_combinations(hazard_types, rcps, risk_type, val_type)
    all_columns = [c for c in df.columns.values.tolist() if c in all_columns]

    risk_rcp_columns = []
    for rt in risk_type:
        for ri, (rcp, rcp_c) in enumerate(list(zip(rcps, rcp_columns))):
            rcp_c = [c for c in rcp_c if c in df.columns.values.tolist()]
            if len(rcp_c) > 0:
                for vt in val_type:
                    vt_cols = [c for c in rcp_c if f"{rt}_{vt}" in c]
                    risk_rcp_columns.append(f"{hazard}__rcp_{rcp}__{rt}_{vt}")
                    column_name = f"{hazard}__rcp_{rcp}__{rt}_{vt}"
                    df[column_name] = df[vt_cols].sum(axis=1)

    return df[[asset_id] + risk_rcp_columns], risk_rcp_columns


def bcr_estimates(
    asset_id,
    option_cost_df,
    risk_df,
    hazard_thresholds_column_name,
    adapt_benefit_columns,
):
    option_cost_df["cost_units"] = "J$"
    risk_df = pd.merge(
        option_cost_df[
            [
                asset_id,
                "adaptation_option",
                hazard_thresholds_column_name,
                "cost_units",
                "adapt_cost_npv",
            ]
        ],
        risk_df,
        how="left",
        on=[asset_id],
    ).fillna(0)
    risk_df = risk_df[risk_df["adapt_cost_npv"] > 0]
    bcr_columns = [c.replace("avoided_risk", "BCR") for c in adapt_benefit_columns]
    risk_df[bcr_columns] = risk_df[adapt_benefit_columns].div(risk_df["adapt_cost_npv"], axis=0)

    return risk_df, bcr_columns


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
    hazard_thresholds,
    cost_multiplication_factors,
    hazard_thresholds_column_name,
    protection_type_name,
    protection_asset_dict,
    flood_params,
    protection_asset_breakdown,
    days=10,
):
    if isinstance(hazard_thresholds, str):
        # Added to handle coasatal adaptation option
        folder_name = protection_type_name

        risk_file = os.path.join(
            results_path,
            folder_name,
            risk_filepath,
        )

        # breakpoint()
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

    else:
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
    protection_asset_dict,
    flood_params,
    protection_asset_breakdown,
):
    if isinstance(hazard_thresholds, str):
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
    else:
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


def write_empty(asset_id: str, output_bcr: str, output_ead: str) -> None:
    """
    Write out CSV files with no data and a simple header.
    """
    for path in (output_bcr, output_ead):
        pd.DataFrame(
            [],
            columns=(asset_id, "adaptation_option", "protection_level", "adapt_cost_npv")
        ).to_csv(path, index=False)
    return


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
    protection_asset_breakdown,
    protection_asset_dict,
    hazard_label,
    asset_gpkg,
    asset_layer,
    proj_end_year,
    rcp,
    rp,
    days,
    flood_thresholds,
    cyclone_damage_curve_change,
    output_dir,
):
    if not os.path.isfile(cost_file):
        raise FileNotFoundError(f"Cost file {cost_file} does not exist.")

    adapt_hazards = [
        {"hazard": "flooding", "hazard_type": ["coastal", "fluvial", "surface"]},
        {"hazard": "coastal", "hazard_type": ["coastal"]},
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
    flood_params = [rcp, rp, proj_end_year]

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
                flood_thresholds,
                flood_thresholds,
                "flood_depth_protection_level",
                "flood_threshold",
                protection_asset_dict,
                flood_params,
                protection_asset_breakdown,
            )
        elif hazard_label == "coastal" and asset_adaptation_cost == "J$/m":
            flood_thresholds = "coastal"
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
                flood_thresholds,
                flood_thresholds,
                "flood_depth_protection_level",
                "coastal_adaptation",
                protection_asset_dict,
                flood_params,
                protection_asset_breakdown,
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
                cyclone_damage_curve_change,
                [1],
                "cyclone_damage_curve_reduction",
                "cyclone_damage_curve_change",
                protection_asset_dict,
                flood_params,
                protection_asset_breakdown,
            )
        else:
            option_df["flood_protection_level"] = "All"
            adapt_benefit_columns = [c.replace("risk", "avoided_risk") for c in risk_columns]
            risk_df = no_adapt_risk_df.copy()
            risk_df.rename(
                columns=dict(list(zip(risk_columns, adapt_benefit_columns))),
                inplace=True,
            )
            bcr_r, bcr_cols = bcr_estimates(
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
