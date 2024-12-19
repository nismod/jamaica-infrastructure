"""Estimate damage and loss timeseries and NPVs

"""

import datetime
import os
import sys
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import geopandas as gpd
import numpy as np
from tqdm import tqdm

from analysis_utils import *

tqdm.pandas()


def estimate_time_series(
    summarised_damages,
    asset_id,
    index_columns,
    risk_type,
    val_type,
    baseline_year,
    projection_end_year,
    growth_rates,
    discounting_rate,
    discounted_values,
    timeseries_results,
):
    years = sorted(list(set(summarised_damages.epoch.values.tolist())))
    start_year = years[0]
    end_year = years[-1]

    if start_year < baseline_year:
        summarised_damages.loc[summarised_damages.epoch == start_year, "epoch"] = (
            baseline_year
        )
        start_year = baseline_year
    if end_year < projection_end_year:
        end_year = projection_end_year

    dsc_rate = calculate_discounting_rate_factor(
        discount_rate=discounting_rate,
        start_year=start_year,
        end_year=end_year,
        maintain_period=1,
    )
    timeseries = np.arange(start_year, end_year + 1, 1)
    hazard_rcp = list(
        set(
            zip(
                summarised_damages.hazard.values.tolist(),
                summarised_damages.rcp.values.tolist(),
            )
        )
    )
    print(hazard_rcp)
    hazard_rcp = [hz_rcp for hz_rcp in hazard_rcp if hz_rcp[1] != "baseline"]

    defence_column = [
        c
        for c in summarised_damages.columns.values.tolist()
        if f"{risk_type}_" in c and f"_{val_type}" in c
    ][0]
    damages_time_series = []
    for ix, (haz, rcp) in enumerate(hazard_rcp):
        haz_rcp_damages = summarised_damages[
            (summarised_damages["hazard"] == haz)
            & (summarised_damages["rcp"].isin(["baseline", rcp]))
        ]
        years = sorted(list(set(haz_rcp_damages.epoch.values.tolist())))
        df = (
            haz_rcp_damages.set_index(index_columns)
            .pivot(columns="epoch")[defence_column]
            .reset_index()
            .rename_axis(None, axis=1)
        ).fillna(0)
        series = np.array([list(timeseries) * len(df.index)]).reshape(
            len(df.index), len(timeseries)
        )
        df["rcp"] = rcp
        if len(years) > 1:
            df[series[0]] = interp1d(
                years, df[years], fill_value="extrapolate", bounds_error=False
            )(series[0])
        else:
            t = [s for s in series[0] if s not in years]
            df[t] = 0
        df[series[0]] = df[series[0]].clip(lower=0.0)
        if risk_type == "EAEL":
            gr_rates = extract_growth_rate_info(
                growth_rates,
                "year",
                f"gdp_{val_type}",
                start_year=start_year,
                end_year=end_year,
            )
            gr_rates = calculate_growth_rate_factor(gr_rates, start_year, end_year)
            df[series[0]] = np.multiply(df[series[0]], gr_rates)
        damages_time_series.append(df)
        df_copy = df.copy()
        df_copy[series[0]] = np.multiply(df_copy[series[0]], dsc_rate)
        df_copy[f"{haz}__rcp_{rcp}__{risk_type}_{val_type}"] = df_copy[series[0]].sum(
            axis=1
        )
        discounted_values.append(
            df_copy[[asset_id, f"{haz}__rcp_{rcp}__{risk_type}_{val_type}"]]
        )
        del df, df_copy

    damages_time_series = pd.concat(damages_time_series, axis=0, ignore_index=False)
    index_columns = [
        c for c in damages_time_series.columns.values.tolist() if c not in timeseries
    ]

    return damages_time_series[index_columns + list(timeseries)], discounted_values


def main(
    config,
    summary_results_folder,
    timeseries_results_folder,
    discounted_results_folder,
    network_csv,
    baseline_year=2019,
    projection_end_year=2100,
    discounting_rate=10,
):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]

    summary_results = os.path.join(output_data_path, summary_results_folder)

    timeseries_results = os.path.join(output_data_path, timeseries_results_folder)
    if os.path.exists(timeseries_results) == False:
        os.mkdir(timeseries_results)

    discounted_results = os.path.join(output_data_path, discounted_results_folder)
    if os.path.exists(discounted_results) == False:
        os.mkdir(discounted_results)

    asset_data_details = pd.read_csv(network_csv)

    growth_rates = pd.read_excel(
        os.path.join(
            processed_data_path, "macroeconomic_data", "gdp_growth_rates.xlsx"
        ),
        sheet_name="Sheet1",
    ).fillna(0)

    for asset_info in asset_data_details.itertuples():
        asset_id = asset_info.asset_id_column
        index_columns = [asset_id, "damage_cost_unit", "hazard"]
        file = os.path.join(
            summary_results,
            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL.csv",
        )
        if os.path.isfile(file) is True:
            summarised_damages = pd.read_csv(file)

            discounted_values = []
            for risk_type in ["EAD", "EAEL"]:
                for val_type in ["amin", "mean", "amax"]:
                    if risk_type == "EAEL":
                        eael_exists = [
                            c
                            for c in summarised_damages.columns.values.tolist()
                            if "EAEL_" in c
                        ]
                        if len(eael_exists) > 0:
                            index_columns = [asset_id, "economic_loss_unit", "hazard"]

                            damages_time_series, discounted_values = (
                                estimate_time_series(
                                    summarised_damages,
                                    asset_id,
                                    index_columns,
                                    risk_type,
                                    val_type,
                                    baseline_year,
                                    projection_end_year,
                                    growth_rates,
                                    discounting_rate,
                                    discounted_values,
                                    timeseries_results,
                                )
                            )
                    else:
                        damages_time_series, discounted_values = estimate_time_series(
                            summarised_damages,
                            asset_id,
                            index_columns,
                            risk_type,
                            val_type,
                            baseline_year,
                            projection_end_year,
                            growth_rates,
                            discounting_rate,
                            discounted_values,
                            timeseries_results,
                        )

                    damages_time_series.to_csv(
                        os.path.join(
                            timeseries_results,
                            f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_{risk_type}_timeseries_{val_type}.csv",
                        ),
                        index=False,
                    )
                    print(
                        f"* Done with {asset_info.asset_gpkg} {asset_info.asset_layer} {risk_type} timeseries {val_type}"
                    )

            dfs = [df.set_index(asset_id) for df in discounted_values]
            discounted_values = pd.concat(dfs, axis=1).fillna(0)
            discounted_values = discounted_values.reset_index()
            discounted_values["damage_cost_unit"] = summarised_damages[
                "damage_cost_unit"
            ].values[0]
            discounted_values["economic_loss_unit"] = summarised_damages[
                "economic_loss_unit"
            ].values[0]
            discounted_values.to_csv(
                os.path.join(
                    discounted_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_EAD_EAEL_npvs.csv",
                ),
                index=False,
            )

            print(f"* Done with {asset_info.asset_gpkg} discounted values")

    # signal to snakemake that the job is complete
    with open(os.path.join(summary_results, "losses_npv.flag"), "w") as fp:
        fp.write(f"{datetime.datetime.now()}")


if __name__ == "__main__":
    CONFIG = load_config()
    # try:
    #     summary_results_folder = str(sys.argv[1])
    #     timeseries_results_folder = str(sys.argv[2])
    #     discounted_results_folder = str(sys.argv[3])
    #     network_csv = str(sys.argv[4])
    # except IndexError:
    #     print("Got arguments", sys.argv)
    #     exit()

    summary_results_folder = str(sys.argv[1])
    timeseries_results_folder = str(sys.argv[2])
    discounted_results_folder = str(sys.argv[3])
    network_csv = str(sys.argv[4])
    main(
        CONFIG,
        summary_results_folder,
        timeseries_results_folder,
        discounted_results_folder,
        network_csv,
    )
