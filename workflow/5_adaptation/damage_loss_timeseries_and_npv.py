import logging
import os
import warnings

import click
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm

from jamaica_infrastructure.adaptation import (
    calculate_discounting_rate_factor,
    calculate_growth_rate_factor,
    extract_growth_rate_info,
)

warnings.simplefilter(action="ignore", category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

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
):
    years = sorted(list(set(summarised_damages.epoch.values.tolist())))
    start_year = years[0]
    end_year = years[-1]

    if start_year < baseline_year:
        summarised_damages.loc[
            summarised_damages.epoch == start_year, "epoch"
        ] = baseline_year
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
    logging.debug(hazard_rcp)
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
            gr_rates = calculate_growth_rate_factor(
                gr_rates, start_year, end_year
            )
            df[series[0]] = np.multiply(df[series[0]], gr_rates)
        damages_time_series.append(df)
        df_copy = df.copy()
        df_copy[series[0]] = np.multiply(df_copy[series[0]], dsc_rate)
        df_copy[
            f"{haz}__rcp_{rcp}__{risk_type}_{val_type}"
        ] = df_copy[series[0]].sum(axis=1)
        discounted_values.append(
            df_copy[[asset_id, f"{haz}__rcp_{rcp}__{risk_type}_{val_type}"]]
        )
        del df, df_copy

    damages_time_series = pd.concat(
        damages_time_series, axis=0, ignore_index=False
    )
    index_columns = [
        c
        for c in damages_time_series.columns.values.tolist()
        if c not in timeseries
    ]

    return (
        damages_time_series[index_columns + list(timeseries)],
        discounted_values,
    )


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        readable=True
    ),
    help="Path to the asset layer data",
)
@click.option(
    "--growth-rates-xls",
    "-r",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        file_okay=True,
        readable=True
    ),
    help="Path to growth rates data",
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
    "--baseline-year",
    "-b",
    default=2019,
    required=False,
    type=int,
    help="Baseline year",
)
@click.option(
    "--projection-end-year",
    "-p",
    default=2100,
    required=False,
    type=int,
    help="Projection end year",
)
@click.option(
    "--discounting-rate",
    "-d",
    default=10,
    required=False,
    type=float,
    help="Discounting rate",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        readable=True
    ),
    help="Path to the output directory",
)
def damage_loss_timeseries_and_npv(
    network_csv,
    growth_rates_xls,
    asset_gpkg,
    asset_layer,
    baseline_year,
    projection_end_year,
    discounting_rate,
    output_path
):
    asset_df = pd.read_csv(network_csv)
    asset_data_details = asset_df[
        (asset_df["asset_gpkg"] == asset_gpkg) &
        (asset_df["asset_layer"] == asset_layer)
    ]
    asset_info = asset_data_details.squeeze()
    asset_prefix = f"{asset_gpkg}_{asset_layer}"

    growth_rates = pd.read_excel(
        growth_rates_xls,
        sheet_name="Sheet1",
    ).fillna(0)

    summarised_damages_csv = os.path.join(
        output_path,
        "direct_damages_summary",
        f"{asset_prefix}_EAD_EAEL.csv"
    )
    summarised_damages = pd.read_csv(summarised_damages_csv)

    discounted_results = os.path.join(
        output_path,
        "loss_damage_npvs"
    )
    timeseries_results = os.path.join(
        output_path,
        "loss_damage_timeseries"
    )

    if os.path.exists(discounted_results) is False:
        os.mkdir(discounted_results)
    if os.path.exists(timeseries_results) is False:
        os.mkdir(timeseries_results)

    asset_id = asset_info.asset_id_column
    index_columns = [asset_id, "damage_cost_unit", "hazard"]

    discounted_values = []
    for risk_type in ["EAD", "EAEL"]:
        for val_type in ["min", "mean", "max"]:
            if risk_type == "EAEL":
                eael_exists = [
                    c
                    for c in summarised_damages.columns.values.tolist()
                    if "EAEL_" in c
                ]
                if len(eael_exists) > 0:
                    index_columns = [
                        asset_id,
                        "economic_loss_unit",
                        "hazard",
                    ]

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
                        )
                    )
            else:
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
                    )
                )

            timeseries_csv = os.path.join(
                timeseries_results,
                f"{asset_prefix}_{risk_type}_timeseries_{val_type}.csv"
            )
            damages_time_series.to_csv(
                timeseries_csv,
                index=False,
            )
            logging.info(timeseries_csv)

    dfs = [df.set_index(asset_id) for df in discounted_values]
    discounted_values = pd.concat(dfs, axis=1).fillna(0)
    discounted_values = discounted_values.reset_index()
    discounted_values["damage_cost_unit"] = summarised_damages[
        "damage_cost_unit"
    ].values[0]
    discounted_values["economic_loss_unit"] = summarised_damages[
        "economic_loss_unit"
    ].values[0]
    discounted_values_csv = os.path.join(
        discounted_results,
        f"{asset_prefix}_EAD_EAEL_npvs.csv"
    )
    discounted_values.to_csv(
        discounted_values_csv,
        index=False,
    )

    logging.info(discounted_values_csv)


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(process)d %(filename)s %(message)s",
        level=logging.INFO,
    )
    damage_loss_timeseries_and_npv()
