import math

import numpy as np
import pandas as pd


def prod(val):
    res = 1
    for ele in val:
        res *= ele
    return res


def calculate_growth_rate_factor(
    growth_rate=5.0,
    start_year=2020,
    end_year=2050,
):
    growth_rates = []

    if isinstance(growth_rate, float):
        for year in range(start_year, end_year + 1):
            growth_rates.append(
                1.0 * math.pow(1.0 + 1.0 * growth_rate / 100.0, year - start_year)
            )
    else:
        for i, (year, rate) in enumerate(growth_rate):
            if year > start_year:
                growth_rates.append(prod([1 + v[1] / 100.0 for v in growth_rate[:i]]))
            else:
                growth_rates.append(1)

    return np.array(growth_rates)


def calculate_discounting_rate_factor(
    discount_rate=4.5,
    start_year=2020,
    end_year=2050,
    maintain_period=4,
    skip_year_one=False,
):
    discount_rates = []
    maintain_years = np.arange(start_year + 1, end_year + 1, maintain_period)
    for year in range(start_year, end_year + 1):
        if year in maintain_years:
            discount_rates.append(
                1.0 / math.pow(1.0 + 1.0 * discount_rate / 100.0, year - start_year)
            )
        else:
            if skip_year_one is True:
                discount_rates.append(0)
            else:
                discount_rates.append(1)

    return np.array(discount_rates)


def extract_growth_rate_info(
    growth_rates, time_column, rate_column, start_year=2000, end_year=2100
):
    growth_year_rates = []
    growth_rates_times = list(sorted(growth_rates[time_column].values.tolist()))
    # And create parameter values
    for y in range(start_year, end_year + 1):
        if y in growth_rates_times:
            growth_year_rates.append(
                (
                    y,
                    growth_rates.loc[
                        growth_rates[time_column] == y, rate_column
                    ].values[0],
                )
            )
        elif y < growth_rates_times[0]:
            growth_year_rates.append(
                (
                    y,
                    growth_rates.loc[
                        growth_rates[time_column] == growth_rates_times[0], rate_column
                    ].values[0],
                )
            )
        elif y > growth_rates_times[-1]:
            growth_year_rates.append(
                (
                    y,
                    growth_rates.loc[
                        growth_rates[time_column] == growth_rates_times[-1], rate_column
                    ].values[0],
                )
            )

    return growth_year_rates


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