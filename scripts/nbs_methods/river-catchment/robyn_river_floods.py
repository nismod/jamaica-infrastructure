import re

import numpy as np
import pandas
import scipy
from scipy.integrate import trapezoid


# Keys in these dictionaries must match the "Classify" values in landuse shapefile
# N.B. deliberate matching of typos and double spaces
LANDUSE_FOREST_PROPORTION = {
    "Bamboo": 1,
    "Bamboo and Fields": 0.5,  # 50% Ag, 50% Bamboo (if you consider bamboo here)
    "Bamboo and Secondary Forest": 1,  # 50% Forest, 50% Ag (if you want some ag fraction)
    "Bare Rock": 0,
    "Bauxite Extraction": 0,
    "Buildings and other infrastructures": 0,
    "Closed broadleaved forest (Primary Forest)": 1,
    "Disturbed broadleaved forest (Secondary Forest)": 1,
    "Fields  and Bamboo": 0.5,  # 50% Ag, 50% Bamboo
    "Fields and Secondary Forest": 0.5,  # 50% Ag, 50% Forest
    "Fields or Secondary Forest/Pine Plantation": 0.5,  # 50% Ag, 50% Secondary Forest/Pine
    "Fields: Bare Land": 0,
    "Fields: Herbaceous crops, fallow, cultivated vegetables": 0,
    "Fields: Pasture,Human disturbed, grassland": 0,
    "Hardwood Plantation: Euculytus": 1,
    "Hardwood Plantation: Mahoe": 1,
    "Hardwood Plantation: Mahogany": 1,
    "Hardwood Plantation: Mixed": 1,
    "Herbaceous Wetland": 0,
    "Mangrove Forest": 0,
    "Open dry forest - Short": 1,
    "Open dry forest - Tall (Woodland/Savanna)": 1,
    "Plantation: Tree crops, shrub crops, sugar cane, banana": 1,
    "Quarry": 0,
    "Secondary Forest": 1,
    "Swamp Forest": 0,
    "Water Body": 0,
}


LANDUSE_AFFORESTABLE_PROPORTION = {
    "Bamboo": 0,
    "Bamboo and Fields": 0.5,  # 50% Ag, 50% Bamboo (if you consider bamboo here)
    "Bamboo and Secondary Forest": 0,  # 50% Forest, 50% Ag (if you want some ag fraction)
    "Bare Rock": 0,
    "Bauxite Extraction": 1,
    "Buildings and other infrastructures": 0,
    "Closed broadleaved forest (Primary Forest)": 0,
    "Disturbed broadleaved forest (Secondary Forest)": 0,
    "Fields  and Bamboo": 0.5,  # 50% Ag, 50% Bamboo
    "Fields and Secondary Forest": 0.5,  # 50% Ag, 50% Forest
    "Fields or Secondary Forest/Pine Plantation": 0.5,  # 50% Ag, 50% Secondary Forest/Pine
    "Fields: Bare Land": 1,
    "Fields: Herbaceous crops, fallow, cultivated vegetables": 1,
    "Fields: Pasture,Human disturbed, grassland": 1,
    "Hardwood Plantation: Euculytus": 0,
    "Hardwood Plantation: Mahoe": 0,
    "Hardwood Plantation: Mahogany": 0,
    "Hardwood Plantation: Mixed": 0,
    "Herbaceous Wetland": 0,
    "Mangrove Forest": 0,
    "Open dry forest - Short": 0,
    "Open dry forest - Tall (Woodland/Savanna)": 0,
    "Plantation: Tree crops, shrub crops, sugar cane, banana": 0,
    "Quarry": 1,
    "Secondary Forest": 0,
    "Swamp Forest": 0,
    "Water Body": 0,
}


def interpolate_damages(
    rp: float, rp_l: float, value_l: np.ndarray, rp_u: float, value_u: np.ndarray
) -> np.ndarray:
    """Interpolate between two return period ndarrays"""
    rp_factor = (np.log(rp) - np.log(rp_l)) / (np.log(rp_u) - np.log(rp_l))
    value = value_l + ((value_u - value_l) * rp_factor)

    return value


def pick_upper_lower_rps(rp: float, rps: list[float]) -> tuple[float]:
    bin_index = np.searchsorted(rps, rp, side="left")
    rp_l = rps[bin_index - 1]
    rp_u = rps[bin_index]
    return rp_l, rp_u


def get_rp_cols(df):
    rp_cols = [col for col in df.columns if "rp" in col]
    rps = [float(col.replace("rp", "")) for col in rp_cols]
    return rp_cols, rps


def interpolate_rp_damages(rps_to_calculate, damages):
    rp_cols, rps = get_rp_cols(damages)
    interpolated_damages = pandas.DataFrame()
    rps = sorted(rps)

    # print(f"{rps:} {rps_to_calculate:}")

    # Calculate and save interpolated damages for new return periods
    for rp in rps_to_calculate:
        rp_l, rp_u = pick_upper_lower_rps(rp, rps)
        rp_damages = interpolate_damages(
            rp, rp_l, damages[f"rp{rp_l}"], rp_u, damages[f"rp{rp_u}"]
        )
        interpolated_damages[f"rp{rp}"] = rp_damages

    return interpolated_damages


def peak_flow_reduction(forest_percentage_change):
    # set up data frame with values from literature
    lookup = pandas.DataFrame(
        {
            "catchment_forest_percentage": [0, 6, 14, 21, 62, 100],
            "rp5.0": [0, 3, 13, 18, 48, 48],
            "rp100.0": [0, 1, 8, 11, 32, 32],
        }
    )
    data = lookup.catchment_forest_percentage
    rp_cols, rps = get_rp_cols(lookup)

    # set up interpolator
    interpolator = scipy.interpolate.RegularGridInterpolator(
        (rps, data),
        lookup[rp_cols].values.T,
        method="linear",
    )

    # set up return periods to get peak flow reduction values (includes 5 and 100)
    interp_rps = [5.0, 10.0, 20.0, 50.0, 100.0]

    # do the interpolation
    vals = interpolator(
        (
            [interp_rps],
            np.array(forest_percentage_change).reshape(
                len(forest_percentage_change), 1
            ),
        )
    ).T
    data = pandas.DataFrame(
        {
            "forest_percentage_change": forest_percentage_change,
        }
    )
    for rp, val in zip(interp_rps, vals):
        data[f"rp{rp}"] = val

    return data


def rp_change_given_flow_reduction(reduction_percent, interp_rp):
    ccra_flow_reductions = pandas.DataFrame(
        {
            "reduction_percent": [0.0, 5.0, 10.0, 20.0, 40.0, 100.0],
            "rp2.0": [2.0, 2.4, 3.2, 6.8, 169, 169],
            "rp2.3": [2.3, 2.8, 3.7, 7.9, 196, 196],
            "rp5.0": [5.0, 6.5, 8.9, 19, 235, 235],
            "rp10.0": [10.0, 13, 19, 43, 473, 473],
            "rp25.0": [25.0, 35, 51, 123, 1330, 1330],
            "rp50.0": [50.0, 72, 107, 268, 2967, 2967],
            "rp100.0": [100.0, 147, 224, 582, 6648, 6648],
            "rp500.0": [500.0, 772, 1229, 3441, 43094, 43094],
            "rp1000.0": [1000.0, 1571, 2544, 7345, 95943, 95943],
        }
    )

    # proportion of baseline flow
    ccra_flow_reductions["flow"] = 1 - ccra_flow_reductions.reduction_percent / 100

    rp_cols, rps = get_rp_cols(ccra_flow_reductions)

    data = ccra_flow_reductions.flow
    flow_to_interpolate = 1 - reduction_percent / 100

    interpolator = scipy.interpolate.RegularGridInterpolator(
        (rps, data),
        ccra_flow_reductions[rp_cols].values.T,
        method="linear",
    )

    # set up return periods to get peak flow reduction values (includes 5 and 100)
    interp_rps = [interp_rp]

    # do the interpolation
    vals = interpolator(
        (
            [interp_rps],
            np.array(flow_to_interpolate).reshape(len(flow_to_interpolate), 1),
        )
    ).T
    data = pandas.DataFrame(
        {
            "reduction_percent": reduction_percent,
        }
    )
    for rp, val in zip(interp_rps, vals):
        data[f"rp{rp}"] = val
    return data


def select_damages(sector_damages, variant="mean"):
    rps = [20.0, 50.0, 100.0, 200.0, 500.0, 1500.0]
    colnames = [
        f"fluvial__rp_{int(rp)}__rcp_baseline__epoch_2010__conf_None_{variant}"
        for rp in rps
    ]
    to_rename = {colname: f"rp{rp}" for colname, rp in zip(colnames, rps)}
    selected_damages = (
        sector_damages.set_index("HYBAS_ID")[colnames].rename(columns=to_rename).copy()
    )
    selected_damages["rp1000000000.0"] = selected_damages["rp1500.0"].copy()
    selected_damages["rp0.0001"] = 0
    return selected_damages


def process_sector_damage_for_rp(rp, all_sector_damage, hydrobasins):
    """
    Process the damages for a given return period.

    Parameters:
      rp (int or float): the return period (e.g. 20, 50, 100)
      all_sector_damage (pd.DataFrame): DataFrame containing damage data.
      hydrobasins (pd.DataFrame): DataFrame with geometry and HYBAS_ID.

    Returns:
      pd.DataFrame: A DataFrame with standardized columns:
                  ["HYBAS_ID", "geometry", "damages", "damages_with_nbs",
                   "avoided_damages", "rp"]
    """
    # Define column names based on the specified return period
    fluvial_col = f"baseline__fluvial__rp_{rp}__baseline__fluvial__rp_mean"
    future_col = f"future__fluvial__rp_{rp}__future__fluvial__rp_mean"

    # Group the damage data by HYBAS_ID and sum (in case you have duplicate HYBAS_IDs)
    damage_df = (
        all_sector_damage[["HYBAS_ID", fluvial_col, future_col]]
        .groupby("HYBAS_ID")
        .sum()
    )

    # Join with the hydrobasins geometry
    damage_df = (
        hydrobasins[["HYBAS_ID", "geometry"]].set_index("HYBAS_ID").join(damage_df)
    )

    # Calculate avoided damages (baseline minus future)
    damage_df[f"avoided__fluvial__rp_{rp}"] = (
        damage_df[fluvial_col] - damage_df[future_col]
    )

    # Add the return period column
    damage_df["rp"] = rp

    # Reset the index and rename the columns consistently:
    #   "damages" will be taken as the baseline fluvial damage column and
    #   "damages_with_nbs" as the future fluvial damage column,
    #   "avoided_damages" from the calculated column.
    damage_df = damage_df.reset_index()
    damage_df.columns = [
        "HYBAS_ID",
        "geometry",
        "damages",
        "damages_with_nbs",
        "avoided_damages",
        "rp",
    ]

    return damage_df


def calculate_ead(df):
    rp_cols = [col for col in df.columns if re.match(r"^rp\d+(\.\d+)?$", col)]
    rp_vals = sorted([float(col.replace("rp", "")) for col in rp_cols])
    rp_vals = np.array(rp_vals[::-1])
    probabilities = 1 / rp_vals
    rp_damages = df[[f"rp{r}" for r in rp_vals]]
    # print("RP", rp_damages.columns)
    return trapezoid(rp_damages, x=probabilities, axis=1)
