"""Estimate direct damages to physical assets exposed to hazards

"""

import sys
import os
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import geopandas as gpd
from shapely import wkb
import numpy as np
from SALib.sample import morris
import SALib.analyze.morris

from analysis_utils import *
from tqdm import tqdm

tqdm.pandas()


def select_damage_curves(
    asset_df, asset_type_column, damage_curve_lookup_df, assigned_sector, hazard_type
):
    selected_assets = damage_curve_lookup_df[
        (damage_curve_lookup_df["sector"] == assigned_sector)
        & (damage_curve_lookup_df["hazard_type"] == hazard_type)
    ][["asset_name", "asset_sheet"]]
    data = pd.read_excel(
        os.path.join(damage_data_path, f"damage_curves_{sector}_{hazard_type}.xlsx"),
        sheet_name=data_key.asset_sheet,
    )
    asset_df = pd.merge(
        asset_df,
        selected_assets,
        how="left",
        left_on=asset_type_column,
        right_on="asset_name",
    )
    asset_df = asset_df[~asset_df["asset_sheet"].isna()]


def get_damage_data(x, damage_data_path, uplift_factor=0, uncertainty_parameter=0):
    data = pd.read_excel(
        os.path.join(
            damage_data_path, f"damage_curves_{x.sector}_{x.hazard_type}.xlsx"
        ),
        sheet_name=x.asset_sheet,
    )
    if x.hazard_type == "flooding":
        x_data = data.flood_depth
    else:
        x_data = data.wind_speed

    y_data = data.damage_ratio_min + uncertainty_parameter * (
        data.damage_ratio_max - data.damage_ratio_min
    )
    # print(y_data)
    y_data = np.minimum(y_data + uplift_factor, 1.0)

    return x_data.values, y_data.values


def convert_cost_units(x, cost_value, cost_unit, conversion_rate):
    if ("US$" in x[cost_unit]) or ("USD" in x[cost_unit]):
        return conversion_rate * x[cost_value]
    else:
        return x[cost_value]


def modify_cost_units(x, cost_dimension, damage_cost_column="damage_cost"):
    if "/km" in str(x[cost_dimension]):
        return 0.001 * x[damage_cost_column]
    else:
        return x[damage_cost_column]


def add_exposure_dimensions(dataframe, dataframe_type="nodes", epsg=4326):
    geo_dataframe = gpd.GeoDataFrame(
        dataframe, geometry="geometry", crs={"init": f"epsg:{epsg}"}
    )
    if dataframe_type == "edges":
        geo_dataframe["exposure"] = geo_dataframe.geometry.length
        geo_dataframe["exposure_unit"] = "m"
    elif dataframe_type == "areas":
        geo_dataframe["exposure"] = geo_dataframe.geometry.area
        geo_dataframe["exposure_unit"] = "m2"
    else:
        geo_dataframe["exposure"] = 1
        geo_dataframe["exposure_unit"] = "unit"
    geo_dataframe.drop("geometry", axis=1, inplace=True)

    index_columns = [
        c for c in geo_dataframe.columns.values.tolist() if c != "exposure"
    ]
    return (
        geo_dataframe.groupby(index_columns, dropna=False)["exposure"]
        .sum()
        .reset_index()
    )


def create_damage_curves(
    damage_data_path, damage_curve_lookup_df, uplift_factor=0, uncertainty_parameter=0
):
    damage_curve_lookup_df["x_y_data"] = damage_curve_lookup_df.progress_apply(
        lambda x: get_damage_data(
            x, damage_data_path, uplift_factor, uncertainty_parameter
        ),
        axis=1,
    )
    damage_curve_lookup_df[["damage_x_data", "damage_y_data"]] = damage_curve_lookup_df[
        "x_y_data"
    ].apply(pd.Series)
    damage_curve_lookup_df.drop("x_y_data", axis=1, inplace=True)

    return damage_curve_lookup_df


def estimate_direct_damage_costs_and_units(
    dataframe,
    damage_ratio_columns,
    cost_unit_column,
    damage_cost_column="damage_cost",
    dataframe_type="nodes",
):
    if dataframe_type == "nodes":
        dataframe[damage_ratio_columns] = dataframe[damage_ratio_columns].multiply(
            dataframe[damage_cost_column], axis="index"
        )
        dataframe["damage_cost_unit"] = dataframe[cost_unit_column]
    else:
        dataframe[damage_ratio_columns] = dataframe[damage_ratio_columns].multiply(
            dataframe[damage_cost_column] * dataframe["exposure"], axis="index"
        )
        cost_unit = dataframe[cost_unit_column].values.tolist()[0]
        dataframe["damage_cost_unit"] = "/".join(cost_unit.split("/")[:-1])

    return dataframe


def main(config, set_count, cost_uncertainty_parameter, damage_uncertainty_parameter):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    output_data_path = config["paths"]["output"]
    epsg_jamaica = 3448

    jd_to_usd = 0.0067  # Jamaican dollar to USD conversion
    direct_damages_results = os.path.join(
        output_data_path, "direct_damages_no_threshold"
    )
    if os.path.exists(direct_damages_results) == False:
        os.mkdir(direct_damages_results)

    hazard_asset_intersection_path = os.path.join(
        output_data_path, "hazard_asset_intersection"
    )
    asset_data_details = pd.read_csv(
        os.path.join(
            processed_data_path,
            "networks",
            "network_layers_hazard_intersections_details.csv",
        )
    )

    hazard_data_path = os.path.join(
        processed_data_path, "hazards", "hazard_descriptions"
    )
    damage_data_path = os.path.join(processed_data_path, "damage_curves")
    damage_curve_lookup = pd.read_csv(
        os.path.join(damage_data_path, "asset_damage_curve_mapping.csv")
    )[["sector", "hazard_type", "asset_name", "asset_sheet"]]
    # hazard_data_files = []
    # for root, dirs, files in os.walk(hazard_data_path):
    #     for file in files:
    #         if file.endswith(".csv"):
    #             hazard_data_files.append(file)
    # # print (hazard_data_files)
    hazard_data_files = ["hazard_layers.csv"]

    """Step 1: Get all the damage curves into a dataframe
    """
    hazard_attributes = [
        {
            "hazard": "coastal",
            "hazard_type": "flooding",
            "hazard_threshold": 0.0,
            "uplift_factor": 0.12,
        },
        {
            "hazard": "fluvial",
            "hazard_type": "flooding",
            "hazard_threshold": 0.0,
            "uplift_factor": 0.0,
        },
        {
            "hazard": "surface",
            "hazard_type": "flooding",
            "hazard_threshold": 0.0,
            "uplift_factor": 0.0,
        },
        {
            "hazard": "cyclone",
            "hazard_type": "TC",
            "hazard_threshold": 0.0,
            "uplift_factor": 0.0,
        },
    ]
    chosen_hazards = [h["hazard"] for h in hazard_attributes]
    flood_hazards = ["coastal", "fluvial", "surface"]
    hazard_attributes = pd.DataFrame(hazard_attributes)
    damage_curves = []
    for idx, hazard in hazard_attributes.iterrows():
        damage_curve_df = damage_curve_lookup[
            damage_curve_lookup["hazard_type"] == hazard["hazard_type"]
        ]
        damage_curve_df["hazard"] = hazard["hazard"]

        damage_curve_df = create_damage_curves(
            damage_data_path,
            damage_curve_df,
            uplift_factor=hazard["uplift_factor"],
            uncertainty_parameter=damage_uncertainty_parameter,
        )
        damage_curves.append(damage_curve_df)

    damage_curves = pd.concat(damage_curves, axis=0, ignore_index=True)
    del damage_curve_df

    for asset_info in asset_data_details.itertuples():
        asset_sector = asset_info.sector
        asset_id = asset_info.asset_id_column
        asset_min_cost = asset_info.asset_min_cost_column
        asset_max_cost = asset_info.asset_max_cost_column
        asset_cost_unit = asset_info.asset_cost_unit_column

        asset_df = gpd.read_file(
            os.path.join(processed_data_path, asset_info.path),
            layer=asset_info.asset_layer,
        )
        asset_df[asset_min_cost] = asset_df.apply(
            lambda x: convert_cost_units(
                x, asset_min_cost, asset_cost_unit, 1.0 / jd_to_usd
            ),
            axis=1,
        )
        asset_df[asset_max_cost] = asset_df.apply(
            lambda x: convert_cost_units(
                x, asset_max_cost, asset_cost_unit, 1.0 / jd_to_usd
            ),
            axis=1,
        )
        asset_df[asset_cost_unit] = asset_df[asset_cost_unit].replace(
            ["USD", "US$"], "J$", regex=True
        )
        # We just need to modify and correct energy asset costs of nodes to J$ and costs of of edges to $J/m
        if asset_sector == "energy" and asset_info.asset_layer == "edges":
            asset_df[asset_min_cost] = asset_df[asset_min_cost] / asset_df["length"]
            asset_df[asset_max_cost] = asset_df[asset_max_cost] / asset_df["length"]
            asset_df[asset_cost_unit] = "J$/m"

        asset_df["damage_cost"] = asset_df[
            asset_min_cost
        ] + cost_uncertainty_parameter * (
            asset_df[asset_max_cost] - asset_df[asset_min_cost]
        )
        asset_df["damage_cost"] = asset_df.progress_apply(
            lambda x: modify_cost_units(x, asset_cost_unit), axis=1
        )
        hazard_damages = []
        for hazard_file in hazard_data_files:
            hazard_intersection_file = os.path.join(
                hazard_asset_intersection_path,
                f"{asset_info.asset_gpkg}_splits__hazard_layers__{asset_info.asset_layer}.geoparquet",
            )
            hazard_data_details = pd.read_csv(
                os.path.join(hazard_data_path, hazard_file), encoding="latin1"
            )
            if os.path.isfile(hazard_intersection_file) is True:
                hazard_df = gpd.read_parquet(hazard_intersection_file)
                hazard_df = hazard_df.to_crs(epsg=epsg_jamaica)
                hazard_df = add_exposure_dimensions(
                    hazard_df, dataframe_type=asset_info.asset_layer, epsg=epsg_jamaica
                )
                for hazard_info in hazard_attributes.itertuples():
                    if (
                        getattr(
                            asset_info,
                            f"{hazard_info.hazard}_asset_damage_lookup_column",
                        )
                        != "none"
                    ):
                        asset_hazard = getattr(
                            asset_info,
                            f"{hazard_info.hazard}_asset_damage_lookup_column",
                        )
                        hazard_keys = hazard_data_details[
                            hazard_data_details["hazard"] == hazard_info.hazard
                        ]["key"].values.tolist()
                        hazard_effect_df = hazard_df[
                            [asset_id, "exposure", "exposure_unit"] + hazard_keys
                        ]
                        damages_df = damage_curves[
                            (damage_curves["sector"] == asset_sector)
                            & (damage_curves["hazard"] == hazard_info.hazard)
                        ]
                        damaged_assets = list(
                            set(damages_df["asset_name"].values.tolist())
                        )
                        affected_assets_df = asset_df[
                            asset_df[asset_hazard].isin(damaged_assets)
                        ][[asset_id, asset_hazard, asset_cost_unit, "damage_cost"]]
                        damaged_assets = list(
                            set(affected_assets_df[asset_hazard].values.tolist())
                        )
                        damages_df = damages_df[
                            damages_df["asset_name"].isin(damaged_assets)
                        ]
                        affected_assets = list(
                            set(affected_assets_df[asset_id].values.tolist())
                        )
                        hazard_effect_df["hazard_threshold"] = (
                            hazard_info.hazard_threshold
                        )
                        if hazard_info.hazard in flood_hazards:
                            hazard_effect_df[hazard_keys] = (
                                hazard_effect_df[hazard_keys]
                                - hazard_info.hazard_threshold
                            )
                            hazard_effect_df = hazard_effect_df[
                                (hazard_effect_df[hazard_keys] > 0).any(axis=1)
                            ]
                            hazard_effect_df = hazard_effect_df[
                                hazard_effect_df[asset_id].isin(affected_assets)
                            ]
                        else:
                            hazard_effect_df[hazard_keys] = np.where(
                                hazard_effect_df[hazard_keys]
                                <= hazard_info.hazard_threshold,
                                0,
                                hazard_effect_df[hazard_keys],
                            )
                            hazard_effect_df = hazard_effect_df[
                                (
                                    hazard_effect_df[hazard_keys]
                                    > hazard_info.hazard_threshold
                                ).any(axis=1)
                            ]
                            hazard_effect_df = hazard_effect_df[
                                hazard_effect_df[asset_id].isin(affected_assets)
                            ]

                        if len(hazard_effect_df.index) == 0:
                            print(
                                f"* No {hazard_info.hazard} intersections with {asset_info.asset_gpkg} {asset_info.asset_layer}"
                            )
                        else:
                            hazard_effect_df = pd.merge(
                                hazard_effect_df,
                                affected_assets_df,
                                how="left",
                                on=[asset_id],
                            )
                            # print (hazard_info.key)
                            for damage_info in damages_df.itertuples():
                                hazard_asset_effect_df = hazard_effect_df[
                                    hazard_effect_df[asset_hazard]
                                    == damage_info.asset_name
                                ]
                                if len(hazard_asset_effect_df.index) > 0:
                                    hazard_asset_effect_df[hazard_keys] = interp1d(
                                        damage_info.damage_x_data,
                                        damage_info.damage_y_data,
                                        fill_value=(
                                            min(damage_info.damage_y_data),
                                            max(damage_info.damage_y_data),
                                        ),
                                        bounds_error=False,
                                    )(hazard_asset_effect_df[hazard_keys])
                                    hazard_asset_effect_df = (
                                        estimate_direct_damage_costs_and_units(
                                            hazard_asset_effect_df,
                                            hazard_keys,
                                            asset_cost_unit,
                                            dataframe_type=asset_info.asset_layer,
                                        )
                                    )

                                    sum_dict = dict([(hk, "sum") for hk in hazard_keys])
                                    hazard_asset_effect_df = (
                                        hazard_asset_effect_df.groupby(
                                            [
                                                asset_id,
                                                "exposure_unit",
                                                "damage_cost_unit",
                                                "exposure",
                                            ],
                                            dropna=False,
                                        )
                                        .agg(sum_dict)
                                        .reset_index()
                                    )

                                    hazard_asset_effect_df[
                                        "damage_uncertainty_parameter"
                                    ] = damage_uncertainty_parameter
                                    hazard_asset_effect_df[
                                        "cost_uncertainty_parameter"
                                    ] = cost_uncertainty_parameter
                                    hazard_damages.append(hazard_asset_effect_df)

                                del hazard_asset_effect_df
                            del hazard_effect_df
                    else:
                        print(
                            f"* {asset_info.asset_gpkg} {asset_info.asset_layer} not affected by {hazard_info.hazard}"
                        )
        if len(hazard_damages) > 0:
            asset_damages_results = os.path.join(
                direct_damages_results,
                f"{asset_info.asset_gpkg}_{asset_info.asset_layer}",
            )
            if os.path.exists(asset_damages_results) == False:
                os.mkdir(asset_damages_results)
            hazard_damages = pd.concat(
                hazard_damages, axis=0, ignore_index=True
            ).fillna(0)
            hazard_damages.to_csv(
                os.path.join(
                    asset_damages_results,
                    f"{asset_info.asset_gpkg}_{asset_info.asset_layer}_direct_damages_parameter_set_{set_count}.csv",
                ),
                index=False,
            )


if __name__ == "__main__":
    CONFIG = load_config()
    try:
        set_count = str(sys.argv[1])
        cost_uncertainty_parameter = float(sys.argv[2])
        damage_uncertainty_parameter = float(sys.argv[3])
    except IndexError:
        print("Got arguments", sys.argv)
        exit()

    main(CONFIG, set_count, cost_uncertainty_parameter, damage_uncertainty_parameter)
