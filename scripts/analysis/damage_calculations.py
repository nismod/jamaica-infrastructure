"""Estimate direct damages to physical assets exposed to hazards

"""

import logging
import os
import warnings

import click
import pandas as pd
import geopandas as gpd
import numpy as np
from tqdm import tqdm

from jamaica_infrastructure.direct_damages import interp1d

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
pd.options.mode.chained_assignment = None  # default='warn'

tqdm.pandas()

epsg_jamaica = 3448
jd_to_usd = 0.0067  # Jamaican dollar to USD conversion


def get_damage_data(x, damage_data_path, uplift_factor=0, uncertainty_parameter=0):
    data = pd.read_excel(
        os.path.join(damage_data_path, f"damage_curves_{x.sector}_{x.hazard_type}.xlsx"),
        sheet_name=x.asset_sheet,
    )
    if x.hazard_type == "flooding":
        x_data = data.flood_depth
    else:
        x_data = data.wind_speed

    y_data = data.damage_ratio_min + uncertainty_parameter * (data.damage_ratio_max - data.damage_ratio_min)
    y_data = np.minimum(y_data * (1 + uplift_factor), 1.0)

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
    geo_dataframe = gpd.GeoDataFrame(dataframe, geometry="geometry", crs={"init": f"epsg:{epsg}"})
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

    index_columns = [c for c in geo_dataframe.columns.values.tolist() if c != "exposure"]
    return geo_dataframe.groupby(index_columns, dropna=False)["exposure"].sum().reset_index()


def create_damage_curves(damage_data_path, damage_curve_lookup_df, uplift_factor=0, uncertainty_parameter=0):
    damage_curve_lookup_df["x_y_data"] = damage_curve_lookup_df.progress_apply(
        lambda x: get_damage_data(x, damage_data_path, uplift_factor, uncertainty_parameter),
        axis=1,
    )
    damage_curve_lookup_df[["damage_x_data", "damage_y_data"]] = damage_curve_lookup_df["x_y_data"].apply(pd.Series)
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
        dataframe[damage_ratio_columns] = dataframe[damage_ratio_columns].multiply(dataframe[damage_cost_column], axis="index")
        dataframe["damage_cost_unit"] = dataframe[cost_unit_column]
    else:
        dataframe[damage_ratio_columns] = dataframe[damage_ratio_columns].multiply(dataframe[damage_cost_column] * dataframe["exposure"], axis="index")
        cost_unit = dataframe[cost_unit_column].values.tolist()[0]
        dataframe["damage_cost_unit"] = "/".join(cost_unit.split("/")[:-1])

    return dataframe


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv",
    "-n",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the asset layer table",
)
@click.option(
    "--hazard-csv",
    "-h",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the hazard table",
)
@click.option(
    "--sensitivity-csv",
    "-s",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to sensitivity parameters table",
)
@click.option(
    "--sensitivity-id",
    "-id",
    required=True,
    type=int,
    help="Sensitivity parameter set ID",
)
@click.option(
    "--asset-gpkg-file",
    "-g",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to asset GPKG file",
)
@click.option(
    "--asset-gpkg-label",
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
    "--damage-curve-mapping-csv",
    "-m",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to table mapping hazard and asset to damage curve reference",
)
@click.option(
    "--damage-threshold-uplift-csv",
    "-t",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to the table of hazard specific thresholds and uplift factors",
)
@click.option(
    "--damage-curves-dir",
    "-d",
    required=True,
    type=click.Path(exists=True, dir_okay=True, file_okay=False, readable=True),
    help="Path to the directory containing damage curves in xlsx format",
)
@click.option(
    "--intersection",
    "-i",
    required=True,
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
    help="Path to geoparquet assets split on raster grids",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
    help="Path to save direct damages to",
)
def direct_damages(
    network_csv,
    hazard_csv,
    sensitivity_csv,
    sensitivity_id,
    asset_gpkg_file,
    asset_gpkg_label,
    asset_layer,
    damage_curve_mapping_csv,
    damage_threshold_uplift_csv,
    damage_curves_dir,
    intersection,
    output_path,
):

    sensitivity = pd.read_csv(sensitivity_csv)
    cost_uncertainty_parameter = sensitivity.loc[sensitivity_id, "cost_uncertainty_parameter"].squeeze()
    damage_uncertainty_parameter = sensitivity.loc[sensitivity_id, "damage_uncertainty_parameter"].squeeze()

    asset_data_details = pd.read_csv(network_csv)
    hazard_data_details = pd.read_csv(hazard_csv, encoding="latin1")
    damage_curve_lookup = pd.read_csv(damage_curve_mapping_csv)[["sector", "hazard_type", "asset_name", "asset_sheet"]]

    hazard_attributes = pd.read_csv(damage_threshold_uplift_csv)
    flood_hazards = hazard_attributes[hazard_attributes["hazard_type"] == "flooding"]["hazard"].values.tolist()

    logging.info("Read damage curves")
    damage_curves = []
    for idx, hazard in hazard_attributes.iterrows():

        damage_curve_df = damage_curve_lookup[damage_curve_lookup["hazard_type"] == hazard["hazard_type"]]
        damage_curve_df["hazard"] = hazard["hazard"]

        damage_curve_df = create_damage_curves(
            damage_curves_dir,
            damage_curve_df,
            uplift_factor=hazard["uplift_factor"],
            uncertainty_parameter=damage_uncertainty_parameter,
        )
        damage_curves.append(damage_curve_df)

    damage_curves = pd.concat(damage_curves, axis=0, ignore_index=True)
    del damage_curve_df

    logging.info("Preprocessing cost data")
    asset_info = asset_data_details.loc[(asset_data_details.asset_gpkg == asset_gpkg_label) & (asset_data_details.asset_layer == asset_layer)].squeeze()
    asset_sector = asset_info.sector
    asset_id = asset_info.asset_id_column
    asset_min_cost = asset_info.asset_min_cost_column
    asset_max_cost = asset_info.asset_max_cost_column
    asset_cost_unit = asset_info.asset_cost_unit_column

    asset_df = gpd.read_file(asset_gpkg_file, layer=asset_layer)
    asset_df[asset_min_cost] = asset_df.apply(
        lambda x: convert_cost_units(x, asset_min_cost, asset_cost_unit, 1.0 / jd_to_usd),
        axis=1,
    )
    asset_df[asset_max_cost] = asset_df.apply(
        lambda x: convert_cost_units(x, asset_max_cost, asset_cost_unit, 1.0 / jd_to_usd),
        axis=1,
    )
    asset_df[asset_cost_unit] = asset_df[asset_cost_unit].replace(["USD", "US$"], "J$", regex=True)
    # We just need to modify and correct energy asset costs of nodes to J$ and costs of of edges to $J/m
    if asset_sector == "energy" and asset_info.asset_layer == "edges":
        asset_df[asset_min_cost] = asset_df[asset_min_cost] / asset_df["length"]
        asset_df[asset_max_cost] = asset_df[asset_max_cost] / asset_df["length"]
        asset_df[asset_cost_unit] = "J$/m"

    asset_df["damage_cost"] = asset_df[asset_min_cost] + cost_uncertainty_parameter * (asset_df[asset_max_cost] - asset_df[asset_min_cost])
    asset_df["damage_cost"] = asset_df.progress_apply(lambda x: modify_cost_units(x, asset_cost_unit), axis=1)
    hazard_damages = []

    logging.info("Reading intersection (splits) data")
    hazard_df = gpd.read_parquet(intersection)
    hazard_df = hazard_df.to_crs(epsg=epsg_jamaica)
    hazard_df = add_exposure_dimensions(hazard_df, dataframe_type=asset_info.asset_layer, epsg=epsg_jamaica)

    logging.info("Calculating direct damages")
    for hazard_info in hazard_attributes.itertuples():
        logging.info(hazard_info.hazard)

        if getattr(asset_info, f"{hazard_info.hazard}_asset_damage_lookup_column") != "none":
            asset_hazard = getattr(asset_info, f"{hazard_info.hazard}_asset_damage_lookup_column")
            hazard_keys = hazard_data_details[hazard_data_details["hazard"] == hazard_info.hazard]["key"].values.tolist()
            hazard_effect_df = hazard_df[[asset_id, "exposure", "exposure_unit"] + hazard_keys]
            damages_df = damage_curves[(damage_curves["sector"] == asset_sector) & (damage_curves["hazard"] == hazard_info.hazard)]
            damaged_assets = list(set(damages_df["asset_name"].values.tolist()))
            affected_assets_df = asset_df[asset_df[asset_hazard].isin(damaged_assets)][[asset_id, asset_hazard, asset_cost_unit, "damage_cost"]]
            damaged_assets = list(set(affected_assets_df[asset_hazard].values.tolist()))
            damages_df = damages_df[damages_df["asset_name"].isin(damaged_assets)]
            affected_assets = list(set(affected_assets_df[asset_id].values.tolist()))

            hazard_effect_df["hazard_threshold"] = hazard_info.hazard_threshold
            if hazard_info.hazard in flood_hazards:
                hazard_effect_df[hazard_keys] = hazard_effect_df[hazard_keys] - hazard_info.hazard_threshold
                hazard_effect_df = hazard_effect_df[(hazard_effect_df[hazard_keys] > 0).any(axis=1)]
                hazard_effect_df = hazard_effect_df[hazard_effect_df[asset_id].isin(affected_assets)]
            else:
                hazard_effect_df[hazard_keys] = np.where(
                    hazard_effect_df[hazard_keys] <= hazard_info.hazard_threshold,
                    0,
                    hazard_effect_df[hazard_keys],
                )
                hazard_effect_df = hazard_effect_df[(hazard_effect_df[hazard_keys] > hazard_info.hazard_threshold).any(axis=1)]
                hazard_effect_df = hazard_effect_df[hazard_effect_df[asset_id].isin(affected_assets)]

            if len(hazard_effect_df.index) == 0:
                logging.info(f"No {hazard_info.hazard} intersections with {asset_info.asset_gpkg} {asset_info.asset_layer}")
            else:
                hazard_effect_df = pd.merge(
                    hazard_effect_df,
                    affected_assets_df,
                    how="left",
                    on=[asset_id],
                )

                for damage_info in damages_df.itertuples():
                    hazard_asset_effect_df = hazard_effect_df[hazard_effect_df[asset_hazard] == damage_info.asset_name]
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
                        hazard_asset_effect_df = estimate_direct_damage_costs_and_units(
                            hazard_asset_effect_df,
                            hazard_keys,
                            asset_cost_unit,
                            dataframe_type=asset_info.asset_layer,
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

                        hazard_asset_effect_df["damage_uncertainty_parameter"] = damage_uncertainty_parameter
                        hazard_asset_effect_df["cost_uncertainty_parameter"] = cost_uncertainty_parameter
                        hazard_damages.append(hazard_asset_effect_df)

                    del hazard_asset_effect_df
                logging.info(f"\n{hazard_effect_df.iloc[:, ::5].sum(numeric_only=True)}")
                del hazard_effect_df
        else:
            logging.info(f"{asset_info.asset_gpkg} {asset_info.asset_layer} not affected by {hazard_info.hazard}")

    hazard_damages = pd.concat(hazard_damages, axis=0, ignore_index=True).fillna(0)
    logging.info(f"Writing to {output_path}\n")
    hazard_damages.to_parquet(output_path, index=False)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    direct_damages()
