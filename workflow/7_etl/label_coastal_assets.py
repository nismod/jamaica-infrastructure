"""
Prepare coastal defence assets for export to visualisation tool.

Vector features are given a globally unique ID.
Avoided costs are saved as parquet format.
"""

import itertools
import logging
import os
import pathlib

import click
import geopandas
import numpy
import pandas

from jamaica_infrastructure.utils import get_asset


def get_id_lookups_fname(layer: pandas.Series, processed_data_dir: str) -> pathlib.Path:
    return pathlib.Path(
        f"{processed_data_dir}/networks_uids/id_lookups/{layer.asset_gpkg}_{layer.asset_layer}_ids.parquet"
    )


def get_coastal_costs(results_dir: str, layer: pandas.Series) -> pandas.DataFrame:
    costs_path = (
        f"{results_dir}/direct_damages_summary/"
        f"{layer.asset_gpkg}_{layer.asset_layer}_EAD_EAEL.csv"
    )
    costs = pandas.read_csv(costs_path, dtype={"rcp": "str"}).set_index(layer.asset_id_column)
    return costs.loc[costs.hazard == "coastal"]


def get_protector_protectee_map(results_dir: str, layer: pandas.Series, protector_id_col: str) -> pandas.DataFrame:
    map_path: str = (
        f"{results_dir}/coastal_protection_assets/network_protection_mappings/"
        f"{layer.asset_gpkg}_{layer.asset_layer}_coastal_filtered.parquet"
    )
    id_map = pandas.read_parquet(map_path).set_index(layer.asset_id_column)
    id_map = id_map[~id_map[protector_id_col].isna()].astype({protector_id_col: int})
    return id_map


def process_EAD_EAEL(
    coastal_layer: pandas.Series,
    network_csv_path: str,
    processed_data_dir: str,
    results_dir: str
) -> None:
    """
    Sum EAD and EAEL for protected asset classes and write these to disk.

    N.B. We are overloading the meaning of the `output_path` file here.
    Normally, for other assets, it contains EAD and EAEL as a result of hazards
    to that asset. In this case, the data is the value of avoided EAD and EAEL
    at protected assets, as a result of the protecting, coastal asset. This
    allows us the use the existing visualisation infrastructure in irv-jamaica
    to show a choropleth map of total avoided damages and losses, ascribed to
    each defence asset.
    """

    # Name of column holding integer ID of coastal defence asset
    # TODO: Extract magic string?
    protector_id_col = "flood_id_rcp_852100_rp_100"

    logging.info("Finding avoided costs")

    # Loop through protectee layers, find their avoided costs (damages and losses)
    # and which coastal defence asset protects them
    data = []
    non_cost_cols = ["hazard", "rcp", "epoch"]
    cost_cols = [f"{var}_undefended_{agg}" for (var, agg) in itertools.product(("EAD", "EAEL"), ("amin", "mean", "amax"))]
    layers: pandas.DataFrame = pandas.read_csv(network_csv_path)
    for layer in layers.itertuples():
        logging.info(f"Layer: {layer.asset_gpkg}, {layer.asset_layer}")

        costs = get_coastal_costs(results_dir, layer)
        costs = costs.loc[:, non_cost_cols + cost_cols]
        if costs.empty:
            logging.info("Costs empty, skipping...\n")
            continue

        id_map = get_protector_protectee_map(results_dir, layer, protector_id_col)

        join: pandas.DataFrame = costs.join(id_map.loc[:, [protector_id_col]])
        join.index = join.index.rename("asset_id")
        logging.info(f"\n{join}")
        data.append(join)

    # Sum avoided costs for each protector asset, for each scenario
    df = pandas.concat(data)
    df = df.groupby([protector_id_col] + non_cost_cols).sum()
    df = df.reset_index().rename(columns={protector_id_col: "id"}).astype({"id": int})

    # Tag protector assets with unique database id
    id_to_uid = pandas.read_parquet(get_id_lookups_fname(coastal_layer, processed_data_dir))
    df = df.set_index("id").join(id_to_uid.set_index("id")).reset_index()

    output_path: str = (
        f"{results_dir}/direct_damages_summary_uids/"
        f"{coastal_layer.asset_gpkg}_{coastal_layer.asset_layer}_EAD_EAEL.parquet"
    )
    logging.info(f"Writing avoided costs in parquet format: {output_path}")
    logging.info(f"\n{df}")
    df.to_parquet(output_path)

    return


def process_damage_loss_and_exposure(
    layer: pandas.Series,
    processed_data_dir: str,
    results_dir: str
) -> None:
    """
    Read files from direct_damages_summary/
    Tag assets with their unique database ID
    Write to direct_damages_summary_uids/
    """

    id_to_uid = pandas.read_parquet(get_id_lookups_fname(layer, processed_data_dir))

    for kind in ("damages", "losses", "exposures"):
        df = pandas.read_parquet(
            f"{results_dir}/direct_damages_summary/{layer.asset_gpkg}_{layer.asset_layer}_{kind}.parquet"
        )

        # Tag protector assets with unique database id
        df = df.set_index("id").join(id_to_uid.set_index("id")).reset_index()

        output_path: str = (
            f"{results_dir}/direct_damages_summary_uids/"
            f"{layer.asset_gpkg}_{layer.asset_layer}_{kind}.parquet"
        )
        logging.info(f"Writing {kind} in parquet format: {output_path}")
        logging.info(f"\n{df}")
        df.to_parquet(output_path)

    return


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--coastal-network-csv",
    "-c",
    required=True,
    help="Path to the coastal layers asset definition file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--network-csv",
    "-n",
    required=True,
    help="Path to the asset definition file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--processed-data-dir",
    "-p",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False, exists=True),
    help="Designated processed data directory.",
)
@click.option(
    "--results-dir",
    "-r",
    required=True,
    type=click.Path(dir_okay=True, file_okay=False, exists=True),
    help="Designated results directory.",
)
@click.option("--asset-gpkg", "-g", required=True, help="asset_gpkg value in the network CSV")
@click.option("--asset-layer", "-l", required=True, help="asset_layer value in the network CSV")
def preprocess_for_visualisation(
    coastal_network_csv: str,
    network_csv: str,
    processed_data_dir: str,
    results_dir: str,
    asset_gpkg: str,
    asset_layer: str,
) -> None:

    layer = get_asset(coastal_network_csv, asset_gpkg, asset_layer)
    logging.info(f"Processing {layer.sector}, {layer.asset_gpkg}, {layer.asset_layer}")

    layer_data = geopandas.read_file(os.path.join(processed_data_dir, layer.path), layer=layer.asset_layer)
    count = len(layer_data)
    layer_data["uid"] = numpy.arange(layer.base_id, layer.base_id + count)
    layer["count"] = count
    # Add vis_path: gpkg path relative to processed_data/network_uids/network_layers.csv
    layer_data_output_fname = os.path.join(processed_data_dir, layer.path.replace("networks", "networks_uids"))
    layer["vis_path"] = layer_data_output_fname.replace("networks_uids/", "")
    # output an updated fragment of the network_csv table
    pandas.DataFrame(layer).T.to_csv(
        os.path.join(
            processed_data_dir,
            "networks_uids",
            f"network_layer_{asset_gpkg}_{asset_layer}.csv",
        ),
        index=False,
    )

    logging.info(f"Writing assets to disk with UIDs: {layer_data_output_fname}")
    pathlib.Path(os.path.dirname(layer_data_output_fname)).mkdir(parents=True, exist_ok=True)
    layer_data.to_file(layer_data_output_fname, layer=layer.asset_layer, index=False, driver="GPKG")

    uid_fname = get_id_lookups_fname(layer, processed_data_dir)
    uid_fname.parent.mkdir(parents=True, exist_ok=True)
    logging.info(f"Writing ID lookup: {uid_fname}")
    layer_data.loc[:, [layer.asset_id_column, "uid"]].to_parquet(uid_fname, index=False)

    pathlib.Path(f"{results_dir}/direct_damages_summary_uids").mkdir(parents=True, exist_ok=True)
    process_EAD_EAEL(layer, network_csv, processed_data_dir, results_dir)

    process_damage_loss_and_exposure(layer, processed_data_dir, results_dir)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    preprocess_for_visualisation()
