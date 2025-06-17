"""
Prepare network assets (and buildings) along with damage and loss files for
export to visualisation tool.

Vector features are given a globally unique ID.
Risk calculations are saved as parquet format.
"""

import itertools
import logging
import os
import pathlib

import click
import geopandas
import numpy
import pandas
import pyarrow.parquet as pq

from jamaica_infrastructure.utils import get_asset


RISK_FILE_SUFFIXES = ('damages.parquet', 'exposures.parquet', 'losses.parquet', 'EAD_EAEL.csv')
HAZARDS = ('coastal', 'cyclone', 'fluvial', 'surface')
RCPS = ('rcp_2.6', 'rcp_4.5', 'rcp_8.5', 'rcp_baseline')
EPOCHS = ('epoch_2010', 'epoch_2030', 'epoch_2050', 'epoch_2070', 'epoch_2080', 'epoch_2100')


def get_results_fname(layer: pandas.Series, results_dir: str, suffix: str) -> pathlib.Path:
    return pathlib.Path(f"{results_dir}/direct_damages_summary/{layer.asset_gpkg}_{layer.asset_layer}_{suffix}")


def get_results_uid_fname(layer: pandas.Series, results_dir: str, suffix: str) -> pathlib.Path:
    return pathlib.Path(f"{results_dir}/direct_damages_summary_uids/{layer.asset_gpkg}_{layer.asset_layer}_{suffix}")


def get_id_lookups_fname(layer: pandas.Series, processed_data_dir: str) -> pathlib.Path:
    return pathlib.Path(f"{processed_data_dir}/networks_uids/id_lookups/{layer.asset_gpkg}_{layer.asset_layer}_ids.parquet")


def process_subset(layer, fname, base_cols, data_cols, id_lookup, hazard, rcp, epoch, suffix, results_dir):
    data = pandas.read_parquet(
        fname,
        columns=base_cols + data_cols
    )
    linked = data.set_index(layer.asset_id_column).join(id_lookup).reset_index()
    output_fname = get_results_uid_fname(layer, results_dir, f"{hazard}__{rcp}__{epoch}__{suffix}")
    logging.info(f"Writing output to {output_fname}")
    linked.to_parquet(output_fname)


def process_buildings(layer, processed_data_dir, results_dir):
    id_lookup = pandas.read_parquet(get_id_lookups_fname(layer, processed_data_dir)).set_index(layer.asset_id_column)
    for suffix in RISK_FILE_SUFFIXES:
        try:
            fname = get_results_fname(layer, results_dir, suffix)
            if 'parquet' in suffix:
                pf = pq.ParquetFile(fname)
                for hazard, rcp, epoch in itertools.product(HAZARDS, RCPS, EPOCHS):
                    base_cols = ['osm_id'] + [col for col in pf.schema.names if 'unit' in col]
                    data_cols = [col for col in pf.schema.names if hazard in col and rcp in col and epoch in col]
                    if data_cols:
                        logging.info(base_cols, hazard, rcp, epoch, len(data_cols))
                        process_subset(layer, fname, base_cols, data_cols, id_lookup, hazard, rcp, epoch, suffix, results_dir)

            elif 'csv' in suffix:
                data = pandas.read_csv(get_results_fname(layer, results_dir, suffix), dtype={'rcp': object})
                linked = data.set_index(layer.asset_id_column).join(id_lookup).reset_index()
                output_fname = get_results_uid_fname(layer, results_dir, suffix.replace('csv', 'parquet'))
                logging.info(f"Writing output to {output_fname}")
                linked.to_parquet(output_fname)

        except Exception as ex:
            raise ex


def process_layer(layer, processed_data_dir, results_dir):
    id_lookup = pandas.read_parquet(get_id_lookups_fname(layer, processed_data_dir)).set_index(layer.asset_id_column)

    for suffix in RISK_FILE_SUFFIXES:
        try:
            results_fname = get_results_fname(layer, results_dir, suffix)
            if 'parquet' in suffix:
                data = pandas.read_parquet(results_fname)
            elif 'csv' in suffix:
                data = pandas.read_csv(results_fname, dtype={'rcp': object})
            else:
                logging.info(f"WARN Skipping suffix with unhandled filetype: {suffix}")
                continue

            linked = data.set_index(layer.asset_id_column).join(id_lookup).reset_index()
            assert len(data) == len(linked), (len(data), len(linked))

            output_fname = get_results_uid_fname(layer, results_dir, suffix.replace('csv', 'parquet'))
            logging.info(f"Writing output to {output_fname}")
            linked.to_parquet(output_fname)

        except FileNotFoundError as ex:
            logging.info(ex)
            raise ex  # TODO: remove?, previously exceptions merely logged


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--network-csv", "-n", required=True, help="Path to the asset definition file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--processed-data-dir", "-p", required=True,
    type=click.Path(dir_okay=True, file_okay=False, exists=True),
    help="Designated processed data directory.",
)
@click.option(
    "--results-dir", "-r", required=True,
    type=click.Path(dir_okay=True, file_okay=False, exists=True),
    help="Designated results directory.",
)
@click.option("--asset-gpkg", "-g", required=True, help="asset_gpkg value in the network CSV")
@click.option("--asset-layer", "-l", required=True, help="asset_layer value in the network CSV")
def preprocess_for_visualisation(
    network_csv: str,
    processed_data_dir: str,
    results_dir: str,
    asset_gpkg: str,
    asset_layer: str
) -> None:
    
    if asset_gpkg == "coastal_protection_feature" and asset_layer == "edges":
        coastal = True
    else:
        coastal = False
        

    layer = get_asset(network_csv, asset_gpkg, asset_layer)
    logging.info(f"Processing {layer.sector}, {layer.asset_gpkg}, {layer.asset_layer}")

    layer_data = geopandas.read_file(os.path.join(processed_data_dir, layer.path), layer=layer.asset_layer)
    count = len(layer_data)
    layer_data['uid'] = numpy.arange(layer.base_id, layer.base_id + count)

    layer_data_output_fname = os.path.join(processed_data_dir, layer.path.replace("networks", "networks_uids"))
    if "buildings" in layer_data_output_fname:
        layer_data_output_fname = layer_data_output_fname.replace("buildings/", "networks_uids/buildings/")

    layer["count"] = count
    # add vis_path: gpkg path relative to processed_data/network_uids/network_layers.csv
    layer["vis_path"] = layer_data_output_fname.replace("networks_uids/", "")
    # output an updated fragment of the network_csv table
    # these will be concatenated by a snakemake rule if all layers are requested
    pandas.DataFrame(layer).T.to_csv(
        os.path.join(processed_data_dir, "networks_uids", f"network_layer_{asset_gpkg}_{asset_layer}.csv"),
        index=False
    )

    logging.info(f"Writing assets to disk with UIDs: {layer_data_output_fname}")
    pathlib.Path(os.path.dirname(layer_data_output_fname)).mkdir(parents=True, exist_ok=True)
    layer_data.to_file(layer_data_output_fname, layer=layer.asset_layer, index=False, driver='GPKG')

    uid_fname = get_id_lookups_fname(layer, processed_data_dir)
    uid_fname.parent.mkdir(parents=True, exist_ok=True)
    logging.info(f"Writing ID lookup: {uid_fname}")
    layer_data.loc[:, [layer.asset_id_column, 'uid']].to_parquet(uid_fname, index=False)

    if not coastal:
        logging.info("Writing results files in parquet format")
        pathlib.Path(f"{results_dir}/direct_damages_summary_uids").mkdir(parents=True, exist_ok=True)
        if 'buildings' in layer.asset_gpkg:
            process_buildings(layer, processed_data_dir, results_dir)
        else:
            process_layer(layer, processed_data_dir, results_dir)

        return


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    preprocess_for_visualisation()
