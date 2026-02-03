import logging

import click
import geopandas as gpd
import numpy as np
import pandas as pd

from jamaica_infrastructure.utils import is_sole_value
from jamaica_infrastructure.raster import write_indexed_splits_to_tiff


@click.command()
@click.version_option("1.0")
@click.option(
    "--network-csv", "-n", required=True, help="Path to the asset definition file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--splits-path", "-s", required=True, help="Input splits file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--grid-path", "-g", required=True, help="Hotspots raster grid",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option("--asset-gpkg", "-G", required=True, help="asset_gpkg value in the network CSV")
@click.option("--asset-layer", "-l", required=True, help="asset_layer value in the network CSV")
@click.option(
    "--output-path", "-o", required=True, help="Path to write the output exposure file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
)
def exposure(network_csv: str, splits_path: str, grid_path: str, asset_gpkg: str, asset_layer: str, output_path: str) -> None:
    
    logging.info("Reading asset metadata")
    asset_data_details = pd.read_csv(network_csv)
    asset_info = asset_data_details.loc[
        (asset_data_details.asset_gpkg == asset_gpkg) & (asset_data_details.asset_layer == asset_layer)
    ].squeeze()
    if asset_info.empty:
        raise ValueError(f"No asset data found for {asset_gpkg=} and {asset_layer=}")

    logging.info("Reading splits")
    splits: gpd.GeoDataFrame = gpd.read_parquet(splits_path)
    splits.geometry = splits.geometry.to_crs(epsg=3448)  # ensure CRS is projected, in meters

    logging.info("Calculating per-cell exposure value in J$")
    logging.info(f"{asset_layer=}")

    if asset_layer == "nodes":
        assert is_sole_value(
            splits[asset_info.asset_cost_unit_column], "J$"
        ), f"Expected J$ but got {set(splits[asset_info.asset_cost_unit_column].unique())}"
        splits["split_rehab_cost_J$"] = splits[asset_info.asset_mean_cost_column]

    elif asset_layer == "edges":
        assert is_sole_value(
            splits[asset_info.asset_cost_unit_column], "J$/m"
        ), f"Expected J$/m but got {set(splits[asset_info.asset_cost_unit_column].unique())}"
        splits["split_length_m"] = splits.geometry.length
        splits["split_rehab_cost_J$"] = (
            splits["split_length_m"] * splits[asset_info.asset_mean_cost_column]
        )

    elif asset_layer == "areas":
        assert is_sole_value(
            splits[asset_info.asset_cost_unit_column], "J$/m2"
        ), f"Expected J$/m2 but got {set(splits[asset_info.asset_cost_unit_column].unique())}"
        splits["split_area_m2"] = splits.geometry.area
        splits["split_rehab_cost_J$"] = (
            splits["split_area_m2"] * splits[asset_info.asset_mean_cost_column]
        )
    else:
        raise NotImplementedError(
            f"Exposure calculation for {asset_layer=} not implemented"
        )

    write_indexed_splits_to_tiff(
        splits,
        grid_path=grid_path,
        output_path=output_path,
        grid_name=f"{asset_gpkg}_{asset_layer}_rehab_cost_J$",
        value_colname="split_rehab_cost_J$",
        index_x="cell_index_0_x",
        index_y="cell_index_0_y",
        no_data_value=np.nan,
    )

    return


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    exposure()
