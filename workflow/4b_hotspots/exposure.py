import logging

import click
import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray

from jamaica_infrastructure.utils import is_sole_value


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
@click.option("--asset-gpkg", "-g", required=True, help="asset_gpkg value in the network CSV")
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
        assert is_sole_value(splits[asset_info.asset_cost_unit_column], "J$")
        splits["split_rehab_cost_J$"] = splits[asset_info.asset_mean_cost_column]
        exposure = splits.loc[:, ["split_rehab_cost_J$", "cell_index_0_x", "cell_index_0_y"]] \
            .groupby(["cell_index_0_x", "cell_index_0_y"]).sum()

    elif asset_layer == "edges":
        assert is_sole_value(splits[asset_info.asset_cost_unit_column], "J$/m")
        splits["split_length_m"] = splits.geometry.length
        splits["split_rehab_cost_J$"] = splits["split_length_m"] * splits[asset_info.asset_mean_cost_column]
        exposure = splits.loc[:, ["split_rehab_cost_J$", "cell_index_0_x", "cell_index_0_y"]] \
            .groupby(["cell_index_0_x", "cell_index_0_y"]).sum()

    elif asset_layer == "areas":
        assert is_sole_value(splits[asset_info.asset_cost_unit_column], "J$/m2")
        splits["split_area_m2"] = splits.geometry.area
        splits["split_rehab_cost_J$"] = splits["split_area_m2"] * splits[asset_info.asset_mean_cost_column]
        exposure = splits.loc[:, ["split_rehab_cost_J$", "cell_index_0_x", "cell_index_0_y"]] \
            .groupby(["cell_index_0_x", "cell_index_0_y"]).sum()

    else:
        raise NotImplementedError(f"Exposure calculation for {asset_layer=} not implemented")

    exposure = exposure.reset_index()

    # set np.nan for non-positive values, store as nodata in raster
    no_data_value = np.nan
    exposure.loc[exposure["split_rehab_cost_J$"] <= 0, "split_rehab_cost_J$"] = no_data_value

    logging.info("Reading raster grid")
    # use the hotspots grid as a template -- inherit the transform for output
    grid = rioxarray.open_rasterio(grid_path).astype(float)  # promote to float
    grid[
        {
            "band": 0,
            "x": exposure.cell_index_0_x.to_xarray(),
            "y": exposure.cell_index_0_y.to_xarray()
        }
    ] = exposure["split_rehab_cost_J$"]
    grid.name = f"{asset_gpkg}_{asset_layer}_rehab_cost_J$"
    grid = grid.rio.write_nodata(no_data_value)

    logging.info(f"Exposure:\n{grid}")

    logging.info("Write out calculated exposure as GeoTIFF")
    grid.rio.to_raster(output_path)

    return


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    exposure()
