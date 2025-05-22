"""
The rules in this file produce raster maps of 'hotspots' of infrastructure risk.

1) Define a raster hotspots grid
2) Split assets on this grid
2) Calculate the sum of asset value in each grid cell (exposure). This will
    require looking up the cost column from the networks table and multiplying by
    the appropriate exposure dimension, e.g. split length for edges.
3) Sum exposure results over asset classes for each sector
4) Perform a transport criticality assessment where assets within a cell are
    failed and calculate the resulting wider economic losses
"""


checkpoint generate_hotspots_grid:
    """
    Create grid for hotspots analysis.

    Test with:
    snakemake -c1 processed_data/hotspots/grid.tiff
    """
    input:
        script = "workflow/4b_hotspots/generate_grid.py",
        boundary = f"{DATA}/boundaries/jamaica.gpkg",
    params:
        # must be in units of boundary CRS
        cell_length = config["hotspots"]["grid"]["cell_length_meters"],
        boundary_buffer = config["hotspots"]["grid"]["boundary_buffer_meters"],
    output:
        grid = f"{DATA}/hotspots/grid.tiff",
    shell:
        """
        python {input.script} \
            --boundary-path {input.boundary} \
            --cell-length-meters {params.cell_length} \
            --boundary-buffer-meters {params.boundary_buffer} \
            --output-path {output.grid}
        """


rule split_assets_by_hotspots_grid:
    """
    Split an asset class by the hotspots grid. This will split linestrings and
    polygons on cell boundaries, creating new split rows for assets that span
    more than one cell.

    Test with:
    snakemake -c1 results/hotspots/splits/roads_splits__hazard_layers__edges.geoparquet
    """
    input:
        script = "workflow/1_damage/split_networks.py",
        networks = config["paths"]["network_layers"],
        hotspots_grid_metadata = "workflow/hotspots_layers.csv",
        # split_networks.py assumes file paths (in hotspots_grid_metadata) are in DATA
        grid = f"{DATA}/hotspots/grid.tiff",
        gpkg = lambda wildcards: f"{DATA}/{get_asset_metadata(wildcards).path}",
    output:
        splits = f"{OUTPUT}/hotspots/splits/{{gpkg}}_splits__hazard_layers__{{layer}}.geoparquet",
    shell:
        """
        python {input.script} \
            --network-csv {input.networks} \
            --hazard-csv {input.hotspots_grid_metadata} \
            --data-dir {DATA} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --output-path {output.splits}
        """


rule hotspots_exposure:
    """
    Calculate value of asset class within grid cell (multiply exposure
    dimension by rebuild cost). Rebuild cost column names are available per
    asset class in networks table under `asset_mean_cost_column`.

    Test with:
    snakemake -c1 results/hotspots/exposure/roads__edges.tiff
    """
    input:
        script = "workflow/4b_hotspots/exposure.py",
        networks = config["paths"]["network_layers"],
        splits = f"{OUTPUT}/hotspots/splits/{{gpkg}}_splits__hazard_layers__{{layer}}.geoparquet",
        grid = f"{DATA}/hotspots/grid.tiff",
    output:
        exposure = f"{OUTPUT}/hotspots/exposure/{{gpkg}}__{{layer}}.tiff",
    shell:
        """
        python {input.script} \
            --network-csv {input.networks} \
            --splits-path {input.splits} \
            --grid-path {input.grid} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --output-path {output.exposure}
        """


def exposure_paths_for_sector(wildcards) -> list[str]:
    """
    Given e.g. 'energy', return the paths to all matching sector exposre files.
    """
    df = pd.read_csv(config["paths"]["network_layers"])
    asset_classes_in_sector = df[df['sector'] == wildcards.sector]
    if len(asset_classes_in_sector) == 0:
        raise ValueError(f"No assets found for {wildcards.sector=}")
    exposure_paths = []
    for asset_class in asset_classes_in_sector.itertuples():
        # as for `hotspots_exposure` rule
        exposure_paths.append(f"{OUTPUT}/hotspots/exposure/{asset_class.asset_gpkg}__{asset_class.asset_layer}.tiff")
    return exposure_paths

rule hotspots_exposure_by_sector:
    """
    Sum exposed value of all asset classes for each sector

    Test with:
    snakemake -c1 results/hotspots/exposure/water.tiff
    """
    input:
        asset_classes = exposure_paths_for_sector
    output:
        sector_sum = f"{OUTPUT}/hotspots/exposure/{{sector}}.tiff"
    run:
        from jamaica_infrastructure.raster import sum_rasters

        sum_rasters(input.asset_classes, output.sector_sum)


rule hotspots_exposure_all_sectors:
    """
    Sum exposed value of assets from water, energy and transport sectors.

    Test with:
    snakemake -c1 results/hotspots/exposure/all_sectors.tiff
    """
    input:
        water = f"{OUTPUT}/hotspots/exposure/water.tiff",
        energy = f"{OUTPUT}/hotspots/exposure/energy.tiff",
        transport = f"{OUTPUT}/hotspots/exposure/transport.tiff",
    output:
        all_sector_sum = f"{OUTPUT}/hotspots/exposure/all_sectors.tiff"
    run:
        from jamaica_infrastructure.raster import sum_rasters

        sum_rasters([input.water, input.energy, input.transport], output.all_sector_sum)


rule transport_hotspots_economic_loss_chunked_by_latitude_slice:
    """
    Find economic losses associated with the loss of road and rail edges.
    Chunked on hotspots grid row.

    Test with:
    snakemake -c1 results/hotspots/transport/y/12.parquet
    """
    input:
        script = "workflow/4b_hotspots/transport.py",
        road_splits = f"{OUTPUT}/hotspots/splits/roads_splits__hazard_layers__edges.geoparquet",
        rail_splits = f"{OUTPUT}/hotspots/splits/rail_splits__hazard_layers__edges.geoparquet",
        flow_data_dir = f"{OUTPUT}/transport_failures/nominal",
    params:
        labour_cost = config["economics"]["labour_cost_JMD_per_hour"],
        trade_rerouting = config["economics"]["disrupted_trade_fraction"],
    output:
        slice_loss = f"{OUTPUT}/hotspots/transport/y/{{grid_y_index}}.parquet",
    shell:
        """
        python {input.script} \
            --labour-cost-JMD-per-hour {params.labour_cost} \
            --trade-rerouting-fraction {params.trade_rerouting} \
            --grid-y-index {wildcards.grid_y_index} \
            --road-splits-path {input.road_splits} \
            --rail-splits-path {input.rail_splits} \
            --flow-data-dir {input.flow_data_dir} \
            --output-path {output.slice_loss}
        """


def hotspots_grid_latitude_slices(*args, **kwargs) -> list[str]:
    import rasterio
    import numpy as np
    grid_path = checkpoints.generate_hotspots_grid.get().output.grid
    with rasterio.open(grid_path) as grid_dataset:
        arr: np.ndarray[float] = grid_dataset.read()[0]
        y, x = arr.shape
    return expand(f"{OUTPUT}/hotspots/transport/y/{{grid_y_index}}.parquet", grid_y_index=range(0, y))

rule economic_loss_transport_hotspots:
    """
    Output per-cell economic loss hotspots results as single raster.

    Test with:
    snakemake -c1 results/hotspots/transport/economic_loss.tiff
    """
    input:
        latitude_slices = hotspots_grid_latitude_slices,
        labour_flows = f"{OUTPUT}/flow_mapping/labour_trips_and_activity.gpq",
        grid = f"{DATA}/hotspots/grid.tiff",
    params:
        labour_cost = config["economics"]["labour_cost_JMD_per_hour"]
    output:
        economic_loss = f"{OUTPUT}/hotspots/transport/economic_loss.tiff",
    run:
        import numpy as np
        import pandas as pd
        import rasterio

        labour_flows = pd.read_parquet(input.labour_flows, columns=["edge_id", "working_trips", "GDP_to_trips"])
        rerouting = pd.concat([pd.read_parquet(path) for path in input.latitude_slices])
        df = pd.merge(rerouting, labour_flows, how="left", on=["edge_id"]).fillna(0)

        df["mean_labour_rerouting_loss"] = params.labour_cost * df["mean_trip_time_loss"] * df["working_trips"]
        df["labour_gdp_loss"] = df["no_access"] * df["GDP_to_trips"]
        df["rerouting_loss"] = (1 - df["no_access"]) * (
            df["mean_labour_rerouting_loss"] + df["trade_rerouting_loss"]
        )
        df["isolation_loss"] = df["no_access"] * (df["labour_gdp_loss"] + df["trade_loss"])
        df["economic_loss"] = df["rerouting_loss"] + df["isolation_loss"]
        df["loss_unit"] = "J$/day"

        index_cols = ["cell_index_y", "cell_index_x"]
        loss = df.loc[:, index_cols + ["rerouting_loss", "isolation_loss", "economic_loss"]] \
            .groupby(index_cols).sum().reset_index()

        with rasterio.open(input.grid) as grid_dataset:
            arr: np.ndarray[float] = grid_dataset.read()[0].astype(np.float32)
            arr[loss.cell_index_y, loss.cell_index_x] = loss.economic_loss

            with rasterio.open(
                output.economic_loss,
                "w",
                driver="GTiff",
                height=arr.shape[0],
                width=arr.shape[1],
                count=1,
                dtype=rasterio.float32,
                crs=grid_dataset.crs,
                transform=grid_dataset.transform
            ) as output_dataset:
                output_dataset.write(arr, 1)