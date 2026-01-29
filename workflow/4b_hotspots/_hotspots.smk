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
5) Calculate the sum of EAD in each grid cell. Read from direct damages, summarise
    across sensitivity parameter sets, spatial join to hotpots grid. For each asset
    class, summed per sector, and total.
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
        resolution = config["hotspots"]["grid"]["resolution_meters"],
        boundary_buffer = config["hotspots"]["grid"]["boundary_buffer_meters"],
    output:
        grid = f"{DATA}/hotspots/grid.tiff",
    shell:
        """
        python {input.script} \\
            --boundary-path {input.boundary} \\
            --cell-length-meters {params.resolution} \\
            --boundary-buffer-meters {params.boundary_buffer} \\
            --output-path {output.grid}
        """


rule split_assets_by_hotspots_grid:
    """
    Split an asset class by the hotspots grid. This will split linestrings and
    polygons on cell boundaries, creating new split rows for assets that span
    more than one cell.

    NB: hotspots_grid_metadata must include the hazard layers for hotspot damage/risk calculations

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
        python {input.script} \\
            --network-csv {input.networks} \\
            --hazard-csv {input.hotspots_grid_metadata} \\
            --data-dir {DATA} \\
            --asset-gpkg {wildcards.gpkg} \\
            --asset-layer {wildcards.layer} \\
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
        python {input.script} \\
            --network-csv {input.networks} \\
            --splits-path {input.splits} \\
            --grid-path {input.grid} \\
            --asset-gpkg {wildcards.gpkg} \\
            --asset-layer {wildcards.layer} \\
            --output-path {output.exposure}
        """


rule hotspots_damage:
    """
    Calculate direct damages for an asset across all hazards with a given
    parameter set.

    Duplicated from `1_damage/_damage.smk::direct_damage` with different
    hazard_csv, sensitivity_parameters, hazard_intersection_file and
    output.damages.

    Test with: snakemake -c1
    results/hotspots/damages_rp/roads_edges/roads_edges_direct_damages.parquet
    """
    input:
        script = "workflow/1_damage/damage_calculations.py",
        network_csv = config["paths"]["network_layers"],
        hazard_csv = "workflow/hotspots_layers.csv",
        sensitivity_parameters = f"workflow/hotspots_sensitivity.csv",
        asset_gpkg = lambda wildcards: f"{DATA}/{get_asset_metadata(wildcards).path}",
        damage_curve_mapping = f"{DATA}/damage_curves/asset_damage_curve_mapping.csv",
        threshold_and_uplift = f"{OUTPUT}/direct_damages/hazard_damage_parameters.csv",
        damage_curves_dir = f"{DATA}/damage_curves",
        damage_curves = lambda wildcards: expand(
            f"{DATA}/damage_curves/damage_curves_{get_asset_metadata(wildcards).sector}_{{hazard_type}}.xlsx",
            hazard_type = ("TC", "flooding")
        ),
        hazard_intersection_file = f"{OUTPUT}/hotspots/splits/{{gpkg}}_splits__hazard_layers__{{layer}}.geoparquet",
    params:
        USD_per_JMD = config["economics"]["USD_per_JMD"],
        sensitivity_id = 0,
    output:
        damages = f"{OUTPUT}/hotspots/damages_rp/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages.parquet",
    shell:
        """
        python {input.script} \\
            --network-csv {input.network_csv} \\
            --hazard-csv {input.hazard_csv} \\
            --sensitivity-csv {input.sensitivity_parameters} \\
            --sensitivity-id {params.sensitivity_id} \\
            --asset-gpkg-file {input.asset_gpkg} \\
            --asset-gpkg-label {wildcards.gpkg} \\
            --asset-layer {wildcards.layer} \\
            --damage-curve-mapping-csv {input.damage_curve_mapping} \\
            --damage-threshold-uplift-csv {input.threshold_and_uplift} \\
            --damage-curves-dir {input.damage_curves_dir} \\
            --intersection {input.hazard_intersection_file} \\
            --USD-per-JMD {params.USD_per_JMD} \\
            --output-path {output.damages}
        """


rule hotspots_ead:
    """Calculate EAD from direct damages, per hazard, per asset layer. Output summed to the hotspots grid.
    """
    input:
        script = "workflow/4b_hotspots/hotspots_ead.py",
        damages = f"{OUTPUT}/hotspots/damages_rp/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages.parquet",
        grid = f"{DATA}/hotspots/grid.tiff",
    output:
        ead = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}__{{hazard}}.tiff",
    shell:
        """
        python {input.script} \\
            --splits-path {input.damages} \\
            --grid-path {input.grid} \\
            --hazard {wildcards.hazard} \\
            --output-path {output.ead}
        """


rule hotspots_ead_sector_all_flood:
    """Per asset layer, sum all flood EAD
    """
    input:
        ead_coastal = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}__coastal.tiff",
        ead_fluvial = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}__fluvial.tiff",
        ead_surface = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}__surface.tiff",
    output:
        ead_flood = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}__all_flood.tiff",
    run:
        from jamaica_infrastructure.raster import sum_rasters
        sum_rasters([input.ead_coastal, input.ead_fluvial, input.ead_surface], output.ead_flood)


rule hotspots_ead_total:
    """Per asset layer, sum all hazard EAD
    """
    input:
        ead_flood = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}__all_flood.tiff",
        ead_cyclone = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}__cyclone.tiff",
    output:
        ead_total = f"{OUTPUT}/hotspots/EAD/{{gpkg}}__{{layer}}.tiff",
    run:
        from jamaica_infrastructure.raster import sum_rasters
        sum_rasters([input.ead_flood, input.ead_cyclone], output.ead_total)


def hotspot_paths_for_ead_sector(wildcards) -> list[str]:
    """Given a sector, return the paths to all matching sector/hazard EAD files.
    """
    df = pd.read_csv(config["paths"]["network_layers"])
    asset_classes_in_sector = df[df['sector'] == wildcards.sector]
    if len(asset_classes_in_sector) == 0:
        raise ValueError(f"No assets found for {wildcards.sector=}")
    paths = []
    for asset_class in asset_classes_in_sector.itertuples():
        # as for e.g. `hotspots_exposure` rule
        paths.append(f"{OUTPUT}/hotspots/EAD/{asset_class.asset_gpkg}__{asset_class.asset_layer}__{wildcards.hazard}.tiff")
    return paths


rule hotspots_ead_hazard_by_sector:
    """Per sector, for a given hazard, sum over asset layers EAD
    """
    input:
        asset_classes = hotspot_paths_for_ead_sector
    output:
        sector_sum = f"{OUTPUT}/hotspots/EAD/{{sector}}__{{hazard}}.tiff"
    wildcard_constraints:
        sector="(water|transport|energy)",
        hazard="(cyclone|surface|fluvial|coastal|all_flood)",
    run:
        from jamaica_infrastructure.raster import sum_rasters
        sum_rasters(input.asset_classes, output.sector_sum)


rule ead_by_hazard_all_sectors:
    """Per hazard, sum over all sector EAD

    Test with:
    snakemake -c1 results/hotspots/EAD/all_sectors__all_flood.tiff
    """
    input:
        water = f"{OUTPUT}/hotspots/EAD/water__{{hazard_class}}.tiff",
        energy = f"{OUTPUT}/hotspots/EAD/energy__{{hazard_class}}.tiff",
        transport = f"{OUTPUT}/hotspots/EAD/transport__{{hazard_class}}.tiff",
    output:
        all_sector_sum = f"{OUTPUT}/hotspots/EAD/all_sectors__{{hazard_class}}.tiff"
    wildcard_constraints:
        hazard_class="(cyclone|all_flood|fluvial|surface|coastal)"
    run:
        from jamaica_infrastructure.raster import sum_rasters
        sum_rasters([input.water, input.energy, input.transport], output.all_sector_sum)


def hotspot_paths_for_sector(wildcards) -> list[str]:
    """
    Given e.g. 'energy', return the paths to all matching sector hotspot files.
    """
    df = pd.read_csv(config["paths"]["network_layers"])
    asset_classes_in_sector = df[df['sector'] == wildcards.sector]
    if len(asset_classes_in_sector) == 0:
        raise ValueError(f"No assets found for {wildcards.sector=}")
    paths = []
    for asset_class in asset_classes_in_sector.itertuples():
        # as for e.g. `hotspots_exposure` rule
        paths.append(f"{OUTPUT}/hotspots/{wildcards.hotspot_metric}/{asset_class.asset_gpkg}__{asset_class.asset_layer}.tiff")
    return paths


rule hotspots_metric_by_sector:
    """
    Sum exposed value or EAD of all asset classes for each sector

    Test with:
    snakemake -c1 results/hotspots/exposure/water.tiff
    """
    input:
        asset_classes = hotspot_paths_for_sector
    output:
        sector_sum = f"{OUTPUT}/hotspots/{{hotspot_metric}}/{{sector}}.tiff"
    wildcard_constraints:
        hotspot_metric="(exposure|EAD)"
    run:
        from jamaica_infrastructure.raster import sum_rasters
        sum_rasters(input.asset_classes, output.sector_sum)


rule hotspots_metric_all_sectors:
    """
    Sum exposed value or EAD of assets from water, energy and transport sectors.

    Test with:
    snakemake -c1 results/hotspots/exposure/all_sectors.tiff
    """
    input:
        water = f"{OUTPUT}/hotspots/{{hotspot_metric}}/water.tiff",
        energy = f"{OUTPUT}/hotspots/{{hotspot_metric}}/energy.tiff",
        transport = f"{OUTPUT}/hotspots/{{hotspot_metric}}/transport.tiff",
    output:
        all_sector_sum = f"{OUTPUT}/hotspots/{{hotspot_metric}}/all_sectors.tiff"
    wildcard_constraints:
        hotspot_metric="(exposure|EAD)"
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
        python {input.script} \\
            --labour-cost-JMD-per-hour {params.labour_cost} \\
            --trade-rerouting-fraction {params.trade_rerouting} \\
            --grid-y-index {wildcards.grid_y_index} \\
            --road-splits-path {input.road_splits} \\
            --rail-splits-path {input.rail_splits} \\
            --flow-data-dir {input.flow_data_dir} \\
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
        rerouting_loss = f"{OUTPUT}/hotspots/transport/rerouting_loss.tiff",
        isolation_loss = f"{OUTPUT}/hotspots/transport/isolation_loss.tiff",
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

        no_data_value = np.nan
        write_kwargs = {
            "driver": "GTiff",
            "height": arr.shape[0],
            "width": arr.shape[1],
            "count": 1,
            "dtype": rasterio.float32,
            "nodata": no_data_value,
            "crs": grid_dataset.crs,
            "transform": grid_dataset.transform
        }

        for variable in ("economic_loss", "isolation_loss", "rerouting_loss"):
            with rasterio.open(output[variable], "w", **write_kwargs) as output_dataset:
                arr[loss.cell_index_y, loss.cell_index_x] = np.where(loss[variable] > 0, loss[variable], no_data_value)
                output_dataset.write(arr, 1)


rule smooth_raster:
    """
    Apply quantity preserving smoothing Gaussian kernel to hotspots quantities.

    Test with:
    snakemake -c1 results/hotspots/transport/economic_loss_smoothed.tiff,
    """
    input:
        script = "workflow/4b_hotspots/kde.py",
        coarse = "{output_path}/{filename}.tiff",
    params:
        resolution = config["hotspots"]["kernel_density_estimation"]["resolution_meters"],
        bandwidth = config["hotspots"]["kernel_density_estimation"]["bandwidth_meters"],
    output:
        smoothed = "{output_path}/{filename}_smoothed.tiff",
    shell:
        """
        python {input.script} \\
            --input-raster-path {input.coarse} \\
            --output-raster-path {output.smoothed} \\
            --output-resolution {params.resolution} \\
            --bandwidth {params.bandwidth}
        """


rule target_hotspot_tiffs:
    input:
        tiffs = (
            expand(
                f"{OUTPUT}/hotspots/{{hotspot_metric}}/{{sector}}_smoothed.tiff",
                hotspot_metric=["EAD", "exposure"],
                sector=["water", "energy", "transport", "all_sectors"]
            ) +
            expand(
                f"{OUTPUT}/hotspots/EAD/{{sector}}__{{hazard}}_smoothed.tiff",
                sector=["water", "energy", "transport", "all_sectors"],
                hazard=["all_flood", "cyclone"],
            )
        )
