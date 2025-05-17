"""
The rules in this file produce raster maps of 'hotspots' of infrastructure risk.

1) Define a raster hotspots grid
2) Split assets on this grid
2) Calculate the sum of asset value in each grid cell (exposure). This will
    require looking up the cost column from the networks table and multiplying by
    the appropriate exposure dimension, e.g. split length for edges.
3) Perform a criticality assessment where all assets within a cell are failed
    and calculate the resulting wider economic losses
4) Sum results over asset classes for each sector
"""


rule generate_hotspots_grid:
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


rule economic_loss_transport_hotspots_per_cell:
    """
    Find economic losses associated with the loss of road and rail edges

    Test with:
    snakemake -c1 results/hotspots/transport/cell/0.tiff
    """
    input:
        script = "?",  # an adapted version of `workflow/3_criticality/single_link_failure.py`
        airport_areas = f"{OUTPUT}/hotspots/splits/airport_polygon_splits__hazard_layers_areas.geoparquet",
        port_areas = f"{OUTPUT}/hotspots/splits/port_polygon_splits__hazard_layers_areas.geoparquet",
        road_edges = f"{OUTPUT}/hotspots/splits/rail_edges_splits__hazard_layers__edges.geoparquet",
        rail_edges = f"{OUTPUT}/hotspots/splits/roads_edges_splits__hazard_layers__edges.geoparquet",
        flow_data_dir = f"{OUTPUT}/transport_failures/nominal",
    output:
        per_cell_loss = f"{OUTPUT}/hotspots/transport/cell/{{cell_id}}.tiff",
    shell:
        """
        python {input.script} \
            --airport-splits {input.airport_areas} \
            --port-splits {input.port_areas} \
            --road-splits {input.road_edges} \
            --rail-splits {input.rail_edges} \
            --flow-data-dir {input.flow_data_dir} \
            --cell-id {wildcards.cell_id} \
            --output {output.per_cell_loss}
        """


rule economic_loss_transport_hotspots:
    """
    Output per-cell economic loss hotspots results as single raster.
    """
    input:
        script = "?",
        cells = "?",  # expand call over f"{OUTPUT}/hotspots/transport/cell/{cell_id}.tiff", for all cell_id in grid
    output:
        economic_loss = f"{OUTPUT}/hotspots/transport/economic_loss.tiff",
    shell:
        """
        ?
        """

