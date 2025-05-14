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
    snakemake -c1 results/hotspots/grid.tiff
    """
    input:
        script = "workflow/4b_hotspots/generate_grid.py",
        boundary = f"{DATA}/boundaries/jamaica.gpkg",
    params:
        cell_length = 1_000,
        boundary_buffer = 10_000,
    output:
        grid = f"{OUTPUT}/hotspots/grid.tiff",
    shell:
        """
        python {input.script} {input.boundary} {params.cell_length} {params.boundary_buffer} {output.grid}
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
        # hopefully this script can be made to work here without breaking changes
        # if not, duplicate and adapt it
        script = "workflow/1_damage/split_networks.py",
        networks = config["paths"]["network_layers"],
        hotspots_grid_metadata = "workflow/hotspots_layers.csv",
        grid = f"{OUTPUT}/hotspots/grid.tiff",
        gpkg = lambda wildcards: f"{DATA}/{get_asset_metadata(wildcards).path}",
    output:
        splits = f"{OUTPUT}/hotspots/splits/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
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
    snakemake -c1 results/hotspots/exposure/roads_edges_exposed_value.tiff
    """
    input:
        script = "workflow/4b_hotspots/exposure.py",
        networks = config["paths"]["network_layers"],
        splits = f"{OUTPUT}/hotspots/splits/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
    output:
        exposure = f"{OUTPUT}/hotspots/exposure/{gpkg}__{layer}.geoparquet",
    shell:
        """
        python {input.script} \
            --splits {input.splits} \
            --network-csv {input.networks} \
            --exposure {output.exposure}
        """


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
        per_cell_loss = f"{OUTPUT}/hotspots/transport/cell/{cell_id}.tiff",
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

