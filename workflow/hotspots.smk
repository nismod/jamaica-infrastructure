"""
The rules in this file produce raster maps of 'hotspots' of infrastructure risk.
These hotspots are computed from criticality analyses, where network assets are
removed one-at-a-time and the network performance degradation measured.

Maps are produced for each infrastructure sector.
"""


rule generate_hotspots_grid:
    """
    Create grid for hotspots analysis.
    Test with:
    snakemake -c1 results/hotspots/grid.tiff
    """
    input:
        script = "scripts/hotspots/generate_grid.py",
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

rule transport_hotspots:
    """
    Grid linestring economic losses to raster hotspots.
    Test with:
    snakemake -c1 results/hotspots/transport.tiff
    """
    input:
        script = "scripts/hotspots/rasterise_linestring_quantity.py",
        road_and_rail_losses_geometry = f"{OUTPUT}/economic_losses/single_point_failure_road_rail_edges_economic_losses.gpkg",
        grid = f"{OUTPUT}/hotspots/grid.tiff",
    params:
        variable_name = "economic_loss",
        aggregation = "sum",
    output:
        hotspots = f"{OUTPUT}/hotspots/transport.tiff",
    shell:
        """
        python \
            {input.script} \
            {input.road_and_rail_losses_geometry} \
            {input.grid} \
            {output.hotspots} \
            {params.variable_name} \
            {params.aggregation}
        """