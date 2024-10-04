"""
- Rasterising/splitting networks by raster grids
- Calculating direct damage estimates
"""


import pandas


hazard_layers = pandas.read_csv(config["paths"]["hazard_layers"])
network_layers = pandas.read_csv(config["paths"]["network_layers"])
bbox = config['extract_bbox']


rule rasterise_networks:
    """
    Run snail to find the intersection of networks and raster grids
    Test with:
    snakemake -c1 results/hazard_asset_intersection
    """
    input:
        networks = config["paths"]["network_layers"],
        hazards = config["paths"]["hazard_layers"],
        script = "scripts/exposure/split_networks.py",
    output:
        splits = directory(f"{OUTPUT}/hazard_asset_intersection")
    shell:
        f"""
        python {{input.script}} {{input.networks}} {{input.hazards}} {DATA} {OUTPUT}
        """


rule direct_damages:
    """
    Run direct damage calculation for all networks and hazards.
    Test with:
    snakemake -c1 results/direct_damages
    """
    input:
        networks = config["paths"]["network_layers"],
        hazards = config["paths"]["hazard_layers"],
        splits = f"{OUTPUT}/hazard_asset_intersection",
        asset_damage_curve_mapping = f"{DATA}/damage_curves/asset_damage_curve_mapping.csv",
        damage_thresholds = f"{DATA}/damage_curves/hazard_damage_parameters.csv",
        script_driver = "scripts/analysis/damage_loss_setup_script.py",
        script_core = "scripts/analysis/damage_calculations.py",
    threads:
        workflow.cores
    output:
        sensitivity_parameter_set = f"{DATA}/parameter_combinations.txt",
        problem_specification = f"{OUTPUT}/direct_damages/ead_eael_results.txt",
    shell:
        f"python {{input.script_driver}} {{input.networks}} {{input.hazards}} {{threads}}"


rule extract:
    input:
        [f"{DATA}/extract/networks/{path}" for path in network_layers.path.unique()],
        [f"{DATA}/extract/{path}" for path in hazard_layers.path]


rule extract_vector_bbox:
    input:
        f"{DATA}/{{path}}.gpkg"
    output:
        f"{DATA}/extract/{{path}}.gpkg"
    params:
        left=bbox['left'],
        right=bbox['right'],
        top=bbox['top'],
        bottom=bbox['bottom'],
    shell:
        """
        # clipsrc needs: xmin ymin xmax ymax
        ogr2ogr \
            -clipsrc {params.left} {params.bottom} {params.right} {params.top} \
            -f GPKG {output} \
            {input}
        """


rule extract_raster_bbox:
    input:
        f"{DATA}/{{path}}.tif"
    output:
        f"{DATA}/extract/{{path}}.tif"
    params:
        left=bbox['left'],
        right=bbox['right'],
        top=bbox['top'],
        bottom=bbox['bottom'],
    shell:
        """
        # target extent (-te) needs: xmin ymin xmax ymax
        gdalwarp \
            -te {params.left} {params.bottom} {params.right} {params.top} \
            -te_srs "EPSG:3448" \
            {input} \
            {output}
        """
