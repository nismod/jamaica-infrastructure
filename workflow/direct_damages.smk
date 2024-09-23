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
        network_layers = f"{DATA}/network_layers.csv",
        hazard_layers = f"{DATA}/hazard_layers.csv",
        script = "../scripts/exposure/split_networks.py",
    output:
        splits = directory(f"{DATA}/hazard_asset_intersection"),
    shell:
        """
        python {input.script} {input.networks} {input.hazards} {wildcards.DATA_DATA}
        """



#   rule direct_damages:
#       """
#       
#       """
#       input:
#           networks = "network_layers.csv",
#           hazards = "hazard_layers.csv",
#           asset_damage_curve_mapping = f"{DATA}/damage_curves/asset_damage_curve_mapping.csv",
#           damage_thresholds = f"{DATA}/damage_curves/hazard_damage_parameters.csv",
#       output:
#
#       script:
#           "../scripts/analysis/damage_calculations.py"



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
