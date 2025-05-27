"""
What are the costs associated with rebuilding damaged assets?
"""

rule write_hazard_transforms:
    """
    Write hazard transforms to disk alongside hazard metadata. I don't think
    this is actually used?

    Test with:
    snakemake -c1 processed_data/hazards/hazard_layers__with_transforms.csv
    """
    input:
        hazard_csv = lambda wildcards: wildcards.hazard_csv + ".csv",
        data_dir = DATA,
    output:
        hazard_transforms_csv = "{hazard_csv}__with_transforms.csv"
    run:
        from jamaica_infrastructure.raster import read_transforms

        hazards = pd.read_csv(input.hazard_csv)
        hazard_transforms, transforms = read_transforms(hazards, input.data_dir)
        hazard_transforms.to_csv(output.hazard_transforms_csv, index=False)


rule rasterise_asset_layer:
    """
    Split networks into nodes, edges, and areas.
    
    Test with:
    snakemake -c1 results/hazard_asset_intersection/roads_splits__hazard_layers__edges.geoparquet
    """
    input:
        script = "workflow/1_damage/split_networks.py",
        networks = config["paths"]["network_layers"],
        hazards = config["paths"]["hazard_layers"],
        gpkg = lambda wildcards: f"{DATA}/{get_asset_metadata(wildcards).path}",
    output:
        splits = "{output_path}/hazard_asset_intersection/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
    shell:
        """
        python {input.script} \
            --network-csv {input.networks} \
            --hazard-csv {input.hazards} \
            --data-dir {DATA} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --output-path {output.splits}
        """


checkpoint sensitivity_parameters:
    """
    Generate sensitivity parameter combinations for the damage calculations.
    
    Test with:
    snakemake -c1 processed_data/sensitivity_parameters.csv
    """
    input:
        config = "config.yaml"
    output:
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv"
    run:
        import pandas as pd
        from SALib.sample import morris

        variables = ("cost_uncertainty_parameter", "damage_uncertainty_parameter")

        if config["sensitivity_analysis"] == True:
            problem = {
                "num_vars": len(variables),
                "names": variables,
                "bounds": [[0, 1.0] for var in variables],
            }
            values = morris.sample(
                problem, 4, num_levels=4, optimal_trajectories=2, local_optimization=False
            )

        elif config["sensitivity_analysis"] == False:
            values = [(1.0, 1.0)]  # one sample, no variation

        else:
            raise ValueError("config['sensitivity_analysis'] must be boolean")

        df = pd.DataFrame(values, columns=variables)
        df.index.name = "set_id"
        df.to_csv(output.sensitivity_parameters, float_format='%.3f')


rule direct_damage:
    """
    Calculate direct damages for an asset across all hazards with a given parameter set.
    
    Test with:
    snakemake -c1 results/direct_damages/roads_edges/roads_edges_direct_damages_parameter_set_0.parquet
    """
    input:
        script = "workflow/1_damage/damage_calculations.py",
        network_csv = config["paths"]["network_layers"],
        hazard_csv = config["paths"]["hazard_layers"],
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        asset_gpkg = lambda wildcards: f"{DATA}/{get_asset_metadata(wildcards).path}",
        damage_curve_mapping = f"{DATA}/damage_curves/asset_damage_curve_mapping.csv",
        damage_threshold_uplift = f"{DATA}/damage_curves/hazard_damage_parameters.csv",
        damage_curves_dir = f"{DATA}/damage_curves",
        damage_curves = lambda wildcards: expand(
            f"{DATA}/damage_curves/damage_curves_{get_asset_metadata(wildcards).sector}_{{hazard_type}}.xlsx",
            hazard_type = HAZARD_TYPES
        ),
        hazard_intersection_file = "{output_path}/hazard_asset_intersection/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
        protection_asset_dict = f"{OUTPUT}/coastal_protection_assets/network_protection_mappings"
    params:
        USD_per_JMD = config["economics"]["USD_per_JMD"],
        sensitivity_id = sensitivity_id_from_slug,
        flood_threshold = get_flood_threshold_from_path,
    output:
        damages = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_direct_damages_{parameter_set}.parquet",
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --hazard-csv {input.hazard_csv} \
            --sensitivity-csv {input.sensitivity_parameters} \
            --sensitivity-id {params.sensitivity_id} \
            --flood-threshold {params.flood_threshold} \
            --asset-gpkg-file {input.asset_gpkg} \
            --asset-gpkg-label {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --damage-curve-mapping-csv {input.damage_curve_mapping} \
            --damage-threshold-uplift-csv {input.damage_threshold_uplift} \
            --damage-curves-dir {input.damage_curves_dir} \
            --intersection {input.hazard_intersection_file} \
            --USD-per-JMD {params.USD_per_JMD} \
            --protection-asset-dict {input.protection_asset_dict} \
            --output-path {output.damages}
        """
