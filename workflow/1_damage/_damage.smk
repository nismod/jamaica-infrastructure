"""
What are the costs associated with rebuilding damaged assets?
"""

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
        splits = f"{OUTPUT}/hazard_asset_intersection/{{gpkg}}_splits__hazard_layers__{{layer}}.geoparquet",
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
                problem, 10, num_levels=4, optimal_trajectories=4, local_optimization=False
            )

        elif config["sensitivity_analysis"] == False:
            # one sample, no variation
            # perturbation is applied: min + sensitivity parameter * (max - min)
            # therefore a sensitivity parameter of 0.5 is the middle estimate
            values = [[0.5] * len(variables)]

        else:
            raise ValueError("config['sensitivity_analysis'] must be boolean")

        df = pd.DataFrame(values, columns=variables)
        df.index.name = "set_id"
        df.to_csv(output.sensitivity_parameters, float_format='%.3f')


rule hazard_damage_parameters:
    """
    The `direct_damage` rule requires a CSV of hazard intensity modifiers by
    hazard type.
    
    In the case of adaptation, the default values are modified by this rule. In
    both cases, a new file is written in each direct damage output folder and
    used as the source of the parameters.

    Test with:
    snakemake -c1 results/direct_damages/hazard_damage_parameters.csv
    snakemake -c1 results/flood_threshold_1p5/direct_damages/hazard_damage_parameters.csv
    snakemake -c1 results/cyclone_damage_curve_change_0p76/direct_damages/hazard_damage_parameters.csv
    """
    input:
        default_parameters = f"{DATA}/damage_curves/hazard_damage_parameters.csv",
    output:
        parameters = "{output_path}/direct_damages/hazard_damage_parameters.csv",
    run:
        from pathlib import Path 

        import pandas as pd

        df: pd.DataFrame = pd.read_csv(input.default_parameters)

        folder_name: str = str(Path(output.parameters).parent.parent)

        if (flood_match := re.search(r"flood_threshold_(\d+)p(\d+)", folder_name)):
            w, d = flood_match.groups()
            df.loc[df["hazard_type"] == "flooding", "hazard_threshold"] = float(f"{w}.{d}")
        if (cyclone_match := re.search(r"cyclone_damage_curve_change_(\d+)p(\d+)", folder_name)):
            w, d = cyclone_match.groups()
            # see https://github.com/nismod/jamaica-infrastructure/blob/v1.0/scripts/analysis/cyclone_changes_poles.py#L58
            df.loc[df["hazard_type"] == "TC", "uplift_factor"] = -1 * float(f"{w}.{d}")

        df.to_csv(output.parameters, index=False)


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
        threshold_and_uplift = "{output_path}/direct_damages/hazard_damage_parameters.csv",
        damage_curves_dir = f"{DATA}/damage_curves",
        damage_curves = lambda wildcards: expand(
            f"{DATA}/damage_curves/damage_curves_{get_asset_metadata(wildcards).sector}_{{hazard_type}}.xlsx",
            hazard_type = ("TC", "flooding", "coastal")
        ),
        hazard_intersection_file = f"{OUTPUT}/hazard_asset_intersection/{{gpkg}}_splits__hazard_layers__{{layer}}.geoparquet",
    params:
        USD_per_JMD = config["economics"]["USD_per_JMD"],
        sensitivity_id = sensitivity_id_from_slug,
        epoch = config["adaptation_options"]["projection_end_year"],
    output:
        damages = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_direct_damages_{parameter_set}.parquet",
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
