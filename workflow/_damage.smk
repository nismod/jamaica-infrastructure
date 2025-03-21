"""
Generate the files required by irv-jamaica/etl/damage_*_files.csv

REQUIRES:
- f"{DATA}/networks/network_layers_hazard_intersections_details[_#].csv" to be present and populated
- f"{DATA}/networks/transport/multi_modal_network.gpkg" to be present and populated
"""

def get_asset_row(wildcards) -> pd.Series:
    """
    Get the path of an asset by its gpkg and layer strings.
    """
    df = pd.read_csv(f"workflow/network_layers.csv")
    row = df[(df['asset_gpkg'] == wildcards.gpkg) & (df['asset_layer'] == wildcards.layer)]
    if len(row) > 1:
        raise ValueError(f"Multiple assets found for gpkg={wildcards.gpkg} and layer={wildcards.layer}")
    return row.squeeze()


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
        from jamaica_infrastructure.transform import read_transforms

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
        script = "scripts/smk-analysis/split_networks.py",
        networks = config["paths"]["network_layers"],
        hazards = config["paths"]["hazard_layers"],
        gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
    output:
        splits = "{output_path}/hazard_asset_intersection/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
    shell:
        f"""
        python {{input.script}} \
            --network-csv {{input.networks}} \
            --hazard-csv {{input.hazards}} \
            --data-dir {{DATA}} \
            --asset-gpkg {{wildcards.gpkg}} \
            --asset-layer {{wildcards.layer}} \
            --output-path {{output.splits}}
        """


checkpoint sensitivity_parameters:
    """
    Generate sensitivity parameter combinations for the damage calculations.
    
    The original sensitivity_parameters.txt file is created by half a dozen or so scripts.
    They all seem to set the same combinations.
    See scripts/analysis/damage_scripts_setup.py for an example.
    
    Here we'll create one canonical version called sensitivity_parameters.csv.
    
    Test with:
    snakemake -c1 processed_data/sensitivity_parameters.csv
    """
    input:
        config = "config.json"
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


def sensitivity_id_from_slug(wildcards):
    return wildcards.parameter_set.replace("parameter_set_", "")

rule direct_damage:
    """
    Calculate direct damages for an asset across all hazards with a given parameter set.
    
    Test with:
    snakemake -c1 results/direct_damages/roads_edges/roads_edges_direct_damages_parameter_set_0.parquet
    """
    input:
        script = "scripts/analysis/damage_calculations.py",
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        hazard_csv = config["paths"]["hazard_layers"],
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        asset_gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
        damage_curve_mapping = f"{DATA}/damage_curves/asset_damage_curve_mapping.csv",
        damage_threshold_uplift = f"{DATA}/damage_curves/hazard_damage_parameters.csv",
        damage_curves_dir = f"{DATA}/damage_curves",
        damage_curves = lambda wildcards: expand(
            f"{DATA}/damage_curves/damage_curves_{get_asset_row(wildcards).sector}_{{hazard_type}}.xlsx",
            hazard_type = HAZARD_TYPES
        ),
        hazard_intersection_file = "{output_path}/hazard_asset_intersection/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
    params:
        sensitivity_id = sensitivity_id_from_slug,
    output:
        damages = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_direct_damages_{parameter_set}.parquet",
    shell:
        f"""
        python {{input.script}} \
            --network-csv {{input.network_csv}} \
            --hazard-csv {{input.hazard_csv}} \
            --sensitivity-csv {{input.sensitivity_parameters}} \
            --sensitivity-id {{params.sensitivity_id}} \
            --asset-gpkg-file {{input.asset_gpkg}} \
            --asset-gpkg-label {{wildcards.gpkg}} \
            --asset-layer {{wildcards.layer}} \
            --damage-curve-mapping-csv {{input.damage_curve_mapping}} \
            --damage-threshold-uplift-csv {{input.damage_threshold_uplift}} \
            --damage-curves-dir {{input.damage_curves_dir}} \
            --intersection {{input.hazard_intersection_file}} \
            --output-path {{output.damages}}
        """
