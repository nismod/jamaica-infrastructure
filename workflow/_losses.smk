def damage_ensemble_files(wildcards):
    # wait for ensemble parameters to be generated, then read the length of the list
    filepath = checkpoints.sensitivity_parameters.get(**wildcards).output.sensitivity_parameters
    n_ensemble = len(pd.read_csv(filepath))
    return expand(
        "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages_parameter_set_{parameter_set}.parquet",
        parameter_set=range(n_ensemble)
    )

rule collapse_sensitivity:
    """
    Summarise direct damage results (aggregate over sensitivity analysis)

    Test with:
    snakemake -c1 results/direct_damages_summary/roads_edges_damages.csv
    """
    input:
        script = "scripts/analysis/direct_damage_summarise.py",
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        damages = damage_ensemble_files
        # TODO: Need EAD_and_EAEL_parameter_set_\d+ files
    output:
        damages = "{output_path}/direct_damages_summary/{gpkg}_{layer}_damages.csv",
        exposure = "{output_path}/direct_damages_summary/{gpkg}_{layer}_exposure.csv",
    shell:
        f"""
        # build list of input damage files: --damage <damage_file_0> --damage <damage_file_1> --damage <damage_file_n>
        DAMAGE_FILES=""
        for FILE in {{input.damages}}; do
            DAMAGE_FILES="$DAMAGE_FILES --damages $FILE"
        done

        python {{input.script}} \
            --network-csv {{input.network_csv}} \
            --sensitivity-csv {{input.sensitivity_parameters}} \
            --asset-gpkg {{wildcards.gpkg}} \
            --asset-layer {{wildcards.layer}} \
            $DAMAGE_FILES \
            --output-damages-path {{output.damages}} \
            --output-exposure-path {{output.exposure}}
        """


def get_single_failure_scenarios(wildcards):
    row = get_asset_row(wildcards)
    sfs = row.single_failure_scenarios
    if sfs == "None" or sfs == "none":
        return []
    if row.sector == "buildings":
        return f"{DATA}/{sfs}"
    return f"{wildcards.output_path}/{sfs}"

rule EAD_EAEL_results:
    """
    Calculate Estimated Annual Damages for an asset across all hazards with a given parameter set.
    
    scripts/analysis/ead_eael_calculations.py
    
    This script is called by scripts/analysis/flood_changes_setup.py which assigns different input args to it for each run.
    
    TODO: Update the script to accept parameter set id rather than parameter set values directly.
    TODO: Script should default flood_protection_name to None
    
    Test with:
    snakemake -c1 results/direct_damages/roads_edges_EAD_EAEL_parameter_set_0.csv
    """
    input:
        script = "scripts/analysis/ead_eael_calculations.py",
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        hazard_csv = config["paths"]["hazard_layers"],
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
        damage_file = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_direct_damages_parameter_set_{parameter_set}.parquet",
        single_failure_scenarios = get_single_failure_scenarios,
    output:
        "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_EAD_EAEL_parameter_set_{parameter_set}.csv",
    shell:
        """
        touch {output}
        """


rule RENAME_EAD_EAEL_FILE:
    """
    There's a potential mismatch between 
    {{output_path}}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv created in loss_summary and
    {{output_path}}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.parquet required as output targets
    
    This rule renames the former to the latter.
    """
    input:
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_EAD_EAEL.csv",
    output:
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_EAD_EAEL.parquet",
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Renaming {input} to {output}"
            cp {input} {output}
        fi
        """


rule loss_summary:
    """
    Summarise all the loss files for an asset.
    
    scripts/analysis/damage_loss_summarised.py
    
    Test with:
    snakemake -c1 results/direct_damages_summary/roads_edges_losses.parquet
    """
    input:
        network_csv = config["paths"]["network_layers"],
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        direct_damage_results = expand(
            "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages_parameter_set_{parameters}.parquet",
            parameters=PARAMETER_SET_IDS,
        ),
        EAD_EAEL_damage_results = expand(
            "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_EAD_EAEL_parameter_set_{parameters}.csv",
            parameters=PARAMETER_SET_IDS,
        ),
        single_failure_scenarios = get_single_failure_scenarios,
    output:
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_losses.parquet",
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_exposures.parquet",
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_damages.parquet",
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_EAD_EAEL.csv",
    shell:
        """
        python scripts/smk-analysis/damage_loss_summarised.py \
            --network_csv {input.network_csv} \
            --sensitivity_parameters {input.sensitivity_parameters} \
            --direct_damage_results {input.direct_damage_results} \
            --EAD_EAEL_damage_results {input.EAD_EAEL_damage_results} \
            --single_failure_scenarios {input.single_failure_scenarios} \
            --output_path {wildcards.output_path}/direct_damages_summary \
            --gpkg {wildcards.gpkg} \
            --layer {wildcards.layer}
        """


