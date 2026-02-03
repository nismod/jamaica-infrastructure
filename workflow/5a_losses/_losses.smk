"""
What are the wider economic losses associated with failure of assets?
"""


def get_single_failure_scenario_file(wildcards):
    row = get_asset_metadata(wildcards)
    sfs = row.single_failure_scenarios
    if sfs.lower() == "none":
        return []
    if row.sector == "buildings":
        return f"{DATA}/{sfs}"
    return f"{OUTPUT}/{sfs}"

rule EAD_EAEL:
    """
    Calculate Estimated Annual Damages and Expected Annual Economic Losses for
    assets across all hazards with a given parameter set.
    
    Test with:
    snakemake -c1 results/direct_damages/airport_polygon_areas/airport_polygon_areas_EAD_EAEL_parameter_set_0.csv
    """
    input:
        script = "workflow/5a_losses/expected_damages_losses_calculations.py",
        network_csv = config["paths"]["network_layers"],
        hazard_csv = config["paths"]["hazard_layers"],
        gpkg = lambda wildcards: f"{DATA}/{get_asset_metadata(wildcards).path}",
        damage_file = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_direct_damages_{parameter_set}.parquet",
        single_failure_scenarios = get_single_failure_scenario_file,
    params:
        bridge_flood_design_RP_years = config["damages"]["bridge_flood_design_RP_years"]
    output:
        EAD_EAEL = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_EAD_EAEL_{parameter_set}.csv",
    shell:
        """
        python {input.script} \\
            --network-csv {input.network_csv} \\
            --hazard-csv {input.hazard_csv} \\
            --asset-gpkg-label {wildcards.gpkg} \\
            --asset-layer {wildcards.layer} \\
            --damage-file {input.damage_file} \\
            --single-failure-scenarios {input.single_failure_scenarios} \\
            --bridge-flood-design-rp-years {params.bridge_flood_design_RP_years} \\
            --output-path {output.EAD_EAEL}
        """


def damage_ensemble_files(wildcards):
    return expand(
        "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages_parameter_set_{parameter_set}.parquet",
        parameter_set=range(get_n_ensemble(wildcards))
    )

def EAD_EAEL_ensemble_files(wildcards):
    return expand(
        "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_EAD_EAEL_parameter_set_{parameters}.csv",
        parameters=range(get_n_ensemble(wildcards))
    )

rule collapse_sensitivity:
    """
    Aggregate over sensitivity ensemble members.

    This rule can require a lot of memory. For the buildings asset class (~1M
    assets), you'll need ~200GB RAM.
    
    Test with:
    snakemake -c1 results/direct_damages_summary/roads_edges_losses.parquet
    """
    input:
        script = "workflow/5a_losses/damage_loss_summarised.py",
        network_csv = config["paths"]["network_layers"],
        direct_damages = damage_ensemble_files,
        EAD_EAEL = EAD_EAEL_ensemble_files,
        single_failure_scenarios = get_single_failure_scenario_file,
    output:
        losses = "{output_path}/direct_damages_summary/{gpkg}_{layer}_losses.parquet",
        exposures = "{output_path}/direct_damages_summary/{gpkg}_{layer}_exposures.parquet",
        damages = "{output_path}/direct_damages_summary/{gpkg}_{layer}_damages.parquet",
        EAD_EAEL = "{output_path}/direct_damages_summary/{gpkg}_{layer}_EAD_EAEL.csv",
    shell:
        """
        DAMAGE_FILES=""
        for FILE in {input.direct_damages}; do
            DAMAGE_FILES="$DAMAGE_FILES --damages $FILE"
        done

        EAD_EAEL_FILES=""
        for FILE in {input.EAD_EAEL}; do
            EAD_EAEL_FILES="$EAD_EAEL_FILES --ead-eael $FILE"
        done

        python {input.script} \\
            --network-csv {input.network_csv} \\
            $DAMAGE_FILES \\
            $EAD_EAEL_FILES \\
            --single-failure-scenarios {input.single_failure_scenarios} \\
            --asset-gpkg {wildcards.gpkg} \\
            --asset-layer {wildcards.layer} \\
            --output-exposures {output.exposures} \\
            --output-damages {output.damages} \\
            --output-losses {output.losses} \\
            --output-ead-eael {output.EAD_EAEL}
        """


