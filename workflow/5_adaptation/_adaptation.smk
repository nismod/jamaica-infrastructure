"""
Generate the files required by irv-jamaica/etl/adaptation_files.csv
"""


rule adaptation_options_costs:
    """
    Generate the adaptation options for each asset.
    
    Test with:
    snakemake -c1 results/adaptation_costs/flooding_costs/waste_water_facilities_NWC_nodes_adaptation_timeseries_and_npvs.csv
    """
    input:
        script = "workflow/5_adaptation/adaptation_options_costs.py",
        cost_file = f"{DATA}/adaptation/adaptation_options_and_costs_jamaica.xlsx",
        network_csv = config["paths"]["network_layers"],
        asset_file = lambda wildcards: f"{DATA}/{get_asset_metadata(wildcards).path}",  # gpkg
    params:
        baseline_year = config["adaptation_options"]["baseline_year"],
        projection_end_year = config["adaptation_options"]["projection_end_year"],
        discounting_rate = config["adaptation_options"]["discounting_rate"],
        epsg = config["adaptation_options"]["epsg_jamaica"]
    output:
        npv = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        unit_costs = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_unit_costs.csv",
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --asset-file {input.asset_file} \
            --cost-file {input.cost_file} \
            --hazard-label {wildcards.hazard} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --output-dir {OUTPUT} \
            --baseline-year {params.baseline_year} \
            --projection-end-year {params.projection_end_year} \
            --discounting-rate {params.discounting_rate} \
            --epsg {params.epsg}
        """


rule damage_loss_timeseries_and_NPV:
    """
    Estimate the damage loss timeseries and NPV for an asset with an adaptation.
    
    Test with:
    snakemake -c1 results/flood_threshold_1p0/loss_damage_npvs/waste_water_facilities_NWC_nodes_EAD_EAEL_npvs.csv
    """
    input:
        script = "workflow/5_adaptation/damage_loss_timeseries_and_npv.py",
        network_csv = config["paths"]["network_layers"],
        growth_rates = f"{DATA}/macroeconomic_data/gdp_growth_rates.xlsx",
        summarised_damages = f"{{output_path}}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv",
    params:
        baseline_year = config["adaptation_options"]["baseline_year"],
        projection_end_year = config["adaptation_options"]["projection_end_year"],
        discounting_rate = config["adaptation_options"]["discounting_rate"]
    output:
        NPV = f"{{output_path}}/loss_damage_npvs/{{gpkg}}_{{layer}}_EAD_EAEL_npvs.csv",
        timeseries = expand(
            "{{output_path}}/loss_damage_timeseries/{{gpkg}}_{{layer}}_{variable}_timeseries_{aggregate}.csv",
            variable=["EAD", "EAEL"],
            aggregate=["amin", "mean", "amax"],
        ),
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --growth-rates-xls {input.growth_rates} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --baseline-year {params.baseline_year} \
            --projection-end-year {params.projection_end_year} \
            --discounting-rate {params.discounting_rate} \
            --output-path {wildcards.output_path}
        """


def get_hazard_thresholds(hazard: str) -> list[str]:
    """Stringify the floats from config. They will be used as paths."""
    if hazard == "TC":
        return ["p".join(str(x).split(".")) for x in config["adaptation_options"]["TC_winds_damage_curve_squash"]]
    elif hazard == "flooding":
        return ["p".join(str(x).split(".")) for x in config["adaptation_options"]["flood_thresholds_meters"]]
    elif hazard == "coastal":
        return ["coastal_adaptation"]
    else:
        raise ValueError(f"Unknown hazard type: {hazard}")

def get_risk_dir(hazard: str, threshold: str) -> str:
    """Get the risk directory for a hazard and threshold."""
    if hazard == "coastal":
        return f"{OUTPUT}/{threshold}"
    else:
        return f"{OUTPUT}/{'flood_threshold' if hazard == 'flooding' else 'cyclone_damage_curve_change'}_{threshold}"

def get_risk_files(wildcards) -> list[str]:
    """Get the paths to risk files for nominal and adjusted cases for a given hazard."""
    hazard = wildcards.hazard
    risk_dirs = [*[OUTPUT], *[get_risk_dir(hazard, threshold) for threshold in get_hazard_thresholds(hazard)]]
    paths = [f"{directory}/loss_damage_npvs/{wildcards.gpkg}_{wildcards.layer}_EAD_EAEL_npvs.csv" for directory in risk_dirs]
    return paths

rule benefit_cost_ratio:
    """
    Generate the benefit cost ratio for each adaptation option.
    
    Test with:
    snakemake -c1 results/adaptation_benefits_costs_bcr/flooding_waste_water_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv
    """
    input:
        script = "workflow/5_adaptation/benefit_cost_ratio_estimations.py",
        network_csv = config["paths"]["network_layers"],
        cost_file = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        protection_asset_dict = f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/{{gpkg}}_{{layer}}_coastal_filtered.parquet",
        protection_asset_breakdown = f"{OUTPUT}/coastal_protection_assets/coastal_protection_assets_breakdown.csv",
        risk_files = get_risk_files
    params:
        projection_end_year = config["adaptation_options"]["projection_end_year"],
        disruption_duration = config["adaptation_options"]["disruption_duration_days"],
        flood_thresholds = config["adaptation_options"]["flood_thresholds_meters"],
        TC_factors = config["adaptation_options"]["TC_winds_damage_curve_squash"],
    output:
        bcr = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_benefits_costs_bcr.csv",
        EAD = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_costs_avoided_EAD_EAEL.csv",
    shell:
        """
        FLOOD_THRESHOLDS=""
        for VALUE in {params.flood_thresholds}; do
            FLOOD_THRESHOLDS="$FLOOD_THRESHOLDS --flood-defence-threshold $VALUE"
        done

        TC_FACTORS=""
        for VALUE in {params.TC_factors}; do
            TC_FACTORS="$TC_FACTORS --tc-damage-curve-factor $VALUE"
        done

        python {input.script} \
            --network-csv {input.network_csv} \
            --cost-file {input.cost_file} \
            --protection-asset-dict {input.protection_asset_dict} \
            --hazard-label {wildcards.hazard} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --proj-end-year {params.projection_end_year} \
            --disruption-duration-days {params.disruption_duration} \
            $FLOOD_THRESHOLDS \
            $TC_FACTORS \
            --output-dir {OUTPUT}
        """
