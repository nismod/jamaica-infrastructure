"""
Generate the files required by irv-jamaica/etl/adaptation_files.csv
"""

from typing import List
import pandas


timeseries_files = []
for risk_type in ["EAD", "EAEL"]:
        timeseries_files = [
            *timeseries_files,
            *[f"{risk_type}_timeseries_{val_type}" for val_type in ["min", "mean", "max"]]
        ]

rule damage_loss_timeseries_and_NPV:
    """
    Estimate the damage loss timeseries and NPV for an asset with an adaptation.
    
    scripts/analysis/damage_loss_timeseries_and_npv.py
    
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
        timeseries = [
            f"{{output_path}}/loss_damage_timeseries/{{gpkg}}_{{layer}}_{file_suffix}.csv"
            for file_suffix in timeseries_files
        ],
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


rule adaptation_options_costs:
    """
    Generate the adaptation options for each asset.
    
    scripts/analysis/adaptation_options_costs.py
    
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


def get_hazard_thresholds(hazard: str) -> List[str]:
    if hazard == "TC":
        return ["0p76"]
    elif hazard == "flooding":
        return ["1p0", "1p5", "2p0", "2p5"]
    else:
        raise ValueError(f"Unknown hazard type: {hazard}")

def get_risk_dir(hazard: str, threshold: str) -> str:
    """
    Get the risk directory for a hazard and threshold.
    """
    return f"{OUTPUT}/{'flood_threshold' if hazard == 'flooding' else 'cyclone_damage_curve_change'}_{threshold}"

risk_dirs = [OUTPUT]
for hazard in HAZARD_TYPES:
    risk_dirs = [
        *risk_dirs,
        *[get_risk_dir(hazard, threshold) for threshold in get_hazard_thresholds(hazard)]
    ]

rule benefit_cost_ratio:
    """
    Generate the benefit cost ratio for each adaptation option.

    scripts/analysis/benefit_cost_ratio_estimations.py
    
    Test with:
    snakemake -c1 results/adaptation_benefits_costs_bcr/flooding_waste_water_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv
    """
    input:
        script = "workflow/5_adaptation/benefit_cost_ratio_estimations.py",
        network_csv = config["paths"]["network_layers"],
        cost_file = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        risk_files = lambda wildcards: [f"{dir}/loss_damage_npvs/{wildcards.gpkg}_{wildcards.layer}_EAD_EAEL_npvs.csv" for dir in risk_dirs],
    output:
        bcr = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_benefits_costs_bcr.csv",
        EAD = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_costs_avoided_EAD_EAEL.csv",
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --cost-file {input.cost_file} \
            --hazard-label {hazard} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --output-dir {OUTPUT}
        """
