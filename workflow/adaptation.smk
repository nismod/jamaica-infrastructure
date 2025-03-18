"""
Generate the files required by irv-jamaica/etl/adaptation_files.csv
"""
from typing import List
import pandas

hazard_types = [
    "TC",
    "flooding",
]

#TODO: deduplicate from damage.smk
def get_asset_row(wildcards) -> pandas.Series:
    """
    Get the path of an asset by its gpkg and layer strings.
    """
    df = pandas.read_csv(f"{DATA}/networks/network_layers_hazard_intersections_details.csv")
    row = df[(df['asset_gpkg'] == wildcards.gpkg) & (df['asset_layer'] == wildcards.layer)]
    if len(row) == 0:
        # Try the other files, and raise a warning if the value is in them
        for f in [
            f"{DATA}/networks/network_layers_hazard_intersections_details_0.csv",
            f"{DATA}/networks/network_layers_hazard_intersections_details_1.csv",
        ]:
            df = pandas.read_csv(f)
            row = df[(df['asset_gpkg'] == wildcards.gpkg) & (df['asset_layer'] == wildcards.layer)]
            if len(row) > 0:
                # print(f"WARNING: Found asset in {f}")
                break
        else:
            raise ValueError(f"No asset found for gpkg={wildcards.gpkg} and layer={wildcards.layer}")
    elif len(row) > 1:
        raise ValueError(f"Multiple assets found for gpkg={wildcards.gpkg} and layer={wildcards.layer}")
    return row.squeeze()

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

risk_dirs = []
for hazard in hazard_types:
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
        asset_data = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        cost_df = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        risk_files = lambda wildcards: [f"{dir}/loss_damage_npvs/{wildcards.gpkg}_{wildcards.layer}_EAD_EAEL_npvs.csv" for dir in risk_dirs],
    output:
        bcr = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_benefits_costs_bcr.csv",
        EAD = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_costs_avoided_EAD_EAEL.csv",
    shell:
        """
        touch {output.bcr}
        touch {output.EAD}
        """

rule adaptation_options:
    """
    Generate the adaptation options for each asset.
    
    scripts/analysis/adaptation_options.py
    
    Test with:
    snakemake -c1 results/adaptation_options/flooding_waste_water_facilities_NWC_nodes_adaptation_timeseries_and_npvs.csv
    """
    input:
        cost_df = f"{DATA}/adaptation/adaptation_options_and_costs_jamaica.xlsx",
        asset_data = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        asset_df = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",  # gpkg
    output:
        npv = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        unit_costs = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_unit_costs.csv",
    shell:
        """
        touch {output.npv}
        touch {output.unit_costs}
        """

rule damage_loss_timeseries_and_NPV:
    """
    Estimate the damage loss timeseries and NPV for an asset with an adaptation.
    
    scripts/analysis/damage_loss_timeseries_and_npv.py
    
    Test with:
    snakemake -c1 results/flood_threshold_1p0/loss_damage_npvs/waste_water_facilities_NWC_nodes_EAD_EAEL_npvs.csv
    """
    input:
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        growth_rates = f"{DATA}/macroeconomic_data/gdp_growth_rates.xlsx",
        summarised_damages = f"{OUTPUT}/{{protection_type}}_{{threshold}}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv",
    output:
        NPV = f"{OUTPUT}/{{protection_type}}_{{threshold}}/loss_damage_npvs/{{gpkg}}_{{layer}}_EAD_EAEL_npvs.csv"
    shell:
        """
        touch {output.NPV}
        """
