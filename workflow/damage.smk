"""
Generate the files required by irv-jamaica/etl/adaptation_files.csv
"""
from typing import List

TARGETS = [
    "TC_electricity_network_v3.1_nodes_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_airport_polygon_areas_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_electricity_network_v3.1_nodes_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_irrigation_assets_NIC_nodes_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_port_polygon_areas_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_potable_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_rail_edges_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_roads_edges_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_roads_nodes_adaptation_costs_avoided_EAD_EAEL.csv",
    "flooding_waste_water_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv",
]

def get_hazard_thresholds(hazard: str) -> List[str]:
    if hazard == "flooding":
        return ["1p0", "1p5", "2p0", "2p5"]
    elif hazard == "TC":
        return ["0p76"]
    else:
        raise ValueError(f"Unknown hazard: {hazard}")

rule benefit_cost_ratio_estimation:
    """
    Determine the benefit-cost ratio for an asset/hazard combination.
    
    TODO: The script will need altering to handle a single hazard and asset at a time to fit Snakemake's style.
        We may also have to break up the script into smaller parts for BCR/EAD_EAEL calculations.
    
    scripts/analysis/benefit_cost_ratio_estimation.py
    
    Test with:
    snakemake -c1 results/adaptation_benefits_costs_bcr/flooding_roads_edges_adaptation_costs_avoided_EAD_EAEL.csv
    """
    input:
        asset_exposure = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        adaptation_costs = f"{OUTPUT}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        unadapted_risks = f"{OUTPUT}/loss_damage_npvs/{{gpkg}}_{{layer}}_EAD_EAEL_npvs.csv",
        adapted_risks = lambda wildcards: expand(
            f"{OUTPUT}/{{{{hazard}}}}_threshold_{{threshold}}/loss_damage_npvs/{{{{gpkg}}}}_{{{{layer}}}}_EAD_EAEL_npvs.csv",
            threshold=get_hazard_thresholds(wildcards.hazard),
        )
    output:
        bcr = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_benefits_costs_bcr.csv",
        ead_eael = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_costs_avoided_EAD_EAEL.csv",
    shell:
        """
        touch {output.bcr}
        touch {output.ead_eael}
        """

rule adaptataion_cost:
    """
    Generate the costs of adapting an asset to a hazard.
    
    scripts/analysis/adaptation_options_costs.py
    
    Test with:
    snakemake -c1 results/adaptation_costs/flooding_costs/roads_edges_adaptation_unit_costs.csv
    """
    input:
        costs = f"{DATA}/adaptation/adaptation_options_and_costs_jamaica.xlsx",
        network_layers_intersection = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
    output:
        unit_costs = f"{DATA}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_unit_costs.csv",
        timeseries_costs = f"{DATA}/adaptation_costs/{{hazard}}_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
    shell:
        """
        touch {output.costs}
        """

rule unadapted_risk:
    """
    Generate the losses if an unprotected asset is damaged or destroyed.
    
    scripts/analysis/damage_loss_timeseries_and_npv.py
    called by:
    - scripts/analysis/flooding_changes_setup.py
    - scripts/analysis/cyclone_changes_poles.py
    
    Test with:
    snakemake -c1 results/loss_damage_npvs/flooding_roads_edges_EAD_EAEL_npvs.csv
    """
    input:
        []
    output:
        risks = f"{DATA}/{{hazard}}_threshold_{{threshold}}/loss_damage_npvs/{{gpkg}}_{{layer}}_EAD_EAEL_npvs.csv",
    shell:
        """
        touch {output.risks}
        """
