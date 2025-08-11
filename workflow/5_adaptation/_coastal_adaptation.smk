"""
Generate the required files for the coastal adaptation section
"""

from typing import List
import pandas

rule coastal_flood_protection_assets:
    """
    Generates the coastal protection areas and coastal protection segments

    Test with:
    snakemake -c1 results/coastal_protection_assets/Jamaica_coastal_protection_areas.gpkg
    """
    input:
        script = "workflow/5_adaptation/generate_coastal_adaptations.py",
        island_inputs = f"{DATA}/adaptation/coastal_adaptation/Foundation_Data.gpkg"
    params:
        rcp = config["coastal_adaptation"]["max_rcp"],
        rp = config["coastal_adaptation"]["max_rp"],
        epoch = config["adaptation_options"]["projection_end_year"],
        inland_buffer = config["coastal_adaptation"]["inland_buffer"],
        flood_threshold = config["coastal_adaptation"]["flood_threshold"],
        eps = config["coastal_adaptation"]["eps"],
        minpts = config["coastal_adaptation"]["minpts"],
        coast_length = config["coastal_adaptation"]["coast_length"],
        overlap = config["coastal_adaptation"]["overlap"],
    output:
        cpa_processing = f"{OUTPUT}/coastal_protection_assets/coastal_protection_processing_layers.gpkg",
        cpa_processed_data = f"{DATA}/networks/coastal_infrastructure/Jamaica_coastal_protection_areas.gpkg",
        cpa = f"{OUTPUT}/coastal_protection_assets/Jamaica_coastal_protection_areas.gpkg",
    shell:
        """
        python {input.script} \
            --island-inputs {input.island_inputs} \
            --output-dir {OUTPUT} \
            --data-dir {DATA} \
            --inland-buffer {params.inland_buffer} \
            --flood-threshold {params.flood_threshold} \
            --rcp {params.rcp} \
            --rp {params.rp} \
            --epoch {params.epoch} \
            --eps {params.eps} \
            --minpts {params.minpts} \
            --coast-length {params.coast_length} \
            --overlap {params.overlap}
        """

rule combine_coastal_protection:
    """
    Combines all the coastal protection areas for all the RCP'S & RP's into one GeoPackage

    Test with:
    snakemake -c1 results/coastal_protection_assets/combined_coastal_protection_area.gpkg
    """
    input:
        script = "workflow/5_adaptation/combine_coastal_protection.py",
        coastal_adaptation_assets = f"{OUTPUT}/coastal_protection_assets/Jamaica_coastal_protection_areas.gpkg"
    params:
        rcp = config["coastal_adaptation"]["max_rcp"],
        rp = config["coastal_adaptation"]["max_rp"],
        epoch = config["adaptation_options"]["projection_end_year"],
    output:
        combined_cpa = f"{OUTPUT}/coastal_protection_assets/combined_coastal_protection_area.gpkg"
    shell:
        """
        python {input.script} \
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \
            --rcp {params.rcp} \
            --rp {params.rp} \
            --epoch {params.epoch} \
            --output-dir {OUTPUT} \
        """


rule map_networks_to_coastal_protection:
    """
    Maps each network layer's assets to the coastal protection area its in and the height of the 
    required coastal protection feature'

    Test with:
    snakemake -c1 results/coastal_protection_assets/network_protection_mappings/electricity_network_v3.1_nodes_coastal_filtered.parquet
    """
    input:
        script = "workflow/5_adaptation/asset_to_coastal_protection_mapping.py",
        network_csv = config["paths"]["network_layers"],
        processed_data_path = config["paths"]["data"],
        coastal_adaptation_assets = f"{OUTPUT}/coastal_protection_assets/Jamaica_coastal_protection_areas.gpkg",
        combine_coastal_protection = f"{OUTPUT}/coastal_protection_assets/combined_coastal_protection_area.gpkg"
    params:
        rcp = config["coastal_adaptation"]["max_rcp"],
        rp = config["coastal_adaptation"]["max_rp"],
        epoch = config["adaptation_options"]["projection_end_year"],
    output:
        flood_asset_dict = f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/{{gpkg}}_{{layer}}_coastal_filtered.parquet"
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --processed-data-path {input.processed_data_path} \
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \
            --combine-coastal-protection {input.combine_coastal_protection} \
            --rcp {params.rcp} \
            --rp {params.rp} \
            --epoch {params.epoch} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --output-dir {OUTPUT} \
        """


def generate_coastal_protection_mapping_paths(wildcards) -> list[str]:
    metadata = pd.read_csv(config["paths"]["network_layers"])
    paths = []
    for row in metadata.itertuples():
        paths.append(
            f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/{row.asset_gpkg}_{row.asset_layer}_coastal_filtered.parquet"
        )
    return paths

rule map_networks_to_coastal_protection_all_layers:
    """
    Expand map_networks_to_coastal_protection for all network layers.

    Test with:
    snakemake -c1 results/coastal_protection_assets/network_protection_mappings/all_layers.flag
    """
    input:
        generate_coastal_protection_mapping_paths
    output:
        flag = f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/all_layers.flag"


rule count_assets_per_coastal_protection:
    """
    Counts the number of assets for each network layer that intersect with each coastal
    flood protection area and sums the costs of the network assets.

    Test with:
    snakemake -c1 results/coastal_protection_assets/coastal_protection_assets_breakdown.csv
    """
    input:
        script = "workflow/5_adaptation/count_coastal_protection_assets.py",
        network_csv = config["paths"]["network_layers"],
        processed_data_path = config["paths"]["data"],
        coastal_adaptation_assets = f"{OUTPUT}/coastal_protection_assets/Jamaica_coastal_protection_areas.gpkg"
    params:
        rcp = config["coastal_adaptation"]["max_rcp"],
        rp = config["coastal_adaptation"]["max_rp"],
        epoch = config["adaptation_options"]["projection_end_year"],
    output:
        flood_asset_count = f"{OUTPUT}/coastal_protection_assets/coastal_protection_assets_breakdown.csv"
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --processed-data-path {input.processed_data_path} \
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \
            --rcp {params.rcp} \
            --rp {params.rp} \
            --epoch {params.epoch} \
            --output-dir {OUTPUT} \
        """


rule benefit_cost_ratio_coastal_protection:
    """
    Generate the benefit cost ratio for each coastal defence adaptation option.
    
    Test with:
    snakemake -c1 results/adaptation_benefits_costs_bcr/coastal_waste_water_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv
    """
    input:
        script = "workflow/5_adaptation/benefit_cost_ratio_estimations_coastal.py",
        network_csv = config["paths"]["network_layers"],
        cost_file = f"{OUTPUT}/adaptation_costs/coastal_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        protection_asset_dict = f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/{{gpkg}}_{{layer}}_coastal_filtered.parquet",
        protection_asset_breakdown = f"{OUTPUT}/coastal_protection_assets/coastal_protection_assets_breakdown.csv",
        risk_files = get_risk_files
    params:
        projection_end_year = config["adaptation_options"]["projection_end_year"],
        disruption_duration = config["adaptation_options"]["disruption_duration_days"],
        rcp = config["coastal_adaptation"]["max_rcp"],
        rp = config["coastal_adaptation"]["max_rp"],
    wildcard_constraints:
        hazard=r"(coastal)"
    output:
        bcr = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_benefits_costs_bcr.csv",
        EAD = f"{OUTPUT}/adaptation_benefits_costs_bcr/{{hazard}}_{{gpkg}}_{{layer}}_adaptation_costs_avoided_EAD_EAEL.csv",
    shell:
        """
        python {input.script} \\
            --network-csv {input.network_csv} \\
            --cost-file {input.cost_file} \\
            --protection-asset-breakdown {input.protection_asset_breakdown} \\
            --protection-asset-dict {input.protection_asset_dict} \\
            --asset-gpkg {wildcards.gpkg} \\
            --asset-layer {wildcards.layer} \\
            --proj-end-year {params.projection_end_year} \\
            --rcp {params.rcp} \\
            --rp {params.rp} \\
            --disruption-duration-days {params.disruption_duration} \\
            --output-dir {OUTPUT}
        """


def generate_coastal_protection_BCR_paths(wildcards) -> list[str]:
    metadata = pd.read_csv(config["paths"]["network_layers"])
    paths = []
    for row in metadata.itertuples():
        paths.append(
            f"{OUTPUT}/adaptation_benefits_costs_bcr/coastal_{row.asset_gpkg}_{row.asset_layer}_adaptation_benefits_costs_bcr.csv",
        )
    return paths

rule benefit_cost_ratio_coastal_protection_all_assets:
    """
    Coastal adaptation BCR files for all protected network layers.
    """
    input:
        generate_coastal_protection_BCR_paths
    output:
        flag = f"{OUTPUT}/adaptation_benefits_costs_bcr/coastal.flag",
    shell:
        """
        touch {output.flag}
        """

