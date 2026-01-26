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
        script = "workflow/5b_coastal_adaptation/1_generate_defence_zones.py",
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
        python {input.script} \\
            --island-inputs {input.island_inputs} \\
            --output-dir {OUTPUT} \\
            --data-dir {DATA} \\
            --inland-buffer {params.inland_buffer} \\
            --flood-threshold {params.flood_threshold} \\
            --rcp {params.rcp} \\
            --rp {params.rp} \\
            --epoch {params.epoch} \\
            --eps {params.eps} \\
            --minpts {params.minpts} \\
            --coast-length {params.coast_length} \\
            --overlap {params.overlap}
        """


rule union_coastal_protection_flood_areas:
    """
    Combines all the coastal protection areas for all the RCP'S & RP's into one GeoPackage

    Test with:
    snakemake -c1 results/coastal_protection_assets/combined_coastal_protection_area.gpkg
    """
    input:
        script = "workflow/5b_coastal_adaptation/2_union_flood_areas.py",
        coastal_adaptation_assets = f"{OUTPUT}/coastal_protection_assets/Jamaica_coastal_protection_areas.gpkg"
    params:
        rcp = config["coastal_adaptation"]["max_rcp"],
        rp = config["coastal_adaptation"]["max_rp"],
        epoch = config["adaptation_options"]["projection_end_year"],
    output:
        combined_cpa = f"{OUTPUT}/coastal_protection_assets/combined_coastal_protection_area.gpkg"
    shell:
        """
        python {input.script} \\
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \\
            --rcp {params.rcp} \\
            --rp {params.rp} \\
            --epoch {params.epoch} \\
            --output-path {output.combined_cpa}
        """


rule map_networks_to_coastal_protection:
    """
    Maps each network layer's assets to the coastal protection area its in and the height of the 
    required coastal protection feature'

    Test with:
    snakemake -c1 results/coastal_protection_assets/network_protection_mappings/electricity_network_v3.1_nodes_coastal_filtered.parquet
    """
    input:
        script = "workflow/5b_coastal_adaptation/3_asset_to_zone_mapping.py",
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
        python {input.script} \\
            --network-csv {input.network_csv} \\
            --processed-data-path {input.processed_data_path} \\
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \\
            --combine-coastal-protection {input.combine_coastal_protection} \\
            --rcp {params.rcp} \\
            --rp {params.rp} \\
            --epoch {params.epoch} \\
            --asset-gpkg {wildcards.gpkg} \\
            --asset-layer {wildcards.layer} \\
            --output-dir {OUTPUT}
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
        script = "workflow/5b_coastal_adaptation/4_count_protected_assets.py",
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
        python {input.script} \\
            --network-csv {input.network_csv} \\
            --processed-data-path {input.processed_data_path} \\
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \\
            --rcp {params.rcp} \\
            --rp {params.rp} \\
            --epoch {params.epoch} \\
            --output-dir {OUTPUT}
        """


rule benefit_cost_ratio_coastal_protection:
    """
    Generate the asset-wise (e.g. road edge 7, road edge 8, airport 2, etc.)
    benefit cost ratio for every asset defended by a coastal defence adaptation
    option.
    
    Test with:
    snakemake -c1 results/adaptation_benefits_costs_bcr/coastal_waste_water_facilities_NWC_nodes_adaptation_costs_avoided_EAD_EAEL.csv
    """
    input:
        script = "workflow/5b_coastal_adaptation/5_benefit_cost_ratio.py",
        network_csv = config["paths"]["network_layers"],
        cost_file = f"{OUTPUT}/adaptation_costs/coastal_costs/{{gpkg}}_{{layer}}_adaptation_timeseries_and_npvs.csv",
        protection_asset_dict = f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/{{gpkg}}_{{layer}}_coastal_filtered.parquet",
        protection_asset_breakdown = f"{OUTPUT}/coastal_protection_assets/coastal_protection_assets_breakdown.csv",
        no_adapt_risk = f"{OUTPUT}/loss_damage_npvs/{{gpkg}}_{{layer}}_EAD_EAEL_npvs.csv",
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
            --no-adapt-risk-file {input.no_adapt_risk} \\
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

def generate_coastal_protection_cost_paths(wildcards) -> list[str]:
    metadata = pd.read_csv(config["paths"]["network_layers"])
    paths = []
    for row in metadata.itertuples():
        paths.append(
            f"{OUTPUT}/adaptation_costs/coastal_costs/{row.asset_gpkg}_{row.asset_layer}_adaptation_timeseries_and_npvs.csv",
        )
    return paths

def generate_coastal_protection_avoided_EAD_EAEL_paths(wildcards) -> list[str]:
    metadata = pd.read_csv(config["paths"]["network_layers"])
    paths = []
    for row in metadata.itertuples():
        paths.append(
            f"{OUTPUT}/adaptation_benefits_costs_bcr/coastal_{row.asset_gpkg}_{row.asset_layer}_adaptation_costs_avoided_EAD_EAEL.csv",
        )
    return paths

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


rule benefit_cost_ratio_coastal_protection_aggregate_adaptation:
    """
    Generate the aggregate benefit cost ratio for each coastal defence
    adaptation option. That is, what is the BCR for each coastal defence zone
    given:
        1) Coastal protection build costs
        2) Avoided damages and losses to all protected assets

    Test with:
    snakemake -c1 results/adaptation_benefits_costs_bcr/coastal_coastal_protection_feature_areas_adaptation_costs_avoided_EAD_EAEL.csv
    """
    input:
        script = "workflow/5b_coastal_adaptation/6a_aggregate_adaptation.py",
        costs = generate_coastal_protection_cost_paths,
        EAD_EAEL = generate_coastal_protection_avoided_EAD_EAEL_paths,
        maps = generate_coastal_protection_mapping_paths,
        network = config["paths"]["network_layers"],
    output:
        BCR = f"{OUTPUT}/adaptation_benefits_costs_bcr/coastal_coastal_protection_feature_areas_adaptation_benefits_costs_bcr.csv",
        EAD_EAEL = f"{OUTPUT}/adaptation_benefits_costs_bcr/coastal_coastal_protection_feature_areas_adaptation_costs_avoided_EAD_EAEL.csv",
    shell:
        """
        COST=""
        for VALUE in {input.costs}; do
            COST="$COST --cost $VALUE"
        done

        EAD_EAEL=""
        for VALUE in {input.EAD_EAEL}; do
            EAD_EAEL="$EAD_EAEL --EAD-EAEL $VALUE"
        done

        LOOKUPS=""
        for VALUE in {input.maps}; do
            LOOKUPS="$LOOKUPS --lookup $VALUE"
        done

        python {input.script} \\
            $COST \\
            $EAD_EAEL \\
            $LOOKUPS \\
            --network {input.network} \\
            --output-BCR {output.BCR} \\
            --output-EAD-EAEL {output.EAD_EAEL}
        """


def generate_coastal_protection_avoided_damage_paths(wildcards) -> list[str]:
    metadata = pd.read_csv(config["paths"]["network_layers"])
    paths = []
    for row in metadata.itertuples():
        paths.append(
            f"{OUTPUT}/direct_damages_summary/{row.asset_gpkg}_{row.asset_layer}_damages.parquet",
        )
    return paths

def generate_coastal_protection_avoided_loss_paths(wildcards) -> list[str]:
    metadata = pd.read_csv(config["paths"]["network_layers"])
    paths = []
    for row in metadata.itertuples():
        paths.append(
            f"{OUTPUT}/direct_damages_summary/{row.asset_gpkg}_{row.asset_layer}_losses.parquet",
        )
    return paths

rule benefit_cost_ratio_coastal_protection_aggregate_rp:
    """
    Generate the aggregate avoided return period damages and losses for each
    coastal defence adaptation option.

    Test with:
    snakemake -c1 results/direct_damages_summary/coastal_protection_feature_areas_damages.parquet
    """
    input:
        script = "workflow/5b_coastal_adaptation/6b_aggregate_rp.py",
        damage = generate_coastal_protection_avoided_damage_paths,
        loss = generate_coastal_protection_avoided_loss_paths,
        maps = generate_coastal_protection_mapping_paths,
        network = config["paths"]["network_layers"],
    output:
        exposure = f"{OUTPUT}/direct_damages_summary/coastal_protection_feature_areas_exposures.parquet",
        damage = f"{OUTPUT}/direct_damages_summary/coastal_protection_feature_areas_damages.parquet",
        loss = f"{OUTPUT}/direct_damages_summary/coastal_protection_feature_areas_losses.parquet",
    shell:
        """
        DAMAGE=""
        for VALUE in {input.damage}; do
            DAMAGE="$DAMAGE --damage $VALUE"
        done

        LOSS=""
        for VALUE in {input.loss}; do
            LOSS="$LOSS --loss $VALUE"
        done

        LOOKUPS=""
        for VALUE in {input.maps}; do
            LOOKUPS="$LOOKUPS --lookup $VALUE"
        done

        python {input.script} \\
            $DAMAGE \\
            $LOSS \\
            $LOOKUPS \\
            --network {input.network} \\
            --output-exposure {output.exposure} \\
            --output-damage {output.damage} \\
            --output-loss {output.loss}
        """
