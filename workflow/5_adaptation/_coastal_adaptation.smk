"""
Generate the required files for the coastal adaptation section
"""

from typing import List
import pandas

rule coastal_flood_protection_assets:
    """
        Generates the coastal protection areas and coastal protection segments

        Test with:
        snakemake -c1 coastal_flood_protection_assets
    """
    input:
        script = "workflow/5_adaptation/generate_coastal_adaptations.py",
        island_inputs = f"{DATA}/adaptation/coastal_adaptation/Foundation_Data.gpkg"
    params:
        inland_buffer = config["coastal_adaptation"]["inland_buffer"],
        flood_threshold = config["coastal_adaptation"]["flood_threshold"],
        eps = config["coastal_adaptation"]["eps"],
        minpts = config["coastal_adaptation"]["minpts"],
        coast_length = config["coastal_adaptation"]["coast_length"],
        overlap = config["coastal_adaptation"]["overlap"],
    output:
        cpa_processing = f"{OUTPUT}/coastal_protection_assets/coastal_protection_processing_layers.gpkg",
        cpa = f"{OUTPUT}/coastal_protection_assets/Jamaica_coastal_protection_areas.gpkg"
    shell:
        """
        python {input.script} \
            --island-inputs {input.island_inputs} \
            --output-dir {OUTPUT} \
            --data-dir {DATA} \
            --inland-buffer {params.inland_buffer} \
            --flood-threshold {params.flood_threshold} \
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
    output:
        combined_cpa = f"{OUTPUT}/coastal_protection_assets/combined_coastal_protection_area.gpkg"
    shell:
        """
        python {input.script} \
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \
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
    output:
        flood_asset_dict = f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/{{gpkg}}_{{layer}}_coastal_filtered.parquet"
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --processed-data-path {input.processed_data_path} \
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \
            --combine-coastal-protection {input.combine_coastal_protection} \
            --asset-gpkg {wildcards.gpkg} \
            --asset-layer {wildcards.layer} \
            --output-dir {OUTPUT} \
        """

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
    output:
        flood_asset_count = f"{OUTPUT}/coastal_protection_assets/coastal_protection_assets_breakdown.csv"
    shell:
        """
        python {input.script} \
            --network-csv {input.network_csv} \
            --processed-data-path {input.processed_data_path} \
            --coastal-adaptation-assets {input.coastal_adaptation_assets} \
            --output-dir {OUTPUT} \
        """


