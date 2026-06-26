"""
Prepare data for ingestion by irv-jamaica repository.
"""

BUILDING_HAZARDS = ["coastal", "cyclone", "fluvial", "surface"]
BUILDING_RCPS = ["rcp_2.6", "rcp_4.5", "rcp_8.5", "rcp_baseline"]
BUILDING_EPOCHS = ["epoch_2010", "epoch_2030", "epoch_2050", "epoch_2070", "epoch_2080", "epoch_2100"]
BUILDING_METRICS = ["exposures", "damages", "losses"]


rule preprocess_for_visualisation:
    """
    Tag asset data with unique IDs and reserialise risk data into parquet format.
    
    Test with:
    snakemake -c1 results/direct_damages_summary_uids/roads_edges_losses.parquet
    """
    input:
        script = "workflow/7_etl/label_assets.py",
        network_csv = config["paths"]["network_layers"],
        losses = f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_losses.parquet",
        damages = f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_damages.parquet",
        exposures = f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_exposures.parquet",
        EAD_EAEL = f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv",
    output:
        # network_uids/ has sub-directories for sector, then _multi-layered_ geopackage files with uids
        # this is hard to accommodate as a listed output because more than one
        # invocation of this per-layer rule will want to write to the same file
        # network_with_uids = f"{DATA}/networks_uids/{{sector}}/<geopackage_path>"?,
        id_lookups = f"{DATA}/networks_uids/id_lookups/{{gpkg}}_{{layer}}_ids.parquet",
        losses = f"{OUTPUT}/direct_damages_summary_uids/{{gpkg}}_{{layer}}_losses.parquet",
        damages = f"{OUTPUT}/direct_damages_summary_uids/{{gpkg}}_{{layer}}_damages.parquet",
        exposures = f"{OUTPUT}/direct_damages_summary_uids/{{gpkg}}_{{layer}}_exposures.parquet",
        EAD_EAEL = f"{OUTPUT}/direct_damages_summary_uids/{{gpkg}}_{{layer}}_EAD_EAEL.parquet",
        network_csv_fragment = f"{DATA}/networks_uids/network_layer_{{gpkg}}_{{layer}}.csv",
    shell:
        """
        python {{input.script}} \\
            --network-csv {{input.network_csv}} \\
            --processed-data-dir {data} \\
            --results-dir {output} \\
            --asset-gpkg {{wildcards.gpkg}} \\
            --asset-layer {{wildcards.layer}}
        """.format(data=DATA, output=OUTPUT)


rule preprocess_buildings_for_visualisation:
    """
    Buildings are a special case with a different set of output files.

    Test with:
    snakemake -c1 processed_data/networks_uids/id_lookups/buildings_assigned_economic_activity_areas_ids.parquet
    """
    input:
        script = "workflow/7_etl/label_assets.py",
        network_csv = config["paths"]["network_layers"],
        losses = f"{OUTPUT}/direct_damages_summary/buildings_assigned_economic_activity_areas_losses.parquet",
        damages = f"{OUTPUT}/direct_damages_summary/buildings_assigned_economic_activity_areas_damages.parquet",
        exposures = f"{OUTPUT}/direct_damages_summary/buildings_assigned_economic_activity_areas_exposures.parquet",
        EAD_EAEL = f"{OUTPUT}/direct_damages_summary/buildings_assigned_economic_activity_areas_EAD_EAEL.csv",
    output:
        id_lookups = f"{DATA}/networks_uids/id_lookups/buildings_assigned_economic_activity_areas_ids.parquet",
        EAD_EAEL = f"{OUTPUT}/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_EAD_EAEL.parquet",
        network_csv_fragment = f"{DATA}/networks_uids/network_layer_buildings_assigned_economic_activity_areas.csv",
        # Additionally, buildings produce losses, damages and exposures for a
        # combination of RCPS, HAZARDS, EPOCHS (see script header).
        building_outputs = expand(
            f"{OUTPUT}/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_{{hazard}}__{{rcp}}__{{epoch}}__{{metric}}.parquet",
            hazard=BUILDING_HAZARDS,
            rcp=BUILDING_RCPS,
            epoch=BUILDING_EPOCHS,
            metric=BUILDING_METRICS,
        ),
    shell:
        """
        python {{input.script}} \\
            --asset-gpkg buildings_assigned_economic_activity \\
            --asset-layer areas \\
            --network-csv {{input.network_csv}} \\
            --processed-data-dir {data} \\
            --results-dir {output} \\
            {hazards} \\
            {rcps} \\
            {epochs}
        """.format(
            data=DATA,
            output=OUTPUT,
            hazards=" ".join([f"--hazards {h}" for h in BUILDING_HAZARDS]),
            rcps=" ".join([f"--rcps {r}" for r in BUILDING_RCPS]),
            epochs=" ".join([f"--epochs {e}" for e in BUILDING_EPOCHS]),
        )


def get_coastal_defence_to_asset_map_paths(wildcards: dict) -> list[str]:
    layers: pd.DataFrame = pd.read_csv(config["paths"]["network_layers"])
    return [
        f"{OUTPUT}/coastal_protection_assets/network_protection_mappings/{row.asset_gpkg}_{row.asset_layer}_coastal_filtered.parquet"
        for row in layers.itertuples()
    ]

def get_coastal_defence_avoided_cost_paths(wildcards: dict) -> list[str]:
    layers: pd.DataFrame = pd.read_csv(config["paths"]["network_layers"])
    return [
        f"{OUTPUT}/coastal_adaptation/direct_damages_summary/{row.asset_gpkg}_{row.asset_layer}_EAD_EAEL.csv"
        for row in layers.itertuples()
    ]

rule preprocess_coastal_features_for_visualisation:
    """
    Tag coastal feature asset data with unique IDs, find sum of avoided damages
    & losses for each coastal defence feature and save data in parquet format.
    
    Test with:
    snakemake -c1 results/direct_damages_summary_uids/coastal_protection_feature_areas_EAD_EAEL.parquet
    """
    input:
        script = "workflow/7_etl/label_coastal_assets.py",
        coastal_network_csv = config["coastal_adaptation"]["coastal_layer"],
        network_csv = config["paths"]["network_layers"],
        defence_to_asset_map = get_coastal_defence_to_asset_map_paths,
        avoided_costs_by_asset_class = get_coastal_defence_avoided_cost_paths,
        exposure = f"{OUTPUT}/direct_damages_summary/coastal_protection_feature_areas_exposures.parquet",
        damage = f"{OUTPUT}/direct_damages_summary/coastal_protection_feature_areas_damages.parquet",
        loss = f"{OUTPUT}/direct_damages_summary/coastal_protection_feature_areas_losses.parquet",
    output:
        id_lookups = f"{DATA}/networks_uids/id_lookups/coastal_protection_feature_areas_ids.parquet",
        network_csv_fragment = f"{DATA}/networks_uids/network_layer_coastal_protection_feature_areas.csv",
        avoided_exposure = f"{OUTPUT}/direct_damages_summary_uids/coastal_protection_feature_areas_exposures.parquet",
        avoided_damage = f"{OUTPUT}/direct_damages_summary_uids/coastal_protection_feature_areas_damages.parquet",
        avoided_loss = f"{OUTPUT}/direct_damages_summary_uids/coastal_protection_feature_areas_losses.parquet",
        avoided_EAD_EAEL = f"{OUTPUT}/direct_damages_summary_uids/coastal_protection_feature_areas_EAD_EAEL.parquet",
    shell:
        """
        python {{input.script}} \\
            --coastal-network-csv {{input.coastal_network_csv}} \\
            --network-csv {{input.network_csv}} \\
            --processed-data-dir {data} \\
            --results-dir {output} \\
            --asset-gpkg coastal_protection_feature \\
            --asset-layer areas
        """.format(data=DATA, output=OUTPUT)


def network_layer_csv_fragments_for_all_layers(wildcards) -> list[str]:
    assets: pd.DataFrame = pd.read_csv(config["paths"]["network_layers"])
    return [f"{DATA}/networks_uids/network_layer_{asset.asset_gpkg}_{asset.asset_layer}.csv" for asset in assets.itertuples()]

rule preprocess_all_for_visualisation:
    """
    Target rule to generate visualisation ready files for all layers.
    """
    input:
        network_layers = config["paths"]["network_layers"],
        all_layers = network_layer_csv_fragments_for_all_layers,
    output:
        network_layers_for_visualisation = f"{DATA}/networks_uids/network_layers.csv"
    run:
        import pandas as pd

        fragments: list[pd.DataFrame] = []
        for fragment_path in input.all_layers:
            fragments.append(pd.read_csv(fragment_path))
        pd.concat(fragments).sort_values("base_id").to_csv(output.network_layers_for_visualisation, index=False)
