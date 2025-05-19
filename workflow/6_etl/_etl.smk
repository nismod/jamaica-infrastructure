"""
Prepare data for ingestion by irv-jamaica repository.
"""

rule preprocess_for_visualisation:
    """
    Tag asset data with unique IDs and reserialise risk data into parquet format.
    
    Test with:
    snakemake -c1 results/direct_damages_summary_uids/roads_edges_losses.parquet
    """
    input:
        script = "workflow/6_etl/preprocess_for_visualisation.py",
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
        f"""
        python {{input.script}} \
            --network-csv {{input.network_csv}} \
            --processed-data-dir {DATA} \
            --results-dir {OUTPUT} \
            --asset-gpkg {{wildcards.gpkg}} \
            --asset-layer {{wildcards.layer}} \
        """


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


rule BUILD_RCP_EPOCH:
    """
    Several target files are in the pattern
    results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_{hazard}__rcp_{rcp}__epoch_{epoch}__{output_path}.parquet
    
    Nowhere generates the underlying results/direct_damages_summary/buildings_assigned_economic_activity_areas_{hazard}__rcp_{rcp}__epoch_{epoch}__{dimension}.parquet
        files.
    
    scripts/preprocess/hazard_metadata.py generates similar files, but not the ones required here, and outputs them to processed_data/ not results/
    """
    output:
        "{output_path}/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_{hazard}__rcp_{rcp}__epoch_{epoch}__{dimension}.parquet"
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Faking buildings economic activity file {output}"
            touch {output}
        fi
        """


rule BUILDINGS_ECONOMIC_ACTIVITY:
    """
    Fake a mapping of buildings to economic activity.
    
    results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_EAD_EAEL.parquet
     requires generation of results/buildings/buildings_assigned_economic_activity.gpkg,
     but it isn't generated anywhere.
    
    processed_data/buildings/buildings_assigned_economic_activity.gpkg is generated, 
        but not results/buildings/buildings_assigned_economic_activity.gpkg
    """
    output:
        "{output_path}/buildings/buildings_assigned_economic_activity.gpkg"
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Faking buildings economic activity file {output}"
            touch {output}
        fi
        """