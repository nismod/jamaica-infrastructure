rule add_uids:
    """
    Add UIDs to a file.
    
    Test with:
    snakemake -c1 results/direct_damages_summary_uids/roads_edges_losses.parquet
    """
    wildcard_constraints:
        dir = r"[\w_]+",
    input:
        "{output_path}/{dir}/{path}"
    output:
        "{output_path}/{dir}_uids/{path}"
    shell:
        """
        cp {input} {output}
        """


rule FAKE_BUILD_RCP_EPOCH:
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


rule FAKE_BUILDINGS_ECONOMIC_ACTIVITY:
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