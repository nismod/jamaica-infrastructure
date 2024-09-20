rule labour_to_work_flow_mapping:
    """
    Identify commuter flows from settlements to sites of work across transport network
    Test with:
    snakemake -c1 processed_data/flow_mapping/mines_flows.gpkg
    """
    input:
        multi_modal_network = "{PROCESSED_DATA}/networks/transport/multi_modal_network.gpkg",
        buildings = "{PROCESSED_DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        population = "{PROCESSED_DATA}/population/population_projections.gpkg",
    output:
        economic_activity = "{PROCESSED_DATA}/flow_mapping/road_nodes_labour_economic_activity_aggregations.gpkg",
        labour_to_sector_activity = "{PROCESSED_DATA}/flow_mapping/labour_to_sectors_trips_and_activity.csv",
        node_activity = "{PROCESSED_DATA}/flow_mapping/origins_destinations_labour_economic_activity.csv",
        labour_flows = "{PROCESSED_DATA}/flow_mapping/labour_to_sectors_flow_paths.csv",
    script:
        "./labour_to_work_flow_mapping.py"
