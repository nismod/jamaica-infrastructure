rule mine_flow_mapping:
    """
    Identify mineral flows from mines across transport network
    Test with:
    snakemake -c1 processed_data/flow_mapping/mines_flows.gpkg
    """
    input:
        land_use = "{PROCESSED_DATA}/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg",
        ports = "{PROCESSED_DATA}/networks/transport/port_polygon.gpkg",
        multi_modal_network = "{PROCESSED_DATA}/networks/transport/multi_modal_network.gpkg",
    output:
        mining_areas = "{PROCESSED_DATA}/mining_data/mining_gdp.gpkg",
        mining_flow_paths = "{PROCESSED_DATA}/flow_mapping/mines_flow_paths.csv",
        mining_flows = "{PROCESSED_DATA}/flow_mapping/mines_flows.gpkg",
    script:
        "./mines_flow_mapping.py"
