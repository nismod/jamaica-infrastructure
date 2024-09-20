rule trade_activity_flow_mapping:
    """
    Identify trade flows from areas of economic activity across the transport network
    Test with:
    snakemake -c1 processed_data/flow_mapping/sector_imports_exports_to_ports_flows.gpkg
    """
    input:
        ports = "{PROCESSED_DATA}/networks/transport/port_polygon.gpkg",
        multi_modal_network = "{PROCESSED_DATA}/networks/transport/multi_modal_network.gpkg",
        buildings = "{PROCESSED_DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        agriculture = "{PROCESSED_DATA}/agriculture_data/agriculture_gdp.gpkg",
        mining_areas = "{PROCESSED_DATA}/mining_data/mining_gdp.gpkg",
    output:
        od = "{PROCESSED_DATA}/flow_mapping/sector_to_ports_flow_paths.csv",
        trade_flows = "{PROCESSED_DATA}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
        economic_activity_by_node_sector = "{PROCESSED_DATA}/flow_mapping/origins_destinations_trade_economic_activity.csv",
    script:
        "./trade_activity_flow_mapping.py"
