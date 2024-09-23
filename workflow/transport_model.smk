rule create_multi_modal_network:
    """
    Create a multi-modal transport network.
    Test with:
    snakemake -c1 processed_data/networks/transport/multi_modal_network.gpkg
    """
    input:
        airports = f"{DATA}/networks/transport/airport_polygon.gpkg",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        rail = f"{DATA}/networks/transport/rail.gpkg",
        road = f"{DATA}/networks/transport/roads.gpkg",
    output:
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
    script:
        "../scripts/transport_model/multi_modal_network_creation.py"


rule labour_to_work_flow_mapping:
    """
    Identify commuter flows from settlements to sites of work across transport network
    Test with:
    snakemake -c1 results/flow_mapping/labour_to_sectors_flow_paths.csv
    """
    input:
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        population = f"{DATA}/population/population_projections.gpkg",
    output:
        economic_activity = f"{DATA}/flow_mapping/road_nodes_labour_economic_activity_aggregations.gpkg",
        labour_to_sector_activity = f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.csv",
        node_activity = f"{OUTPUT}/flow_mapping/origins_destinations_labour_economic_activity.csv",
        labour_flows = f"{OUTPUT}/flow_mapping/labour_to_sectors_flow_paths.csv",
    script:
        "../scripts/transport_model/labour_to_work_flow_mapping.py"


rule mine_flow_mapping:
    """
    Identify mineral flows from mines across transport network
    Test with:
    snakemake -c1 processed_data/flow_mapping/mines_flows.gpkg
    """
    input:
        land_use = f"{DATA}/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
    output:
        mining_areas = f"{DATA}/mining_data/mining_gdp.gpkg",
        mining_flow_paths = f"{OUTPUT}/flow_mapping/mines_flow_paths.csv",
        mining_flows = f"{OUTPUT}/flow_mapping/mines_flows.gpkg",
    script:
        "../scripts/transport_model/mines_flow_mapping.py"


rule trade_activity_flow_mapping:
    """
    Identify trade flows from areas of economic activity across the transport network
    Test with:
    snakemake -c1 processed_data/flow_mapping/sector_imports_exports_to_ports_flows.gpkg
    """
    input:
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        agriculture = f"{DATA}/agriculture_data/agriculture_gdp.gpkg",
        mining_areas = f"{DATA}/mining_data/mining_gdp.gpkg",
    output:
        od = f"{OUTPUT}/flow_mapping/sector_to_ports_flow_paths.csv",
        trade_flows = f"{OUTPUT}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
        economic_activity_by_node_sector = f"{OUTPUT}/flow_mapping/origins_destinations_trade_economic_activity.csv",
    script:
        "../scripts/transport_model/trade_activity_flow_mapping.py"
