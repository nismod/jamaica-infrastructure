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
        economic_activity = f"{OUTPUT}/flow_mapping/road_nodes_labour_economic_activity_aggregations.gpkg",
        labour_to_sector_activity = f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
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
    snakemake -c1 results/flow_mapping/sector_imports_exports_to_ports_flows.gpkg
    """
    input:
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        agriculture = f"{DATA}/agriculture_data/agriculture_gdp.gpkg",
        mining_areas = f"{DATA}/mining_data/mining_gdp.gpkg",
    output:
        od = f"{OUTPUT}/flow_mapping/sector_to_ports_flow_paths.pq",
        trade_flows = f"{OUTPUT}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
        economic_activity_by_node_sector = f"{OUTPUT}/flow_mapping/origins_destinations_trade_economic_activity.csv",
    script:
        "../scripts/transport_model/trade_activity_flow_mapping.py"


rule transport_failure_flow_allocation_setup:
    """
    Create input data and chunk problem prior to exhaustively failing transport edges.
    Test with:
    snakemake -c1 results/transport_failures/nominal
    """
    input:
        script_setup = "scripts/transport_model/transport_failure_scenario_setup.py",
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        od = f"{OUTPUT}/flow_mapping/sector_to_ports_flow_paths.pq",
        trade_flows = f"{OUTPUT}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
        labour_to_sector_activity = f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
    params:
        n_chunks = 128
    output:
        collated_input_data = directory(f"{OUTPUT}/transport_failures/nominal"),
        failure_scenarios = f"{OUTPUT}/transport_failures/parallel_transport_scenario_selection.txt",
        failure_scenarios_resampled = f"{OUTPUT}/transport_failures/parallel_transport_scenario_selection_resample.txt",
    shell:
        """
        python {input.script_setup} {params.n_chunks}
        """


rule transport_failure_flow_allocation:
    """
    Exhaustively fail transport edges, reallocate economic flows.
    Uses approximately 16GB RAM per core.
    Test with:
    snakemake -c1 results/transport_failures/scenario_results
    """
    input:
        script_core = "scripts/transport_model/transport_failure_analysis.py",
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        collated_input_data = f"{OUTPUT}/transport_failures/nominal",
        failure_scenarios = f"{OUTPUT}/transport_failures/parallel_transport_scenario_selection.txt",
    threads:
        workflow.cores
    resources:
        mem_gb = 16
    output:
        transport_failures_scenarios = directory(f"{OUTPUT}/transport_failures/scenario_results"),
    run:
        import datetime
        import subprocess

        args = [
            "parallel",
            "--halt-on-error",
            "1",
            "--lb",
            "-j",
            str(workflow.cores),
            "--colsep",
            ",",
            "-a",
            input.failure_scenarios,
            "python",
            "scripts/transport_model/transport_failure_analysis.py",
            "{}",
        ]
        print(args)
        subprocess.run(args)


rule transport_accumulate_flows_to_edges:
    """
    Accumulate nominal economic flows to transport edges.
    Test with:
    snakemake -c1 results/flow_mapping/labour_trips_and_activity.pq
    """
    input:
        script = "scripts/transport_model/accumulate_flows_to_edges.py",
        multi_modal_network = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        labour_flows = f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
    threads:
        workflow.cores
    output:
        edge_flows = f"{OUTPUT}/flow_mapping/labour_trips_and_activity.pq",
        edge_flows_geometry = f"{OUTPUT}/flow_mapping/labour_trips_and_activity.gpq",
    shell:
        """
        python \
            {input.script} \
            {workflow.cores} \
            {input.multi_modal_network} \
            {input.labour_flows} \
            {output.edge_flows} \
            {output.edge_flows_geometry}
        """


rule transport_failure_flow_allocation_combine_results:
    """
    Gather step to combine chunked failured analysis.
    Test with:
    snakemake -c1 results/transport_failures/single_point_failure_road_rail_edges_economic_losses.csv
    """
    input:
        script = "scripts/transport_model/transport_failure_results_combine.py",
        edge_flows = f"{OUTPUT}/flow_mapping/labour_trips_and_activity.pq",
        transport_failures_scenarios = f"{OUTPUT}/transport_failures/scenario_results",
    output:
        road_and_rail_losses = f"{OUTPUT}/transport_failures/single_point_failure_road_rail_edges_economic_losses.csv",
    shell:
        """
        python {input.script}
        """

