rule transport_scenario_edge_map:
    """
    We split the transport maps into chunks for parallel processing.
    
    This script determines the min/max edge numbers at the boundary of each split.
    
    TODO: turn into a standalone script
    
    Test with:
    snakemake -c1 results/transport_failures/transport_scenario_edge_map.csv
    """
    input:
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
    output:
        edge_split_map = "{output_path}/transport_failures/transport_scenario_edge_map.csv",
    run:
        # Adapted from scripts/transport_model/transport_failure_scenario_setup.py
        import geopandas
        import numpy

        edges = geopandas.read_file(input.edges, layer="edges")
        rail_edges = edges[(edges["from_mode"] == "rail") & (edges["to_mode"] == "rail")][
            "edge_id"
        ].values.tolist()
        road_edges = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")][
            "edge_id"
        ].values.tolist()
        edge_fail = rail_edges + road_edges

        num_values = numpy.linspace(0, len(edge_fail) - 1, config["single_link_failure_chunk_count"])
        with open(output.edge_split_map, "w+") as f:
            f.write("id,minEdge,maxEdge\n")
            for n in range(len(num_values) - 1):
                f.write(f"{n},{int(num_values[n])},{int(num_values[n + 1])}\n")


rule FAKE_MACROECONOMIC_FILES:
    """
    Create fake macroeconomic files.
    
    These are placeholders for the real files.
    
    Test with:
    snakemake -c1 results/macroecomic_data/domestic_export_by_sector.xlsx
    """
    output:
        [
            "{output_path}/macroeconomic_data/domestic_export_by_sector.xlsx",
            "{output_path}/macroeconomic_data/import_by_industry.xlsx",
        ]
    shell:
        """
        for f in {output}; do
            if [ ! -s "$f" ]; then
                echo "WARNING: Faking macroeconomic file $f"
                touch $f
            fi
        done
        """


rule collate_flow_data:
    """
    Collate flow data for transport failure analysis.
    
    scripts/transport_model/transport_failure_scenario_setup.py:collate_data_flow
    
    Test with:
    snakemake -c1 results/transport_failures/nominal/all_flows.pq
    """
    input:
        labour_flow_edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        trade_flow_edges = "{output_path}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
        trade_flows = "{output_path}/flow_mapping/sector_to_ports_flow_paths.pq",
        labour_flows = "{output_path}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
    output:
        [
            "{output_path}/transport_failures/nominal/labour/network.gpq",
            "{output_path}/transport_failures/nominal/trade/network.gpq",
            "{output_path}/transport_failures/nominal/labour/flows.pq",
            "{output_path}/transport_failures/nominal/trade/flows.pq",
            "{output_path}/transport_failures/nominal/labour/edges.pq",
            "{output_path}/transport_failures/nominal/trade/edges.pq",
            "{output_path}/transport_failures/nominal/all_flows.pq",
            "{output_path}/transport_failures/nominal/trade/trade_sectors.json",
        ]
    shell:
        """
        for f in {input} {output}; do
            touch $f
        done
        """


rule trade_activity_flow_mapping:
    """
    Create a mapping of trade activity to flows.
    
    scripts/transport_model/trade_activity_flow_mapping.py
    
    Test with:
    snakemake -c1 results/flow_mapping/sector_imports_exports_to_ports_flows.gpkg
    """
    input:
        jam_ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        nodes = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        # edges = same file as nodes
        exports = "{output_path}/macroeconomic_data/domestic_export_by_sector.xlsx",
        imports = "{output_path}/macroeconomic_data/import_by_industry.xlsx",
        # fuel_shares reads same file
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        # Files read in the trade_details loop
        trade_agriculture = f"{DATA}/agriculture_data/agriculture_gdp.gpkg",
        trade_mining = f"{DATA}/mining_data/mining_gdp.gpkg",
    output:
        [
            "{output_path}/flow_mapping/sector_to_ports_flow_paths.csv",
            "{output_path}/flow_mapping/sector_to_ports_flow_paths.pq",
            "{output_path}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
            "{output_path}/flow_mapping/origins_destinations_trade_economic_activity.csv"
        ]
    shell:
        """
        python scripts/smk-analysis/trade_activity_flow_mapping.py \
            --jam_ports {input.jam_ports} \
            --nodes {input.nodes} \
            --exports {input.exports} \
            --imports {input.imports} \
            --buildings {input.buildings} \
            --trade_agriculture {input.trade_agriculture} \
            --trade_mining {input.trade_mining} \
            --output_path {wildcards.output_path}/flow_mapping
        """


rule RENAME_LABOUR_FILE:
    """
    There's a potential mismatch between 
    {{output_path}}/flow_mapping/labour_to_sectors_trips_and_activity.pq created in labour_to_work_flow_mapping and
    {{output_path}}/flow_mapping/labour_trips_and_activity.pq required by single_point_failure_road_rail
    
    This rule renames the former to the latter.
    """
    input:
        "{output_path}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
    output:
        "{output_path}/flow_mapping/labour_trips_and_activity.pq",
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Renaming {input} to {output}"
            cp {input} {output}
        fi
        """


rule labour_to_work_flow_mapping:
    """
    Create a mapping of labour to work flows.
    
    scripts/transport_model/labour_to_work_flow_mapping.py
    
    Test with:
    snakemake -c1 results/flow_mapping/labour_to_sectors_trips_and_activity.pq
    """
    input:
        nodes = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        population = f"{DATA}/population/population_projections.gpkg",
    output:
        [
            "{output_path}/flow_mapping/road_nodes_labour_economic_activity_aggregations.gpkg",
            "{output_path}/flow_mapping/labour_to_sectors_flow_paths.csv",
            "{output_path}/flow_mapping/labour_to_sectors_trips_and_activity.csv",
            "{output_path}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
            "{output_path}/flow_mapping/origins_destinations_labour_economic_activity.csv",
        ]
    shell:
        """
        python scripts/smk-analysis/labour_to_work_flow_mapping.py \
            --nodes {input.nodes} \
            --buildings {input.buildings} \
            --population {input.population} \
            --output_path {wildcards.output_path}/flow_mapping
        """


rule FAKE_ELECTRICTY_SINGLE_POINT_FAILURES:
    """
    The electricity single point failures are generated elsewhere.
    """
    output:
        [
            "{output_path}/electricity_failures/single_point_failure_results_nodes.csv",
            "{output_path}/electricity_failures/single_point_failure_results_edges.csv",
        ]
    shell:
        """
        for f in {output}; do
            if [ ! -s "$f" ]; then
                echo "WARNING: Faking electricity single point failure file $f"
                touch $f
            fi
        done
        """


rule single_point_failure_electricity_water:
    """
    Create a single point failure file for electricity and water assets.
    
    This is a placeholder for the real file.
    
    Test with:
    snakemake -c1 results/single_point_failures/electricity_water_single_point_failures.csv
    """
    input:
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        potable_economic_activity_buildings = f"{DATA}/networks_economic_activity/potable_facilities_buildings_economic_activity_mapping.csv",
        potable_economic_activity = f"{DATA}/networks_economic_activity/potable_facilities_dependent_economic_activity.csv",
        pipelines_economic_activity = f"{DATA}/networks_economic_activity/potable_pipelines_dependent_economic_activity.csv",
        irrigation_economic_activity = f"{DATA}/networks_economic_activity/irrigation_nodes_dependent_economic_activity.csv",
        irrigation_edges_economic_activity = f"{DATA}/networks_economic_activity/irrigation_edges_dependent_economic_activity.csv",
        electricity_economic_activity = f"{DATA}/networks_economic_activity/electricity_dependent_economic_activity.csv",
        electricity_nodes_failure_results = "{output_path}/electricity_failures/single_point_failure_results_nodes.csv",
        electricity_edges_failure_results = "{output_path}/electricity_failures/single_point_failure_results_edges.csv",
        electricity_water_mapping = f"{DATA}/networks/energy/mapping_water_to_electricity.csv",
        electricity_economic_activity_buildings = f"{DATA}/networks_economic_activity/electricity_buildings_economic_activity_mapping.csv",
    output:
        potable_facilities = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_potable_facilities_economic_losses.csv",
        potable_pipelines = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_potable_pipelines_economic_losses.csv",
        irrigation_nodes = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_irrigation_nodes_economic_losses.csv",
        irrigation_edges = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_irrigation_edges_economic_losses.csv",
        electricity_nodes_no_water = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_electricity_nodes_no_water.csv",
        electricity_edges_no_water = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_electricity_edges_no_water.csv",
        electricity_nodes = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_electricity_nodes_economic_losses.csv",
        electricity_edges = "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_electricity_edges_economic_losses.csv",
    shell:
        """
        touch {output.potable_facilities}
        touch {output.potable_pipelines}
        touch {output.irrigation_nodes}
        touch {output.irrigation_edges}
        touch {output.electricity_nodes_no_water}
        touch {output.electricity_edges_no_water}
        touch {output.electricity_nodes}
        touch {output.electricity_edges}
        """


rule single_point_failure_road_rail:
    """
    Create a single point failure file for road and rail assets.
    
    scripts/analysis/transport_single_point_failure_results_combine.py
    
    TODO: scripts/transport_model/transport_failure_scenario_setup.py needs adjusting
        so its output files use a pattern of `*_#chunk.csv` rather than
        `*_#minEdge_#maxEdge.csv` so we can know the names of files from a single chunks parameter.
    
    Test with:
    snakemake -c1 results/single_point_failures/roads_edges_single_point_failures.csv
    """
    input:
        # single_link_failures are read into all_failures in the walk through scenario_results/ directory
        single_link_failures = expand(
            "{{output_path}}/transport_failures/scenario_results/single_link_failure_{chunk}.csv",
            chunk=range(config["single_link_failure_chunk_count"]),
        ),
        labour_flows = "{output_path}/flow_mapping/labour_trips_and_activity.pq",
        bridges = f"{DATA}/networks/transport/roads.gpkg",  # bridges are a layer in the roads network
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        bridge_labour_trips = "{output_path}/flow_mapping/origins_destinations_labour_economic_activity.csv",
        od_losses = "{output_path}/flow_mapping/origins_destinations_trade_economic_activity.csv",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
    output:
        [
            "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_road_rail_edges_economic_losses.csv",
            "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_road_bridges_economic_losses.csv",
            "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_ports_economic_losses.csv",
            "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_rail_stations_economic_losses.csv",
            "{output_path}/economic_losses/single_failure_scenarios/single_point_failure_airports_economic_losses.csv",
        ]
    shell:
        """
        for f in {output}; do
            touch $f
        done
        """


rule single_link_failures:
    """
    Create single link failure results.
    
    scripts/transport_model/transport_failure_analysis.py 
        
    Test with:
    snakemake -c1 results/single_point_failures/roads_edges_single_point_failures.csv
    """
    input:
        edge_split_map = "{output_path}/transport_failures/transport_scenario_edge_map.csv",
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",

        read_flow_data = [
            "{output_path}/transport_failures/nominal/labour/network.gpq",
            "{output_path}/transport_failures/nominal/trade/network.gpq",
            "{output_path}/transport_failures/nominal/labour/flows.pq",
            "{output_path}/transport_failures/nominal/trade/flows.pq",
            "{output_path}/transport_failures/nominal/labour/edges.pq",
            "{output_path}/transport_failures/nominal/trade/edges.pq",
            "{output_path}/transport_failures/nominal/all_flows.pq",
            "{output_path}/transport_failures/nominal/trade/trade_sectors.json",
        ]
    output:
        single_link_failures = expand(
            "{{output_path}}/transport_failures/scenario_results/single_link_failure_{chunk}.csv",
            chunk=range(config["single_link_failure_chunk_count"]),
        ),
    shell:
        """
        for f in {output.single_link_failures}; do
            touch $f
        done
        """
