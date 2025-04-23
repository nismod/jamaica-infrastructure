"""
How does a network perform while missing a given link?
"""


rule collate_flow_data:
    """
    Collate flow data for transport failure analysis.
    
    Test with:
    snakemake -c1 results/transport_failures/nominal/all_flows.pq
    """
    input:
        script = "workflow/3_criticality/collate_flow_data.py",
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
            "{output_path}/transport_failures/nominal/labour/edge_indexes.pq",
            "{output_path}/transport_failures/nominal/trade/edge_indexes.pq",
            "{output_path}/transport_failures/nominal/all_flows.pq",
            "{output_path}/transport_failures/nominal/trade/trade_sectors.json",
        ]
    shell:
        f"""
        python {{input.script}} \
            --results-dir {{wildcards.output_path}} \
            --processed-data-dir {DATA}
        """


rule transport_scenario_edge_map:
    """
    We split the transport maps into chunks for parallel processing.
    
    This script determines the min/max edge numbers at the boundary of each split.
    
    Test with:
    snakemake -c1 results/transport_failures/transport_scenario_edge_map.csv
    """
    input:
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
    params:
        # include as a param to trigger re-run on change
        chunk_count = config["single_link_failure_chunk_count"]
    output:
        edge_split_map = temp("{output_path}/transport_failures/transport_scenario_edge_map.csv"),
    run:
        import logging

        import geopandas
        import numpy

        logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

        logging.info("Filtering edges")
        edges = geopandas.read_file(input.edges, layer="edges")
        rail_edges = edges[(edges["from_mode"] == "rail") & (edges["to_mode"] == "rail")]["edge_id"].values.tolist()
        road_edges = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")]["edge_id"].values.tolist()
        edges_to_fail = rail_edges + road_edges

        n_chunk = int(config["single_link_failure_chunk_count"])
        logging.info(f"{len(edges_to_fail)} edges to fail, splitting into {n_chunk} chunks")
        indicies = numpy.round(numpy.linspace(0, len(edges_to_fail), n_chunk + 1)).astype(int)
        logging.info(f"Using following indicies:\n{indicies}")

        logging.info("Write out start and stop indicies to disk")
        with open(output.edge_split_map, "w+") as f:
            f.write("id,minEdge,maxEdge\n")
            for n in range(len(indicies) - 1):
                f.write(f"{n},{indicies[n]},{indicies[n + 1]}\n")


rule single_link_failures:
    """
    Create single link failure results.
    
    scripts/transport_model/transport_failure_analysis.py 
        
    Test with:
    snakemake -c1 results/transport_failures/scenario_results/single_link_failure_0.csv
    """
    input:
        edge_split_map = "{output_path}/transport_failures/transport_scenario_edge_map.csv",
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        read_flow_data = [
            "{output_path}/transport_failures/nominal/labour/network.gpq",
            "{output_path}/transport_failures/nominal/trade/network.gpq",
            "{output_path}/transport_failures/nominal/labour/flows.pq",
            "{output_path}/transport_failures/nominal/trade/flows.pq",
            "{output_path}/transport_failures/nominal/labour/edge_indexes.pq",
            "{output_path}/transport_failures/nominal/trade/edge_indexes.pq",
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


rule ELECTRICTY_SINGLE_POINT_FAILURES:
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
