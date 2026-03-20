"""
How does a network perform while missing a given link?
"""


rule collate_flow_data:
    """
    Collate flow data for transport failure analysis.

    Test with:
    snakemake -c1 results/transport_failures/nominal/
    """
    input:
        script = "workflow/4_criticality/collate_flow_data.py",
        labour_flow_edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        trade_flow_edges = f"{OUTPUT}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
        trade_flows = f"{OUTPUT}/flow_mapping/sector_to_ports_flow_paths.pq",
        labour_flows = f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
    output:
        [
            directory(f"{OUTPUT}/transport_failures/nominal"),
            f"{OUTPUT}/transport_failures/nominal/all_flows.pq",
        ]
    shell:
        f"""
        python {{input.script}} \\
            --results-dir {OUTPUT} \\
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
        edge_split_map = f"{OUTPUT}/transport_failures/transport_scenario_edge_map.csv",
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

    Test with:
    snakemake -c1 results/transport_failures/scenario_results/single_link_failure_0.csv
    """
    input:
        script = "workflow/4_criticality/single_link_failures.py",
        edge_chunk_map_csv = f"{OUTPUT}/transport_failures/transport_scenario_edge_map.csv",
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        flow_data = f"{OUTPUT}/transport_failures/nominal/",
    params:
        # include as a param to trigger re-run on change
        chunk_count = config["single_link_failure_chunk_count"]
    output:
        chunk = protected(f"{OUTPUT}/transport_failures/scenario_results/single_link_failure_{{chunk}}.csv"),
    shell:
        f"""
        python {{input.script}} \\
            --edge-chunk-map-csv {{input.edge_chunk_map_csv}} \\
            --chunk-id {{wildcards.chunk}} \\
            --edges-file {{input.edges}} \\
            --flow-data-dir {{input.flow_data}} \\
            --output-path {{output.chunk}} \\
            --hourly-wage {config["economics"]["labour_cost_JMD_per_hour"]} \\
            --trade-effect {config["economics"]["disrupted_trade_fraction"]}
        """


rule rail_stations_failure_analysis:
    """
    Remove railway stations from multi-modal transport network and estimate the
    arising economic losses.

    Test with:
    snakemake -c1 results/transport_failures/single_station_failures_scenarios.csv
    """
    input:
        script = "workflow/4_criticality/rail_stations_failure_analysis.py",
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        rail_nodes = f"{DATA}/networks/transport/rail.gpkg",
        flow_data_dir = f"{OUTPUT}/transport_failures/nominal/",
    output:
        station_failures = f"{OUTPUT}/transport_failures/single_station_failures_scenarios.csv",
    shell:
        f"""
        python {{input.script}} \\
            --edges-file {{input.edges}} \\
            --rail-nodes-file {{input.rail_nodes}} \\
            --flow-data-dir {{input.flow_data_dir}} \\
            --output-path {{output.station_failures}} \\
            --hourly-wage {config["economics"]["labour_cost_JMD_per_hour"]} \\
            --trade-effect {config["economics"]["disrupted_trade_fraction"]}
        """


rule bridge_failure_analysis:
    """
    Remove road bridges from multi-modal transport network and estimate the
    arising economic losses.

    Test with:
    snakemake -c1 results/transport_failures/single_bridge_failures_scenarios.csv
    """
    input:
        script = "workflow/4_criticality/roads_bridges_failure_analysis.py",
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        road_nodes = f"{DATA}/networks/transport/roads.gpkg",
        flow_data_dir = f"{OUTPUT}/transport_failures/nominal/",
    output:
        bridge_failures = f"{OUTPUT}/transport_failures/single_bridge_failures_scenarios.csv",
    shell:
        f"""
        python {{input.script}} \\
            --edges-file {{input.edges}} \\
            --road-nodes-file {{input.road_nodes}} \\
            --flow-data-dir {{input.flow_data_dir}} \\
            --output-path {{output.bridge_failures}} \\
            --hourly-wage {config["economics"]["labour_cost_JMD_per_hour"]} \\
            --trade-effect {config["economics"]["disrupted_trade_fraction"]}
        """


rule single_point_failure_road_rail:
    """
    Create a single point failure file for road and rail assets.

    Test with:
    snakemake -c1 results/economic_losses/single_failure_scenarios/single_point_failure_road_rail_edges_economic_losses.csv
    """
    input:
        script = "workflow/4_criticality/transport_single_point_failure_results_combine.py",
        single_link_failures = expand(
            f"{OUTPUT}/transport_failures/scenario_results/single_link_failure_{{chunk}}.csv",
            chunk=range(config["single_link_failure_chunk_count"]),
        ),
        station_failures = f"{OUTPUT}/transport_failures/single_station_failures_scenarios.csv",
        bridge_failures = f"{OUTPUT}/transport_failures/single_bridge_failures_scenarios.csv",
        labour_flows = f"{OUTPUT}/flow_mapping/labour_trips_and_activity.gpq",
        bridges = f"{DATA}/networks/transport/roads.gpkg",
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        bridge_labour_trips = f"{OUTPUT}/flow_mapping/origins_destinations_labour_economic_activity.csv",
        od_losses = f"{OUTPUT}/flow_mapping/origins_destinations_trade_economic_activity.csv",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        airports = f"{DATA}/networks/transport/airport_polygon.gpkg",
    params:
        # Include as a param to trigger re-run on change
        chunk_count = config["single_link_failure_chunk_count"]
    output:
        road_rail_edges = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_road_rail_edges_economic_losses.csv",
        road_bridges = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_road_bridges_economic_losses.csv",
        ports = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_ports_economic_losses.csv",
        rail_stations = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_rail_stations_economic_losses.csv",
        airports = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_airports_economic_losses.csv",
    shell:
        f"""
        python {{input.script}} \\
            --results-dir {OUTPUT} \\
            --processed-data-dir {DATA}
        """


rule ELECTRICTY_SINGLE_POINT_FAILURES:
    """
    The electricity single point failures are generated elsewhere.
    """
    output:
        [
            f"{OUTPUT}/electricity_failures/single_point_failure_results_nodes.csv",
            f"{OUTPUT}/electricity_failures/single_point_failure_results_edges.csv",
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
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.geoparquet",
        potable_economic_activity_buildings = f"{DATA}/networks_economic_activity/potable_facilities_buildings_economic_activity_mapping.csv",
        potable_economic_activity = f"{DATA}/networks_economic_activity/potable_facilities_dependent_economic_activity.csv",
        pipelines_economic_activity = f"{DATA}/networks_economic_activity/potable_pipelines_dependent_economic_activity.csv",
        irrigation_economic_activity = f"{DATA}/networks_economic_activity/irrigation_nodes_dependent_economic_activity.csv",
        irrigation_edges_economic_activity = f"{DATA}/networks_economic_activity/irrigation_edges_dependent_economic_activity.csv",
        electricity_economic_activity = f"{DATA}/networks_economic_activity/electricity_dependent_economic_activity.csv",
        electricity_nodes_failure_results = f"{OUTPUT}/electricity_failures/single_point_failure_results_nodes.csv",
        electricity_edges_failure_results = f"{OUTPUT}/electricity_failures/single_point_failure_results_edges.csv",
        electricity_water_mapping = f"{DATA}/networks/energy/mapping_water_to_electricity.csv",
        electricity_economic_activity_buildings = f"{DATA}/networks_economic_activity/electricity_buildings_economic_activity_mapping.csv",
    output:
        potable_facilities = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_potable_facilities_economic_losses.csv",
        potable_pipelines = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_potable_pipelines_economic_losses.csv",
        irrigation_nodes = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_irrigation_nodes_economic_losses.csv",
        irrigation_edges = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_irrigation_edges_economic_losses.csv",
        electricity_nodes_no_water = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_electricity_nodes_no_water.csv",
        electricity_edges_no_water = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_electricity_edges_no_water.csv",
        electricity_nodes = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_electricity_nodes_economic_losses.csv",
        electricity_edges = f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_electricity_edges_economic_losses.csv",
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
