"""
Generate the files required by irv-jamaica/etl/damage_*_files.csv

REQUIRES:
- f"{DATA}/networks/network_layers_hazard_intersections_details.csv" to be present and populated

"""

import pandas

rule add_uids:
    """
    Add UIDs to a file.
    
    TODO: identify script that does it - think it's a Jupyter notebook
    
    Test with:
    snakemake -c1 results/direct_damages_summary_uids/roads_edges_losses.parquet
    """
    wildcard_constraints:
        dir = r"[\w_]+",
    input:
        f"{OUTPUT}/{{dir}}/{{path}}"
    output:
        f"{OUTPUT}/{{dir}}_uids/{{path}}"
    shell:
        """
        touch {output}
        """

rule RENAME_EAD_EAEL_FILE:
    """
    There's a potential mismatch between 
    {OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv created in loss_summary and
    {OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.parquet required as output targets
    
    This rule renames the former to the latter.
    """
    input:
        f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv",
    output:
        f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.parquet",
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Renaming {input} to {output}"
            cp {input} {output}
        fi
        """

rule loss_summary:
    """
    Summarise all the loss files for an asset.
    
    scripts/analysis/damage_loss_summarised.py
    
    Test with:
    snakemake -c1 results/direct_damages_summary/roads_edges_losses.parquet
    """
    input:
        network_csv = config["paths"]["network_layers"],
        parameter_combinations = f"{DATA}/parameter_combinations.csv",
        direct_damage_results = expand(
            f"{OUTPUT}/direct_damages/{{{{gpkg}}}}_{{{{layer}}}}/{{{{gpkg}}}}_{{{{layer}}}}_direct_damages_parameter_set_{{parameters}}.parquet",
            parameters=range(13),
        ),
        EAD_EAEL_damage_results = expand(
            f"{OUTPUT}/direct_damages/{{{{gpkg}}}}_{{{{layer}}}}/{{{{gpkg}}}}_{{{{layer}}}}_EAD_EAEL_parameter_set_{{parameters}}.csv",
            parameters=range(13),
        ),
    output:
        f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_losses.parquet",
        f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_exposures.parquet",
        f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_damages.parquet",
        f"{OUTPUT}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv",
    shell:
        """
        touch {output}
        """

def get_asset_row(wildcards) -> pandas.Series:
    """
    Get the path of an asset by its gpkg and layer strings.
    """
    df = pandas.read_csv(f"{DATA}/networks/network_layers_hazard_intersections_details.csv")
    row = df[(df['asset_gpkg'] == wildcards.gpkg) & (df['asset_layer'] == wildcards.layer)]
    if len(row) == 0:
        # Try the other files, and raise a warning if the value is in them
        for f in [
            f"{DATA}/networks/network_layers_hazard_intersections_details_0.csv",
            f"{DATA}/networks/network_layers_hazard_intersections_details_1.csv",
        ]:
            df = pandas.read_csv(f)
            row = df[(df['asset_gpkg'] == wildcards.gpkg) & (df['asset_layer'] == wildcards.layer)]
            if len(row) > 0:
                print(f"WARNING: Found asset in {f}")
                break
        else:
            raise ValueError(f"Multiple assets found for gpkg={wildcards.gpkg} and layer={wildcards.layer}")
    return row.squeeze()

hazard_types = [
    "TC",
    "flooding",
]

rule direct_damage_results:
    """
    Calculate direct damages for an asset across all hazards with a given parameter set.
    
    scripts/analysis/damage_calculations.py
    
    This script is called by scripts/analysis/damage_loss_setup_script.py which assigns different input args to it for each run.
    
    TODO: Update the script to accept parameter set id rather than parameter set values directly.
    
    Test with:
    snakemake -c1 results/direct_damages/roads_edges_direct_damages_parameter_set_0.parquet
    """
    input:
        # script args
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        hazard_csv = config["paths"]["hazard_layers"],
        damage_curves_csv = f"{DATA}/damage_curves/asset_damage_curve_mapping.csv",
        hazard_damage_parameters_csv = f"{DATA}/damage_curves/hazard_damage_parameters.csv",
        parameter_combinations = f"{DATA}/parameter_combinations.csv",

        # internal reads
        gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
        # damage_curves is required in create_damage_curves
        damage_curves = lambda wildcards: expand(
            f"{DATA}/damage_curves/damage_curves_{get_asset_row(wildcards).sector}_{{hazard_type}}.xlsx",
            hazard_type = hazard_types
        ),
        hazard_intersection_file = f"{OUTPUT}/hazard_asset_intersection/{{gpkg}}_splits__hazard_layers__{{layer}}.geoparquet",
    output:
        f"{OUTPUT}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages_parameter_set_{{parameter_set}}.parquet",
    shell:
        """
        touch {output}
        """

rule extract_network_layer:
    """
    Split networks into nodes, edges, and areas.
    
    scripts/exposure/split_networks.py
    
    TODO: The script must be split to support layer as a wildcard
    TODO: Check the functions below the script to see if they require other files
    
    Test with:
    snakemake -c1 results/hazard_asset_intersection/roads_splits__hazard_layers__edges.geoparquet
    """
    input:
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        hazard_csv = config["paths"]["hazard_layers"],

        gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
    output:
        # hazard_with_transforms = config["paths"]["hazard_layers"].replace(".csv", "__with_transforms.csv"),
        geoparquet = f"{OUTPUT}/hazard_asset_intersection/{{gpkg}}_splits__hazard_layers__{{layer}}.geoparquet",
    shell:
        """
        touch {output[geoparquet]}
        """

def get_single_failure_scenarios(wildcards):
    sfs = get_asset_row(wildcards).single_failure_scenarios
    if sfs == "None" or sfs == "none":
        return []
    return f"{OUTPUT}/{sfs}"

rule EAD_EAEL_results:
    """
    Calculate Estimated Annual Damages for an asset across all hazards with a given parameter set.
    
    scripts/analysis/ead_eael_calculations.py
    
    This script is called by scripts/analysis/flood_changes_setup.py which assigns different input args to it for each run.
    
    TODO: Update the script to accept parameter set id rather than parameter set values directly.
    TODO: Script should default flood_protection_name to None
    
    Test with:
    snakemake -c1 results/direct_damages/roads_edges_EAD_EAEL_parameter_set_0.csv
    """
    input:
        # script args
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        hazard_csv = config["paths"]["hazard_layers"],
        parameter_combinations = f"{DATA}/parameter_combinations.csv",

        # internal reads
        gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
        damage_file = f"{OUTPUT}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages_parameter_set_{{parameter_set}}.parquet",
        # damage_curves is required in create_damage_curves
        single_failure_scenarios = get_single_failure_scenarios,
    output:
        f"{OUTPUT}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_EAD_EAEL_parameter_set_{{parameter_set}}.csv",
    shell:
        """
        touch {output}
        """

rule parameter_combinations:
    """
    Generate parameter combinations for the damage calculations.
    
    The original parameter_combinations.txt file is created by half a dozen or so scripts.
    They all seem to set the same combinations.
    See scripts/analysis/damage_scripts_setup.py for an example.
    
    Here we'll create one canonical version called parameter_combinations.csv.
    
    Test with:
    snakemake -c1 results/direct_damages_summary/parameter_combinations.txt
    """
    output:
        f"{DATA}/parameter_combinations.csv"
    shell:
        """
        touch {output}
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
            f"{OUTPUT}/transport_failures/scenario_results/single_link_failure_{{chunk}}.csv",
            chunk=range(config["single_link_failure_chunk_count"]),
        ),
        labour_flows = f"{OUTPUT}/flow_mapping/labour_trips_and_activity.pq",
        bridges = f"{DATA}/networks/transport/roads.gpkg",  # bridges are a layer in the roads network
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        bridge_labour_trips = f"{OUTPUT}/flow_mapping/origins_destinations_labour_economic_activity.csv",
        od_losses = f"{OUTPUT}/flow_mapping/origins_destinations_trade_economic_activity.csv",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
    output:
        [
            f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_road_rail_edges_economic_losses.csv",
            f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_road_bridges_economic_losses.csv",
            f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_ports_economic_losses.csv",
            f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_rail_stations_economic_losses.csv",
            f"{OUTPUT}/economic_losses/single_failure_scenarios/single_point_failure_airports_economic_losses.csv",
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
        edge_split_map = f"{OUTPUT}/transport_failures/transport_scenario_edge_map.csv",
        edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",

        read_flow_data = [
            f"{OUTPUT}/transport_failures/nominal/labour/network.gpq",
            f"{OUTPUT}/transport_failures/nominal/trade/network.gpq",
            f"{OUTPUT}/transport_failures/nominal/labour/flows.pq",
            f"{OUTPUT}/transport_failures/nominal/trade/flows.pq",
            f"{OUTPUT}/transport_failures/nominal/labour/edges.pq",
            f"{OUTPUT}/transport_failures/nominal/trade/edges.pq",
            f"{OUTPUT}/transport_failures/nominal/all_flows.pq",
            f"{OUTPUT}/transport_failures/nominal/trade/trade_sectors.json",
        ]
    output:
        single_link_failures = expand(
            f"{OUTPUT}/transport_failures/scenario_results/single_link_failure_{{chunk}}.csv",
            chunk=range(config["single_link_failure_chunk_count"]),
        ),
    shell:
        """
        for f in {output.single_link_failures}; do
            touch $f
        done
        """

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
        edge_split_map = "{OUTPUT}/transport_failures/transport_scenario_edge_map.csv",
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

rule collate_flow_data:
    """
    Collate flow data for transport failure analysis.
    
    scripts/transport_model/transport_failure_scenario_setup.py:collate_data_flow
    
    Test with:
    snakemake -c1 results/transport_failures/nominal/all_flows.pq
    """
    input:
        labour_flow_edges = f"{DATA}/networks/transport/multi_modal_network.gpkg",
        trade_flow_edges = f"{OUTPUT}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
        trade_flows = f"{OUTPUT}/flow_mapping/sector_to_ports_flow_paths.pq",
        labour_flows = f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
    output:
        [
            f"{OUTPUT}/transport_failures/nominal/labour/network.gpq",
            f"{OUTPUT}/transport_failures/nominal/trade/network.gpq",
            f"{OUTPUT}/transport_failures/nominal/labour/flows.pq",
            f"{OUTPUT}/transport_failures/nominal/trade/flows.pq",
            f"{OUTPUT}/transport_failures/nominal/labour/edges.pq",
            f"{OUTPUT}/transport_failures/nominal/trade/edges.pq",
            f"{OUTPUT}/transport_failures/nominal/all_flows.pq",
            f"{OUTPUT}/transport_failures/nominal/trade/trade_sectors.json",
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
        exports = f"{OUTPUT}/macroeconomic_data/domestic_export_by_sector.xlsx",
        imports = f"{OUTPUT}/macroeconomic_data/import_by_industry.xlsx",
        # fuel_shares reads same file
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
        # Files read in the trade_details loop
        trade_agriculture = f"{DATA}/agriculture_data/agriculture_gdp.gpkg",
        trade_mining = f"{DATA}/mining_data/mining_gdp.gpkg",
    output:
        [
            f"{OUTPUT}/flow_mapping/sector_to_ports_flow_paths.csv",
            f"{OUTPUT}/flow_mapping/sector_to_ports_flow_paths.pq",
            f"{OUTPUT}/flow_mapping/sector_imports_exports_to_ports_flows.gpkg",
            f"{OUTPUT}/flow_mapping/origins_destinations_trade_economic_activity.csv"
        ]
    shell:
        """
        for f in {output}; do
            touch $f
        done
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
            f"{OUTPUT}/flow_mapping/road_nodes_labour_economic_activity_aggregations.gpkg",
            f"{OUTPUT}/flow_mapping/labour_to_sectors_flow_paths.csv",
            f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.csv",
            f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
            f"{OUTPUT}/flow_mapping/origins_destinations_labour_economic_activity.csv",
        ]
    shell:
        """
        for f in {output}; do
            touch $f
        done
        """

rule RENAME_LABOUR_FILE:
    """
    There's a potential mismatch between 
    {OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq created in labour_to_work_flow_mapping and
    {OUTPUT}/flow_mapping/labour_trips_and_activity.pq required by single_point_failure_road_rail
    
    This rule renames the former to the latter.
    """
    input:
        f"{OUTPUT}/flow_mapping/labour_to_sectors_trips_and_activity.pq",
    output:
        f"{OUTPUT}/flow_mapping/labour_trips_and_activity.pq",
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Renaming {input} to {output}"
            cp {input} {output}
        fi
        """

rule FAKE_MACROECONOMIC_FILES:
    """
    Create fake macroeconomic files.
    
    These are placeholders for the real files.
    
    Test with:
    snakemake -c1 results/macroecomic_data/domestic_export_by_sector.xlsx
    """
    output:
        [
            f"{OUTPUT}/macroeconomic_data/domestic_export_by_sector.xlsx",
            f"{OUTPUT}/macroeconomic_data/import_by_industry.xlsx",
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

rule FAKE_ELECTRICTY_SINGLE_POINT_FAILURES:
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
        f"{OUTPUT}/buildings/buildings_assigned_economic_activity.gpkg"
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Faking buildings economic activity file {output}"
            touch {output}
        fi
        """

rule FAKE_BUILD_RCP_EPOCH:
    """
    Several target files are in the pattern
    results/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_{hazard}__rcp_{rcp}__epoch_{epoch}__{output}.parquet
    
    Nowhere generates the underlying results/direct_damages_summary/buildings_assigned_economic_activity_areas_{hazard}__rcp_{rcp}__epoch_{epoch}__{dimension}.parquet
        files.
    
    scripts/preprocess/hazard_metadata.py generates similar files, but not the ones required here, and outputs them to processed_data/ not results/
    """
    output:
        f"{OUTPUT}/direct_damages_summary_uids/buildings_assigned_economic_activity_areas_{{hazard}}__rcp_{{rcp}}__epoch_{{epoch}}__{{dimension}}.parquet"
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Faking buildings economic activity file {output}"
            touch {output}
        fi
        """
