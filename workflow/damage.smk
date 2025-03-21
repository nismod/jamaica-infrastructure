"""
Generate the files required by irv-jamaica/etl/damage_*_files.csv

REQUIRES:
- f"{DATA}/networks/network_layers_hazard_intersections_details[_#].csv" to be present and populated
- f"{DATA}/networks/transport/multi_modal_network.gpkg" to be present and populated
"""

import pandas as pd

PARAMETER_SET_IDS = range(13)
HAZARD_TYPES = (
    "TC",
    "flooding",
)


def get_asset_row(wildcards) -> pd.Series:
    """
    Get the path of an asset by its gpkg and layer strings.
    """
    df = pd.read_csv(f"workflow/network_layers.csv")
    row = df[(df['asset_gpkg'] == wildcards.gpkg) & (df['asset_layer'] == wildcards.layer)]
    if len(row) > 1:
        raise ValueError(f"Multiple assets found for gpkg={wildcards.gpkg} and layer={wildcards.layer}")
    return row.squeeze()


rule write_hazard_transforms:
    """
    Write hazard transforms to disk alongside hazard metadata. I don't think
    this is actually used?

    Test with:
    snakemake -c1 processed_data/hazards/hazard_layers__with_transforms.csv
    """
    input:
        hazard_csv = lambda wildcards: wildcards.hazard_csv + ".csv",
        data_dir = DATA,
    output:
        hazard_transforms_csv = "{hazard_csv}__with_transforms.csv"
    run:
        from jamaica_infrastructure.transform import read_transforms

        hazards = pd.read_csv(input.hazard_csv)
        hazard_transforms, transforms = read_transforms(hazards, input.data_dir)
        hazard_transforms.to_csv(output.hazard_transforms_csv, index=False)


rule rasterise_asset_layer:
    """
    Split networks into nodes, edges, and areas.
    
    Test with:
    snakemake -c1 results/hazard_asset_intersection/roads_splits__hazard_layers__edges.geoparquet
    """
    input:
        script = "scripts/smk-analysis/split_networks.py",
        networks = config["paths"]["network_layers"],
        hazards = config["paths"]["hazard_layers"],
        gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
    output:
        splits = "{output_path}/hazard_asset_intersection/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
    shell:
        f"""
        python {{input.script}} \
            --network-csv {{input.networks}} \
            --hazard-csv {{input.hazards}} \
            --data-dir {{DATA}} \
            --asset-gpkg {{wildcards.gpkg}} \
            --asset-layer {{wildcards.layer}} \
            --output-path {{output.splits}}
        """


checkpoint sensitivity_parameters:
    """
    Generate sensitivity parameter combinations for the damage calculations.
    
    The original sensitivity_parameters.txt file is created by half a dozen or so scripts.
    They all seem to set the same combinations.
    See scripts/analysis/damage_scripts_setup.py for an example.
    
    Here we'll create one canonical version called sensitivity_parameters.csv.
    
    Test with:
    snakemake -c1 processed_data/sensitivity_parameters.csv
    """
    input:
        config = "config.json"
    output:
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv"
    run:
        import pandas as pd
        from SALib.sample import morris

        variables = ("cost_uncertainty_parameter", "damage_uncertainty_parameter")

        if config["sensitivity_analysis"] == True:
            problem = {
                "num_vars": len(variables),
                "names": variables,
                "bounds": [[0, 1.0] for var in variables],
            }
            values = morris.sample(
                problem, 4, num_levels=4, optimal_trajectories=2, local_optimization=False
            )

        elif config["sensitivity_analysis"] == False:
            values = [(1.0, 1.0)]  # one sample, no variation

        else:
            raise ValueError("config['sensitivity_analysis'] must be boolean")

        df = pd.DataFrame(values, columns=variables)
        df.index.name = "set_id"
        df.to_csv(output.sensitivity_parameters, float_format='%.3f')


def sensitivity_id_from_slug(wildcards):
    return wildcards.parameter_set.replace("parameter_set_", "")

rule direct_damage:
    """
    Calculate direct damages for an asset across all hazards with a given parameter set.
    
    Test with:
    snakemake -c1 results/direct_damages/roads_edges/roads_edges_direct_damages_parameter_set_0.parquet
    """
    input:
        script = "scripts/analysis/damage_calculations.py",
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        hazard_csv = config["paths"]["hazard_layers"],
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        asset_gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
        damage_curve_mapping = f"{DATA}/damage_curves/asset_damage_curve_mapping.csv",
        damage_threshold_uplift = f"{DATA}/damage_curves/hazard_damage_parameters.csv",
        damage_curves_dir = f"{DATA}/damage_curves",
        damage_curves = lambda wildcards: expand(
            f"{DATA}/damage_curves/damage_curves_{get_asset_row(wildcards).sector}_{{hazard_type}}.xlsx",
            hazard_type = HAZARD_TYPES
        ),
        hazard_intersection_file = "{output_path}/hazard_asset_intersection/{gpkg}_splits__hazard_layers__{layer}.geoparquet",
    params:
        sensitivity_id = sensitivity_id_from_slug,
    output:
        damages = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_direct_damages_{parameter_set}.parquet",
    shell:
        f"""
        python {{input.script}} \
            --network-csv {{input.network_csv}} \
            --hazard-csv {{input.hazard_csv}} \
            --sensitivity-csv {{input.sensitivity_parameters}} \
            --sensitivity-id {{params.sensitivity_id}} \
            --asset-gpkg-file {{input.asset_gpkg}} \
            --asset-gpkg-label {{wildcards.gpkg}} \
            --asset-layer {{wildcards.layer}} \
            --damage-curve-mapping-csv {{input.damage_curve_mapping}} \
            --damage-threshold-uplift-csv {{input.damage_threshold_uplift}} \
            --damage-curves-dir {{input.damage_curves_dir}} \
            --intersection {{input.hazard_intersection_file}} \
            --output-path {{output.damages}}
        """


def damage_ensemble_files(wildcards):
    # wait for ensemble parameters to be generated, then read the length of the list
    filepath = checkpoints.sensitivity_parameters.get(**wildcards).output.sensitivity_parameters
    n_ensemble = len(pd.read_csv(filepath))
    return expand(
        "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages_parameter_set_{parameter_set}.parquet",
        parameter_set=range(n_ensemble)
    )

rule collapse_sensitivity:
    """
    Summarise direct damage results (aggregate over sensitivity analysis)

    Test with:
    snakemake -c1 results/direct_damages_summary/roads_edges_damages.csv
    """
    input:
        script = "scripts/analysis/direct_damage_summarise.py",
        network_csv = f"{DATA}/networks/network_layers_hazard_intersections_details.csv",
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        damages = damage_ensemble_files
        # EAD_and_EAEL_files
    output:
        damages = "{output_path}/direct_damages_summary/{gpkg}_{layer}_damages.csv",
        exposure = "{output_path}/direct_damages_summary/{gpkg}_{layer}_exposure.csv",
    shell:
        f"""
        # build list of input damage files: --damage <damage_file>
        DAMAGE_FILES=""
        for FILE in {{input.damages}}; do
            DAMAGE_FILES="$DAMAGE_FILES --damages $FILE"
        done

        python {{input.script}} \
            --network-csv {{input.network_csv}} \
            --sensitivity-csv {{input.sensitivity_parameters}} \
            --asset-gpkg {{wildcards.gpkg}} \
            --asset-layer {{wildcards.layer}} \
            $DAMAGE_FILES \
            --output-damages-path {{output.damages}} \
            --output-exposure-path {{output.exposure}}
        """


def get_single_failure_scenarios(wildcards):
    row = get_asset_row(wildcards)
    sfs = row.single_failure_scenarios
    if sfs == "None" or sfs == "none":
        return []
    if row.sector == "buildings":
        return f"{DATA}/{sfs}"
    return f"{wildcards.output_path}/{sfs}"


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


rule loss_summary:
    """
    Summarise all the loss files for an asset.
    
    scripts/analysis/damage_loss_summarised.py
    
    Test with:
    snakemake -c1 results/direct_damages_summary/roads_edges_losses.parquet
    """
    input:
        network_csv = config["paths"]["network_layers"],
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",
        direct_damage_results = expand(
            "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_direct_damages_parameter_set_{parameters}.parquet",
            parameters=PARAMETER_SET_IDS,
        ),
        EAD_EAEL_damage_results = expand(
            "{{output_path}}/direct_damages/{{gpkg}}_{{layer}}/{{gpkg}}_{{layer}}_EAD_EAEL_parameter_set_{parameters}.csv",
            parameters=PARAMETER_SET_IDS,
        ),
        single_failure_scenarios = get_single_failure_scenarios,
    output:
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_losses.parquet",
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_exposures.parquet",
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_damages.parquet",
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_EAD_EAEL.csv",
    shell:
        """
        python scripts/smk-analysis/damage_loss_summarised.py \
            --network_csv {input.network_csv} \
            --sensitivity_parameters {input.sensitivity_parameters} \
            --direct_damage_results {input.direct_damage_results} \
            --EAD_EAEL_damage_results {input.EAD_EAEL_damage_results} \
            --single_failure_scenarios {input.single_failure_scenarios} \
            --output_path {wildcards.output_path}/direct_damages_summary \
            --gpkg {wildcards.gpkg} \
            --layer {wildcards.layer}
        """


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
        sensitivity_parameters = f"{DATA}/sensitivity_parameters.csv",

        # internal reads
        gpkg = lambda wildcards: f"{DATA}/{get_asset_row(wildcards).path}",
        damage_file = "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_direct_damages_parameter_set_{parameter_set}.parquet",
        # damage_curves is required in create_damage_curves
        single_failure_scenarios = get_single_failure_scenarios,
    output:
        "{output_path}/direct_damages/{gpkg}_{layer}/{gpkg}_{layer}_EAD_EAEL_parameter_set_{parameter_set}.csv",
    shell:
        """
        touch {output}
        """


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


rule RENAME_EAD_EAEL_FILE:
    """
    There's a potential mismatch between 
    {{output_path}}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.csv created in loss_summary and
    {{output_path}}/direct_damages_summary/{{gpkg}}_{{layer}}_EAD_EAEL.parquet required as output targets
    
    This rule renames the former to the latter.
    """
    input:
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_EAD_EAEL.csv",
    output:
        "{output_path}/direct_damages_summary/{gpkg}_{layer}_EAD_EAEL.parquet",
    shell:
        """
        if [ ! -s "{output}" ]; then
            echo "WARNING: Renaming {input} to {output}"
            cp {input} {output}
        fi
        """
