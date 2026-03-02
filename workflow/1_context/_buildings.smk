"""Prepare building footprints with attributes for allocation to economic sectors.
"""

rule download_osm_jamaica:
    """Download country extract
    """
    output:
        pbf = f"{RAW}/osm/jamaica-{config["context"]["osm_extract_date"]}.osm.pbf",
    shell:
        """
        extract_dir=$(dirname {output.pbf})
        extract_date={config[context][osm_extract_date]}

        pushd $extract_dir

            # download extract
            wget http://download.geofabrik.de/central-america/jamaica-${{extract_date}}.osm.pbf
            wget http://download.geofabrik.de/central-america/jamaica-${{extract_date}}.osm.pbf.md5

            # check extract
            md5sum --check jamaica-${{extract_date}}.osm.pbf.md5

        popd
        """


rule filter_osm_buildings:
    """Filter OpenStreetMap extract for all buildings in Jamaica
    """
    input:
        pbf = f"{RAW}/osm/jamaica-{config["context"]["osm_extract_date"]}.osm.pbf",
    output:
        pbf = f"{DATA}/buildings/jamaica-buildings.osm.pbf",
        gpkg = f"{DATA}/buildings/jamaica-buildings.gpkg",
    shell:
        """
        osmium tags-filter \
            {input.pbf} \
            wnr/building \
            --overwrite \
            -o {output.pbf}

        # Extract features from .osm.pbf to .gpkg
        OSM_CONFIG_FILE=workflow/1_context/osmconf.ini ogr2ogr -f GPKG \
            {output.gpkg} \
            {output.pbf} \
            multipolygons
        """

rule building_height_assignments:
    """
    """
    input:
        script = "workflow/1_context/building_height_assignments.py",
        # It also needs the below script to exist to run a subprocess, not sure if that would be stated as an input or param or not at all
        vector_raster_intersections_script = "workflow/1_context/vector_raster_intersections.py",
        building_layers = f"{DATA}/buildings/building_layers.csv",
        building_rasters = f"{DATA}/buildings/building_rasters.csv",
        buildings_econ_activity = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg"
    output:
        results = f"{OUTPUT}/buildings_heights/buildings_heights.gpkg"
    shell:
        """
        python {input.script} \
            --data-dir {DATA} \
            --output-dir {OUTPUT} \
            --vector-raster-intersections-script {input.vector_raster_intersections_script} \
            --building-layers {input.building_layers} \
            --building-rasters {input.building_rasters} \
            --buildings-econ-activity {input.buildings_econ_activity} 
        """

rule land_use_layers_process:
    """
    """
    input:
        script = "workflow/1_context/land_use_layers_process.py",
        tnc_landuse = f"{RAW}/nsdmb/GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb",
        forest_landuse = f"{RAW}/Landuse 2013 data/2013_landuse_Landcover.shp",
        mining_landuse = f"{RAW}/global_mining_areas/global_mining_polygons_v1.gpkg",
        forest_sector_mapping = f"{DATA}/land_type_and_use/forest_classes_with_sector_mapping.csv",
        tnc_sector_mapping = f"{DATA}/land_type_and_use/tnc_classes_with_sector_mapping.csv"
    output:
        input_land_use_layers = f"{DATA}/land_type_and_use/input_land_use_layers.gpkg",
        jamaica_land_use_combined = f"{DATA}/land_type_and_use/jamaica_land_use_combined.gpkg",
        jamaica_land_use_combined_with_sectors = f"{DATA}/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg",
        #forest_classes = f"{DATA}/land_type_and_use/forest_classes.csv",
        #tnc_classes = f"{DATA}/land_type_and_use/tnc_classes.csv",
        #land_use_modified = f"{DATA}/land_type_and_use/land_use_modified.gpkg",
    shell:
       """
        python {input.script} \
            --incoming-data-dir {RAW} \
            --data-dir {DATA} \
            --tnc-landuse-path {input.tnc_landuse} \
            --forest-landuse-path {input.forest_landuse} \
            --mining-landuse-path {input.mining_landuse} \
            --forest-sector-mapping-path {input.forest_sector_mapping} \
            --tnc-sector-mapping-path {input.tnc_sector_mapping} \
        """

rule planning_layers:
    """
    """
    input:
        script = "workflow/1_context/planning_layers.py",
        planning_layer_polygons = f"{RAW}/buildings/nsdmb_planning_layers_polygons.xlsx",
        planning_layers = f"{RAW}/buildings/planning_layers.gpkg",
        manchester_layers_path = f"{RAW}/buildings/land_use_layers_uses_codes.csv",
        jamaica_parishes_path = f"{DATA}/boundaries/admin_boundaries.gpkg",

        #there is an array that also defines the following files to be read:
        clarendon_landuse_classes_with_sector_codes = f"{RAW}/buildings/clarendon_landuse_planning_classes_with_sector_codes.csv",
        machester_landuse_classes_with_sector_codes = f"{RAW}/buildings/manchester_landuse_planning_classes_with_sector_codes.csv",
        landuse_classes_with_sector_codes = f"{RAW}/buildings/landuse_classes_with_sector_codes.csv",
    output:
        #landuse_planning_layers.gpkg = f"{RAW}/buildings/landuse_planning_layers.gpkg"
        landuse_planning_layers_with_sectors = f"{RAW}/buildings/landuse_planning_layers_with_sectors.gpkg"
    shell:
       """
        python {input.script} \
            --incoming-data-dir {RAW} \
            --data-dir {DATA} \
            --planning-layer-polygons {input.planning_layer_polygons} \
            --planning-layers-path {input.planning_layers} \
            --manchester-layers-path {input.manchester_layers_path} \
            --jamaica-parishes-path {input.jamaica_parishes_path} \
            --clarendon-landuse-classes-with-sector-codes {input.clarendon_landuse_classes_with_sector_codes} \
            --machester-landuse-classes-with-sector-codes {input.machester_landuse_classes_with_sector_codes} \
            --landuse-classes-with-sector-codes {input.landuse_classes_with_sector_codes} \
        """

rule building_data_process:
    """
    """
    input:
        script = "workflow/1_context/buildings_data_process.py",
        poi_data = f"{RAW}/hotosm_data/hotosm_jam_points_of_interest_gpkg/hotosm_jam_points_of_interest.gpkg",
        poi_mapping = f"{RAW}/hotosm_data/hotosm_jam_points_of_interest_gpkg/hotosm_mapping_economic_sectors.xlsx",
        buildings_input_map = f"{RAW}/buildings/osm_building_column_mapping_economic_sectors.xlsx",
        buildings_input = f"{RAW}/buildings/jamaica-buildings.gpkg",
        landuse_planning_layers = f"{RAW}/buildings/landuse_planning_layers_with_sectors.gpkg", #generated by planning_layers script
        landuse_types = f"{DATA}/land_type_and_use/jamaica_land_use_combined.gpkg", #generated by land_use_layers_process script
        landuse_types_with_sectors = f"{DATA}/land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg", #generated by land_use_layers_process script
        commercial_buildings = f"{RAW}/JAM_classified_buildings_industries/JAM_classified_buildings.shp",
        nsdmb_layers = f"{RAW}/buildings/nsdmb_layers_economic_sectors.xlsx",
        ports = f"{DATA}/networks/transport/port_polygon.gpkg",
        airports = f"{DATA}/networks/transport/airport_polygon.gpkg",
        residential_buildings = f"{RAW}/JAM_residential/jam_residential.shp",
        known_building_assign = f"{RAW}/buildings/manual_building_assign.shp",
        fishing_locations = f"{DATA}/land_type_and_use/aqua_farms.gpkg",
        population = f"{DATA}/population/population_projections.gpkg",

        #there is a layer_dict defined in the script that also seems to use the following:
        trade_plants = f"{RAW}/buildings/nsdmb_tradeplants_layer_economic_sectors.xlsx",
        exporters = f"{RAW}/buildings/nsdmb_exporters_layer_economic_sectors.xlsx",
        industry = f"{RAW}/buildings/nsdmb_industry_layer_economic_sectors.xlsx",

        #there are some cost shapefiles also in an array:
        building_costs_1 = f"{RAW}/construction_permits/2019-2020/Construction_Permit_Mapping_Jan-March_Q4_2019-2020.shp",
        building_costs_2 = f"{RAW}/construction_permits/2019-2020/Construction_Permits_October_-_December_Q3_2019-2020.shp",
        building_costs_3 = f"{RAW}/construction_permits/2019-2020/Construction_Permits_Q1_2019_2020.shp",
        building_costs_4 = f"{RAW}/construction_permits/2019-2020/Q2-2019-Construction-Permits.shp",

    params:
        residential_min_area = 11.15, #perhaps should be added to config.yaml
    output:
        #there are some intermediary files that get generated (not sure if they have to be in output):
        #incoming_data/buildings/commercial_buildings_assigned_economic_sectors_intermediate.gpkg
        #incoming_data/hotosm_data/hotosm_jam_points_of_interest_gpkg/points_of_interest_assigned_economic_sectors_intermediate.gpkg
        #incoming_data/buildings/building_attributes.csv
        #incoming_data/buildings/buidling_sector_subsectors.csv
        #incoming_data/buildings/building_cost_ranges.csv
        intermediate_file = f"{RAW}/buildings/buildings_assigned_economic_sectors_intermediate.gpkg", #seems to be required by socio economic model scripts
        result = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg"

    shell:
        """
        python {input.script} \\
            --incoming-data-dir {RAW} \\
            --data-dir {DATA} \\
            --residential-min-area {params.residential_min_area} \\
            --poi-data {input.poi_data} \\
            --poi-mapping {input.poi_mapping} \\
            --buildings-input-map {input.buildings_input_map} \\
            --buildings-input-file {input.buildings_input} \\
            --landuse-planning-layers-file {input.landuse_planning_layers} \\
            --landuse-types-file {input.landuse_types} \\
            --landuse-types-with-sectors-file {input.landuse_types_with_sectors} \\
            --commercial-buildings-file {input.commercial_buildings} \\
            --nsdmb-layers-file {input.nsdmb_layers} \\
            --ports-file {input.ports} \\
            --airports-file {input.airports} \\
            --residential-buildings-file {input.residential_buildings} \\
            --known-building-assign-file {input.known_building_assign} \\
            --fishing-locations-file {input.fishing_locations} \\
            --population-file {input.population} \\
            --trade-plants-file {input.trade_plants} \\
            --exporters-file {input.exporters} \\
            --industry-file {input.industry} \\
            --building-costs-1 {input.building_costs_1} \\
            --building-costs-2 {input.building_costs_2} \\
            --building-costs-3 {input.building_costs_3} \\
            --building-costs-4 {input.building_costs_4}
        """

# rule buildings_sector_allocation:
#     """
#     The script associted with this rule SEEMS to functionally do something similar to the buildings data process rule
#     but it seems to be far simpler. Perhaps it was an earlier version?
#     """
#     input:
#         script = "workflow/1_context/buildings_sector_allocations.py",
#         buildings_input = f"{RAW}/buildings/jamaica-buildings.gpkg",
#         buildings_input_map = f"{RAW}/buildings/osm_building_column_mapping_economic_sectors.xlsx",
#         nsdmb_layers = f"{RAW}/buildings/nsdmb_layers_economic_sectors.xlsx",
#         trade_plants = f"{RAW}/buildings/nsdmb_tradeplants_layer_economic_sectors.xlsx",
#         exporters = f"{RAW}/buildings/nsdmb_exporters_layer_economic_sectors.xlsx",
#         industry = f"{RAW}/buildings/nsdmb_industry_layer_economic_sectors.xlsx",
#     params:
#         epsg = config["adaptation_options"]["epsg_jamaica"],
#     output:
#         unique_column_mapping = f"{RAW}/buildings/unique_column_mapping.xlsx",
#         #incoming_data/buildings/buildings_assigned_economic_sectors_intermediate.gpkg
#         #incoming_data/buildings/test_values.csv
#         #incoming_data/buildings/industry_points.csv
#         result = f"{DATA}/buildings/buildings_assigned_economic_sectors.gpkg"
#     shell:
#        """
#         python {input.script} \
#             --incoming-data {RAW} \
#             --data-dir {DATA} \
#             --epsg {epsg} \


#         """



