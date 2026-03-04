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

rule parse_nic_2005_building_sectors:
    input:
        script = "workflow/1_context/parse_nic_2005_building_sectors.py",
        buildings = f"{DATA}/buildings/buildings_assigned_economic_activity.gpkg",
    output:
        buildings = f"{DATA}/buildings/buildings_nic2005.geoparquet",
    shell:
        """
        python {input.script} \
            --buildings {input.buildings} \
            --output {output.buildings}
        """
