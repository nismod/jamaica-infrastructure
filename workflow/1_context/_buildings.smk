"""Prepare building footprints with attributes for allocation to economic sectors.
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
