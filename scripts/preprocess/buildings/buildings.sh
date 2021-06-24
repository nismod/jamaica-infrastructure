#!/usr/bin/env bash
#
# Extract buildings from OSM
#
set -e
set -x

# Extract date string
date="210622"

# Extract features from .osm.pbf to .gpkg
osmium tags-filter \
    incoming_data/osm/jamaica-$date.osm.pbf \
    wnr/building \
    --overwrite \
    -o incoming_data/osm/jamaica-buildings.osm.pbf

OSM_CONFIG_FILE=scripts/preprocess/buildings/osmconf.ini ogr2ogr -f GPKG \
    incoming_data/osm/jamaica-buildings.gpkg \
    incoming_data/osm/jamaica-buildings.osm.pbf \
    multipolygons
