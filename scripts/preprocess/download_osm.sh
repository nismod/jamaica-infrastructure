#!/usr/bin/env bash
#
# Download country extract
#

pushd incoming_data/osm

# download extract
wget http://download.geofabrik.de/central-america/jamaica-210622.osm.pbf

# check extract
md5sum --check jamaica-210622.osm.pbf.md5

popd
