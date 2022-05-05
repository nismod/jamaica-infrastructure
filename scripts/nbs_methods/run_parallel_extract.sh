# split land use on 10km grid
seq 1 149 | parallel -j20 'python ~/jamaica-infrastructure/scripts/preprocess/split_to_extract.py nbs/grid_10km.gpkg {} land_type_and_use/jamaica_land_use_combined_with_sectors.gpkg areas'

# link each to slope
seq 1 149 | parallel -j20 'python   ~/jamaica-infrastructure/scripts/preprocess/link_landuse_slope.py   nbs/slope.tif   land_type_and_use/jamaica_land_use_combined_with_sectors__{}.gpkg   nbs/jamaica_land_use_combined_with_sectors_and_slope__{}.gpkg   forest_id,Classify,LU_CODE,tnc_id,NAME,TNCCODE,global_id,global_LU_type,sector_code_forest,subsector_code_forest,infra_sector_forest,sector_code_tnc,subsector_code_tnc,infra_sector_tnc'

# join to other vector layers
# NB off by one!
seq 2 150 | parallel -j20 'python ~/jamaica-infrastructure/scripts/nbs_methods/join-layers.py {}'

# convert to GeoJSONSeq
find -name "combined*gpkg" | sed 's/.gpkg//' | parallel -j20 ogr2ogr -t_srs "EPSG:4326" -f "GeoJSONSeq" "{}.geojsonl" "{}.gpkg"

# Combine to single GeoJSONSeq file, copy to visualisation tool folder
cat land_use_slope/extracts/combined__*l > land_use_slope/combined_polygons.geojsonl
cp land_use_slope/combined_polygons.geojsonl ~/infra-risk-vis/etl/vector/natural_terrestrial_combined.geojsonl

# Rename and convert to GPKG
mv land_use_slope/combined_polygons.geojsonl land_use_slope/natural_terrestrial_combined.geojsonl
ogr2ogr -t_srs "EPSG:3448" -f "GPKG" land_use_slope/natural_terrestrial_combined.gpkg land_use_slope/natural_terrestrial_combined.geojsonl
