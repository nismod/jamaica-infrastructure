#
# Access hazard metadata and direct damages from Jamaica-phase-3 processed_data directory
#

pushd  ../../../processed_data/nbs-river-catchment

    cp /mnt/linux-filestore/mistral/Jamaica-phase-3/processed_data/hazards/_hazard_layers__with_transforms.csv .
    cp /mnt/linux-filestore/mistral/Jamaica-phase-3/processed_data/sensitivity_parameters.csv  .
    ln -s /mnt/linux-filestore/mistral/Jamaica-phase-3/results/direct_damages direct_damages

popd

#
# Other requirements
#
# 1. Process DEM to identify stream network and catchments
#   - potentially condition DEM to remove small sinks (r.hydrodem)
#   - calculate flow accumulation, drainage direction, the location of streams and watershed basins (r.watershed)
# 2. Process all possible at-risk points
#   - from extreme flood map (RP 1500), take all flooded cells, convert to points
#   - snap to nearest point on the stream network (r.stream.snap)

# Drainage direction derived from DEM using r.watershed
# - drainage_direction_100.tif

# All possible exposure points derived from each flooded cell of the RP1500 flood raster, with flood_i,j indices
# - non_snapped_points_with_flood_indices.parquet

# Exposure points from the previous step, snapped to river stream network using r.stream.snap, with dem_i,j indices
# - snapped_points_with_dem_indices.parquet

#
# Outputs
#

# Set of snapped points corresponding to exposure points where there are infrastructure damages, to find upstream catchments
# [geometry, dem_i, dem_j]
# - points_to_catchment.parquet

# Relation between snapped points and exposure points
# [flood_i, flood_j, dem_i, dem_j]
# - flood_ij_to_dem_ij_relation.parquet

# Set of catchment geometries
# [geometry, dem_i, dem_j]
# - catchments.parquet
