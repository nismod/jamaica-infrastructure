# Flood-based rasters are in degrees, 0.002
python workflow/4b_hotspots/kde.py \
   -i processed_data/nbs-river-catchment/baseline__fluvial__ead_min.tif \
   -o processed_data/nbs-river-catchment/baseline__fluvial__ead_min_smoothed.tif \
   -b 0.005 # deg
python workflow/4b_hotspots/kde.py \
   -i processed_data/nbs-river-catchment/baseline__fluvial__ead_max.tif \
   -o processed_data/nbs-river-catchment/baseline__fluvial__ead_max_smoothed.tif \
   -b 0.005 # deg
python workflow/4b_hotspots/kde.py \
   -i processed_data/nbs-river-catchment/avoided__fluvial__ead_min.tif \
   -o processed_data/nbs-river-catchment/avoided__fluvial__ead_min_smoothed.tif \
   -b 0.005 # deg
python workflow/4b_hotspots/kde.py \
   -i processed_data/nbs-river-catchment/avoided__fluvial__ead_max.tif \
   -o processed_data/nbs-river-catchment/avoided__fluvial__ead_max_smoothed.tif \
   -b 0.005 # deg

# DEM-based rasters are in metres, 30m cell
python workflow/4b_hotspots/kde.py \
    -i processed_data/nbs-river-catchment/damage_reduction_min.tif \
    -o processed_data/nbs-river-catchment/damage_reduction_min_smoothed.tif \
    -b 300 # m
python workflow/4b_hotspots/kde.py \
    -i processed_data/nbs-river-catchment/damage_reduction_max.tif \
    -o processed_data/nbs-river-catchment/damage_reduction_max_smoothed.tif \
    -b 300 # m
