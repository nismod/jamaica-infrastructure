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

#
# Coarsen then smooth avoided damages
#
gdalwarp -overwrite -t_srs EPSG:3448 processed_data/nbs-river-catchment/avoided__fluvial__ead_min.tif processed_data/nbs-river-catchment/avoided__fluvial__ead_min_3448.tif
gdalwarp -overwrite -t_srs EPSG:3448 processed_data/nbs-river-catchment/avoided__fluvial__ead_max.tif processed_data/nbs-river-catchment/avoided__fluvial__ead_max_3448.tif

gdalwarp -overwrite -r sum -tr 300 300 -t_srs EPSG:3448 processed_data/nbs-river-catchment/avoided__fluvial__ead_min_3448.tif processed_data/nbs-river-catchment/avoided__fluvial__ead_min_300m.tif
gdalwarp -overwrite -r sum -tr 300 300 -t_srs EPSG:3448 processed_data/nbs-river-catchment/avoided__fluvial__ead_max_3448.tif processed_data/nbs-river-catchment/avoided__fluvial__ead_max_300m.tif

python workflow/4b_hotspots/kde.py -i processed_data/nbs-river-catchment/avoided__fluvial__ead_min_300m.tif -o processed_data/nbs-river-catchment/avoided__fluvial__ead_min_300m_smoothed.tif -b 300
python workflow/4b_hotspots/kde.py -i processed_data/nbs-river-catchment/avoided__fluvial__ead_max_300m.tif -o processed_data/nbs-river-catchment/avoided__fluvial__ead_max_300m_smoothed.tif -b 300
