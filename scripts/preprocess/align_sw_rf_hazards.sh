#
# Align river flooding and surface water hazards on the same grid
#
# The data share CRS and resolution, just clipped to different height/width by
# a few pixels.
#
# Run from the processed_data/hazards directory
#

ls Global\ Flood\ Map/Jamaica/*/Raw\ Depths/*02.tif | sed 's/.tif//' | \
parallel gdalwarp \
    -t_srs EPSG:4326 \
    -te -78.3690277774130948 17.7129166663936637 -75.9693055359931009 18.5254166728936589 \
    -ts 8639 2925 \
    -r near \
    {}.tif {}-aligned.tif
