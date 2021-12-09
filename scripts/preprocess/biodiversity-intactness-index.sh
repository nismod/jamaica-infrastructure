#
# Download, convert and clip Biodiversity Indicator Index
#

# Download from source
curl https://data.nhm.ac.uk/dataset/17179b71-c5e1-435c-b3db-ebc7a65d980a/resource/8531b4dc-bd44-4586-8216-47b3b8d60e85/download/lbii.zip > lbii.zip

# Unzip and convert
unzip lbii.zip && rm lbii.zip
gdal_translate lbii.asc lbii.tif -co "COMPRESS=LZW" # compressed TIFF about 700M, similar to the downloaded zip
rm lbii.asc # about 5.4G

# Clip to approximate Jamaica bounds (xmin, ymin, xmax, ymax) coordinates in lat/lon
gdalwarp -te -78.5 17.6 -76.1 18.6 -s_srs EPSG:4326 lbii.tif lbii-jamaica.tif
