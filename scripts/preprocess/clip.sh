#!/usr/bin/env bash

# Clip to bounding box
# usage:
#   bash clip.sh <adm3> <xmin> <ymin> <xmax> <ymax>
#
# e.g. for Ghana:
#   bash clip.sh GHA -3.262509 4.737128 1.187968 11.162937

adm3=$1

xmin=$2
ymin=$3
xmax=$4
ymax=$5

# Sense check
echo "Clipping files for $adm3 to bounding box"
echo "  $xmin (xmin)"
echo "  $ymin (ymin)"
echo "  $xmax (xmax)"
echo "  $ymax (ymax)"

# Make container directory
mkdir -p $adm3
# Clip files
ls *.tif | parallel gdalwarp -te $xmin $ymin $xmax $ymax {} $adm3/{}