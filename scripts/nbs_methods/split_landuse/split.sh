python make_valid.py

SNAIL_PROGRESS=TRUE snail -vv split -f 2013_landuse_landcover.gpkg -t 1000 0 605000 0 -1000 708000 --width 300 --height 200 -o split_1k.gpkg
SNAIL_PROGRESS=TRUE snail -vv split -f 2013_landuse_landcover.gpkg -t 10000 0 605000 0 -10000 708000 --width 30 --height 20 -o split_10k.gpkg
