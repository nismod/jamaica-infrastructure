#!/bin/bash

#
# Run surrounding_landcover.py across landcover class/groups and buffer distances
#

## Sense-check parallel command - echo only
# parallel -j 1 'echo class {1} buffer {2}' \
#     ::: primary-forest secondary-forest dry-forest wetland plantations fields agriculture \
#     ::: 5 100 250 500 1000

parallel -j 30 'python surrounding_landcover.py "~/data/jamaica-local/incoming/Terrestrial Land Cover (From Forestry Department)" ~/data/jamaica-local/results/nbs {1} {2}' \
    ::: primary-forest secondary-forest dry-forest wetland plantations fields agriculture \
    ::: 5 100 250 500 1000
