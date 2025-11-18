# launch this script from a GRASS GIS terminal, with micromamba/conda environment activated:
#   cd scripts/nbs_methods/river-catchment
#   # create GRASS project
#   grass -c EPSG:3448 nbs_project
#   # activate conda environment
#   micromamba activate jsrat
#   # run notebook to import drainage_direction_100.tif and set project region
#   # for more parallelism, consider setting up multiple GRASS projects
#   jupyter notebook 2_grass.init.ipynb
#   cd ./nbs_project
#   grass PERMANENT
#   # run this script (in parallel) against the GRASS project - more than ~2 workers can lead to locking
#   # the last two parameters specify the lower and upper bounds of a slice to work on from points_to_catchment.parquet
#   seq 0 2 | parallel --lb -j 2 python  ../2_grass.r.water.outlet.py {} 2 0 91468
import sys
from pathlib import Path

import geopandas as gpd
import grass.script as gs
from grass.script.core import gisenv
from tqdm.auto import tqdm


def process_catchment(row):
    # easting, northing to get coordinates for GRASS r.water.outlet
    x, y = row.geometry.x, row.geometry.y

    # dem i and j to use in output name and filename
    dem_i, dem_j = row.dem_i, row.dem_j

    base = f"basin_{dem_i}_{dem_j}"
    vect_map = f"{base}_vec"
    vect_gpkg = out_dir / f"{base}.gpkg"

    # return early if GPKG exists
    if vect_gpkg.exists():
        return

    try:
        gs.run_command(
            "r.water.outlet",
            input="drainage_direction_100",
            output=base,
            coordinates=(x, y),
            overwrite=True,
            superquiet=True,
        )

        gs.run_command(
            "r.to.vect",
            input=base,
            output=vect_map,
            type="area",
            overwrite=True,
            superquiet=True,
        )

        gs.run_command(
            "v.out.ogr",
            input=vect_map,
            output=str(vect_gpkg),
            format="GPKG",
            overwrite=True,
            superquiet=True,
        )
    except Exception as e:
        print(f"Failed for ({dem_i},{dem_j}) at ({x},{y}): {e}")


if __name__ == "__main__":
    print(
        "* Example usage:  seq 0 2 | parallel --lb -j 2 python  ../2_grass.r.water.outlet.py {} 2 0 91468"
    )
    worker, nworkers = int(sys.argv[1]), int(sys.argv[2])
    lb, ub = int(sys.argv[3]), int(sys.argv[4])
    print("Worker", worker, "of", nworkers)
    print("From", lb, "to", ub)

    base_path = Path("../../../../processed_data/nbs-river-catchment")
    out_dir = base_path / "upstream_basins_100"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(gisenv())

    points = gpd.read_parquet(base_path / "points_to_catchment.parquet")[lb:ub]

    # Loop over jobs
    for i, row in tqdm(
        enumerate(points.itertuples()),
        total=len(points),
        desc=f"worker {worker}",
    ):
        if ((i + worker) % nworkers) == 0:
            process_catchment(row)
        else:
            continue
