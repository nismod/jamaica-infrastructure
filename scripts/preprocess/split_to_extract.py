"""Split a data layer to extract using a grid of cells
"""

import os
import sys

import geopandas


def main(grid_fname, cell_number, data_fname, data_layer):
    grid = geopandas.read_file(grid_fname)
    cell = grid.loc[cell_number]
    data = geopandas.read_file(data_fname, layer=data_layer)
    subset = data.cx[cell.left : cell.right, cell.bottom : cell.top].copy()
    subset.geometry = subset.intersection(cell.geometry)
    base, ext = os.path.splitext(data_fname)
    output_fname = f"{base}__{cell_number}.gpkg"
    subset.to_file(output_fname, layer=data_layer, driver="GPKG")


if __name__ == "__main__":
    try:
        grid_fname, cell_number, data_fname, data_layer = sys.argv[1:]
    except:
        print(sys.argv[1:])
        print("Expected usage:")
        print(f"    python {__file__} grid.gpkg 0 data.gpkg layername")
        exit()

    print(f"Reading grid from {grid_fname}")
    print(f"Reading data from {data_fname}|{data_layer}")
    print(f"Processing cell {cell_number}")

    main(grid_fname, int(cell_number), data_fname, data_layer)
