"""Helper to convert parquet to geopackage outputs
"""
import sys

import geopandas


if __name__ == '__main__':
    try:
        _, in_file, out_file, layername = sys.argv
    except IndexError:
        print("Error. Please provide input and output filenames.")
        print(f"Usage: python {__file__} input.geoparquet output.gpkg layername")

    gdf = geopandas.read_parquet(in_file)
    gdf.to_file(out_file, layer=layername, driver="GPKG")
