"""Preprocess NetCDF tropical cyclone fixed return period files

- fix metadata to comply with CF conventions
- write TIFFs for relative convenience - matches other hazard layers
"""
import os
from glob import glob

import rioxarray
import xarray
from preprocess_utils import load_config


def main(config):
    storm_path = os.path.join(
        config["paths"]["data"], "hazards", "TC_data_fixed_return"
    )
    # fix metadata
    for fname in get_fnames(storm_path):
        with xarray.open_dataset(fname) as ds:
            ds = add_meta(ds)
            for var in list(ds):
                for rp in list(ds["mean"].rp.data):
                    out_fname = fname.replace(".nc", f"_rp{int(rp)}_{var}.tif")
                    ds[var].sel(rp=rp).rio.to_raster(out_fname)



def get_fnames(storm_path):
    return [
        fname for fname in glob("*.nc") if "STORM" in fname and "withmeta" not in fname
    ]


def add_meta(ds):
    ds = ds.transpose("rp", "lat", "lon")
    ds.coords["lon"].attrs = {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
        "axis": "X",
    }
    ds.rio.write_crs("epsg:4326", inplace=True)
    return ds


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
