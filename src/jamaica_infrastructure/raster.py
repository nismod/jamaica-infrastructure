import os
from collections import namedtuple

import numpy as np
import pandas as pd
import rasterio


# Helper class to store a raster transform and CRS
Transform = namedtuple("Transform", ["crs", "width", "height", "transform"])


def read_transforms(hazards, data_dir):
    transforms = []
    transform_id = 0
    hazard_transforms = []
    for hazard in hazards.itertuples():
        hazard_path = hazard.path
        with rasterio.open(os.path.join(data_dir, hazard_path)) as dataset:
            crs = dataset.crs
            width = dataset.width
            height = dataset.height
            transform = Transform(crs, width, height, tuple(dataset.transform))
        # add transform to list if not present
        if transform not in transforms:
            transforms.append(transform)
            transform_id = transform_id + 1

        # record hazard/transform details
        hazard_transform_id = transforms.index(transform)
        hazard_transform = hazard._asdict()
        del hazard_transform["Index"]
        hazard_transform["transform_id"] = hazard_transform_id
        hazard_transform["width"] = transform.width
        hazard_transform["height"] = transform.height
        hazard_transform["crs"] = str(transform.crs)
        hazard_transform["transform_0"] = transform.transform[0]
        hazard_transform["transform_1"] = transform.transform[1]
        hazard_transform["transform_2"] = transform.transform[2]
        hazard_transform["transform_3"] = transform.transform[3]
        hazard_transform["transform_4"] = transform.transform[4]
        hazard_transform["transform_5"] = transform.transform[5]
        hazard_transforms.append(hazard_transform)
    hazard_transforms = pd.DataFrame(hazard_transforms)

    return hazard_transforms, transforms


def sum_rasters(input_paths: list[str], output_path: str, driver="GTiff") -> None:
    """
    Given a list of raster files on the same grid, sum them and write
    the result to disk.
    """

    base_path, *other_paths = input_paths

    with rasterio.open(base_path) as base_dataset:
        arr = base_dataset.read()[0]

        for path in other_paths:
            with rasterio.open(path) as other_raster:
                other = other_raster.read()[0]
                arr = np.nansum(
                    np.dstack((arr, other)), 2
                )  # nansum to avoid propagating nans over values

        with rasterio.open(
            output_path,
            "w",
            driver=driver,
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            crs=base_dataset.crs,
            transform=base_dataset.transform,
            nodata=base_dataset.nodata,
        ) as output_dataset:
            output_dataset.write(arr, 1)
