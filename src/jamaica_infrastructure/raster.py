import rasterio


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
                arr += other_raster.read()[0]

        with rasterio.open(
            output_path,
            "w",
            driver=driver,
            height=arr.shape[0],
            width=arr.shape[1],
            count=1,
            dtype=arr.dtype,
            crs=base_dataset.crs,
            transform=base_dataset.transform
        ) as output_dataset:
            output_dataset.write(arr, 1)
