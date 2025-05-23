import logging

import click
import numpy as np
import rasterio
from rasterio.transform import from_origin


def gaussian_kernel(shape: tuple[int, int], center: tuple[int, int], bandwidth: float, target_integral: float):
    """
    Generate a 2D Gaussian over grid of `shape`, centered at `center` (y, x),
    with given `bandwidth` (std dev in pixels), normalised to sum to
    `target_integral`.
    """

    rows, cols = shape
    y, x = np.indices((rows, cols))
    kernel = np.exp(-((x - center[1]) ** 2 + (y - center[0]) ** 2) / (2 * bandwidth ** 2))

    integral = kernel.sum()
    if integral > 0:
        kernel *= (target_integral / integral)

    return kernel


@click.command()
@click.version_option("1.0")
@click.option(
    "--input-raster-path", "-i", required=True, help="Path to raster of quantity to smooth.",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--output-raster-path", "-o", required=True, help="Path to write output KDE smoothed raster to.",
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--output-resolution", "-r", required=True, type=float,
    help="Grid resolution of output raster."
)
@click.option(
    "--bandwidth", "-b", required=True, type=float,
    help="Kernel width parameter, the standard deviation in meters.",
)
def main(
    input_raster_path: str,
    output_raster_path: str,
    output_resolution: float,
    bandwidth: float,
):
    logging.info("Loading input raster")
    with rasterio.open(input_raster_path) as src:
        input_data = src.read(1)
        input_transform = src.transform
        input_crs = src.crs
        input_shape = input_data.shape
        input_res = src.res[0]

        # Calculate extent
        # N.B. The transform we're using from the hotspots grid has bottom left
        # as the origin
        left = input_transform.c
        bottom = input_transform.f
        right = left + input_shape[1] * input_res
        top = bottom + input_shape[0] * input_res

    logging.info("Define output grid")
    width = int(np.ceil((right - left) / output_resolution))
    height = int(np.ceil((top - bottom) / output_resolution))
    output_transform = from_origin(left, top, output_resolution, output_resolution)
    output_data = np.zeros((height, width), dtype=np.float32)

    logging.info("Adding Gaussian kernels")
    for input_row in range(input_shape[0]):
        for input_col in range(input_shape[1]):

            loss = input_data[input_row, input_col]
            if not np.isfinite(loss) or loss <= 0:
                continue

            # World coordinates of the center of the input cell
            x_world, y_world = rasterio.transform.xy(
                input_transform, input_row, input_col, offset='center'
            )
            # Output grid indices of that point
            output_col_c, output_row_c = ~output_transform * (x_world, y_world)
            center = (output_row_c, output_col_c)

            # Convert bandwidth to output grid pixels
            bw_pixels: float = bandwidth / output_resolution

            # Add full-image kernel
            output_data += gaussian_kernel((height, width), center, bw_pixels, loss)

    logging.info("Asserting input and output raster sums are equal")
    assert np.isclose(input_data.sum(), output_data.sum(), rtol=1E-3, atol=1)

    logging.info(f"Write out to {output_raster_path}")
    write_kwargs = {
        'driver': 'GTiff',
        'height': height,
        'width': width,
        'count': 1,
        'dtype': 'float32',
        'crs': input_crs,
        'transform': output_transform,
    }
    with rasterio.open(output_raster_path, 'w', **write_kwargs) as dst:
        dst.write(output_data, 1)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    main()
