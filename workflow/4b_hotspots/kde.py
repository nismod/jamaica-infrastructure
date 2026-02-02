import logging

import click
import numpy as np
import rasterio
import scipy
import tqdm
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
    "--output-resolution", "-r", required=False, type=float,
    help="Grid resolution of output raster."
)
@click.option(
    "--bandwidth", "-b", required=True, type=float,
    help="Kernel width parameter, the standard deviation in meters.",
)
def main(
    input_raster_path: str,
    output_raster_path: str,
    output_resolution: float | None,
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

    if output_resolution is None or np.isclose(input_res, output_resolution, rtol=1e-3):
        logging.info("Smoothing on grid with unchanged resolution")
        output_transform = input_transform
        # Run relatively fast convolution without changing resolution
        bw_pixels: int = int(bandwidth / input_res)
        assert bw_pixels > 0, "Bandwidth must be >= pixel resolution"
        logging.info(f"Define kernel {bandwidth=}, {bw_pixels=}")
        kernel = np.outer(
            scipy.signal.windows.gaussian(bw_pixels * 8, bw_pixels),
            scipy.signal.windows.gaussian(bw_pixels * 8, bw_pixels),
        )
        input_data = np.nan_to_num(input_data)
        output_data = scipy.signal.fftconvolve(input_data, kernel, mode="same")

        # Normalise back to preserve total
        sum_in = np.nansum(input_data)
        sum_out = np.nansum(output_data)
        output_data *= sum_in / sum_out
        logging.debug(f"{sum_in=}, {sum_out=}, {output_data=}")

    else:
        logging.info("Smoothing on grid with {output_resolution=}")
        width = int(np.ceil((right - left) / output_resolution))
        height = int(np.ceil((top - bottom) / output_resolution))
        output_transform = from_origin(left, top, output_resolution, output_resolution)
        output_data = np.zeros((height, width), dtype=np.float32)

        logging.info("Adding Gaussian kernels")
        for input_row in tqdm.tqdm(range(input_shape[0])):
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
    assert np.isclose(np.nansum(input_data), np.nansum(output_data), rtol=1E-3, atol=1)

    logging.info(f"Write out to {output_raster_path}")
    write_kwargs = {
        'driver': 'GTiff',
        'height': output_data.shape[0],
        'width': output_data.shape[1],
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
