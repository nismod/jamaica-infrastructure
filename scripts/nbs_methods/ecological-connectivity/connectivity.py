#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 15:11:06 2023
Edited on Sun Feb 18 by sarahgall

@author: matthiaswildemeersch
... and nudged along by Fred -- 2024-03-11
... tomalrussell 2026-01-12
"""

import math
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numba import njit, prange


HALF_ROOT_TWO = np.sqrt(2) / 2.0


def _gen_mode_flag(gen_mode: str, number_of_species_gens: int) -> int:
    if gen_mode == "one_generation":
        return 0
    if gen_mode == "multi-generation":
        return 1 if number_of_species_gens < 10 else 2
    raise ValueError(f"Unsupported generation mode: {gen_mode}")


@njit(cache=True)
def _regularized_gamma_upper_integer(k: int, x: float) -> float:
    """
    Regularized upper incomplete gamma Q(k, x) for integer k.

    Notes:
        Q(k, x) = e^{-x} * sum_{n=0}^{k-1} x^n / n!
    """
    if k <= 0:
        return 0.0
    total = 0.0
    term = 1.0
    for n in range(k):
        if n > 0:
            term *= x / n
        total += term
    return math.exp(-x) * total


def identify_sector(row: int, col: int) -> int:
    """
    Find the sector that a given target cell lies in. Calculate forward bearing
    in radians from origin to the cell defined by supplied row and column
    indices. Use a bearing to sector lookup to find the sector.

    Args:
        row: Row index of target cell, increasing southwards
        col: Column index of target cell, increasing eastwards

    Returns:
        Sector identifier
    """

    # -row+col 0-0.5PI  N E
    # +row+col 0.5PI - 1PI  S E
    # +row-col 1PI-1.5 PI  S W
    # -row-col 1.5PI-2PI  N W

    # TODO: fix behaviour whereby rows in the southern half of the column containing bullseye are labelled as north not south
    # for visual explanation see seg_bad.png
    # need a new conditional clause for col == 0, or a better way of computing the bearing
    bearing_rad: float = 0

    if row < 0:
        if col < 0:
            # -row-col 1.5PI-2PI NORTH WEST
            bearing_rad = 2 * math.pi - math.atan(-col / -row)
        elif col > 0:
            # -row+col 0- 0.5PI  NORTH EAST
            bearing_rad = math.atan(col / -row)
        elif col == 0:
            # NORTH
            bearing_rad = 0
    elif row > 0:
        if col < 0:
            # +row-col 1PI-1.5 PI  SOUTH WEST
            bearing_rad = math.pi + math.atan(-col / row)
        elif col > 0:
            # +row+col 0.5PI - 1PI  SOUTH EAST
            bearing_rad = math.pi - math.atan(col / row)
        elif col == 0:
            # SOUTH
            bearing_rad = math.pi
    elif row == 0:
        if col < 0:
            bearing_rad = 1.5 * math.pi  # WEST
        elif col > 0:
            bearing_rad = math.pi / 2  # EAST
        elif col == 0:
            bearing_rad = -1  # Self

    # Cardinal bearings (buffered!) N 0 NE 1 E 2 SE 3 S 4 SW 5 W 6 NW 7
    # Original approach used paired bearings N & S= 0 NE & SW=1 etc for alignment

    sector: int = 0

    if bearing_rad < 0.392699082:
        return 0  # N
    elif bearing_rad < 1.178097245:
        return 1  # NE
    elif bearing_rad < 1.963495408:
        return 2  # E
    elif bearing_rad < 2.748893572:
        return 3  # SE
    elif bearing_rad < 3.534291735:
        return 4  # S
    elif bearing_rad < 4.319689899:
        return 5  # SW
    elif bearing_rad < 5.105088062:
        return 6  # W
    elif bearing_rad < 5.890486225:
        return 7  # NW
    else:
        return 0  # N


def find_ring_and_sector(
    row: int, col: int, ring_radii: np.ndarray, dartboard_radius: int
) -> "tuple[int, int]":
    """
    Determine which part of the dartboard we're in.

    Args:
        row: Row index
        col: Column index
        ring_radii: Iterable of dartboard ring radii

    Returns:
        Index of sector, index of ring
    """
    relative_row = (
        row - dartboard_radius
    )  # i_r_row is the horizintal distance to the dart center cell, center is in row index 30
    relative_col = col - dartboard_radius
    # relative_row = row - 30  # i_r_row is the horizintal distance to the dart center cell, center is in row index 30
    # relative_col = col - 30

    # geometrical distance from the regarded cell to the dart center
    distance: float = np.sqrt(relative_row**2 + relative_col**2)

    if distance > np.max(ring_radii):
        # outside the biggest circle
        # relative area of circle inside (touching) square is pi/4 or 78%, would expect ~22% of iterations to end here
        return (-1, -1)

    # sector index
    sector: int = identify_sector(relative_row, relative_col)

    # ring index
    ring: int = 0
    while distance > ring_radii[ring]:
        ring += 1
        if ring == len(ring_radii):
            break

    return sector, ring


@njit(cache=True)
def path_length(
    f_PCurrSum: np.ndarray,
    f_PLastSum: np.ndarray,
    f_rad_width: np.ndarray,
    f_perm_scale: float,
    f_min_perm: float,
    f_C: np.ndarray,
    sector: int,
    ring: int,
) -> np.ndarray:
    """
    Find the shortest path length from the centre of the dartboard to the edge,
    traversing each ring. Paths cannot travel purely circumferetially. They can
    move radially outwards, or to the sector either side of that which is
    radially outwards.

    Args:
        f_PCurrSum: Shortest path distance, current sum
        f_PLastSum: Shortest path distance, sum from previous (inner) ring
        f_rad_width: Ring thicknesses
        f_perm_scale: ?
        f_min_perm: ?
        f_C: Average zone condition (mean of intersecting cells' condition)
        sector: Index of sector
        ring: Index of ring

    Returns:
        Updated version of input matrix, f_PCurrSum, with new value at `sector`.
    """
    # index of direction sector to the left of regarded sector
    i_left_seg = sector - 1 if sector > 0 else 7

    # index of direction sector to the right of regarded sector
    i_right_seg = (sector + 1) % 8

    # A and B represent two paths to reach a zone
    # following lines match with path at the left of the source cell #goes left
    f_PtestA = f_PLastSum[i_left_seg] + f_rad_width[
        ring
    ] / (  # effective distance of straight path
        (f_perm_scale * f_C[i_left_seg][ring]) + f_min_perm
    )

    f_PtestB = f_PLastSum[sector] + HALF_ROOT_TWO * (  # diagonal path
        f_rad_width[ring] + f_rad_width[ring]
    ) / ((f_perm_scale * f_C[i_left_seg][ring]) + f_min_perm)

    if f_PtestA < f_PtestB:
        f_PCurrSum[i_left_seg] = (
            f_PtestA  # PCurrSum is the length of the shorter of the two paths
        )
    else:
        f_PCurrSum[i_left_seg] = f_PtestB

    # following lines correspond with path at the right of source cell
    # two path going to the right cell, compare which one is better   # formula 4 in paper
    f_PtestA = f_PLastSum[i_right_seg] + f_rad_width[ring] / (
        (f_perm_scale * f_C[i_right_seg][ring]) + f_min_perm
    )  # goes straight outside
    f_PtestB = f_PLastSum[sector] + HALF_ROOT_TWO * (
        f_rad_width[ring] + f_rad_width[ring]
    ) / (
        (f_perm_scale * f_C[i_right_seg][ring]) + f_min_perm
    )  # goes diagonal

    # if path a is shorter, then the current Psum[right_seg]  = effective distance A
    if f_PtestA < f_PtestB:
        f_PCurrSum[i_right_seg] = f_PtestA
    else:
        f_PCurrSum[i_right_seg] = f_PtestB

    # following lines correspond with path in the centre of source cell #centre
    f_PtestA = f_PLastSum[sector] + f_rad_width[ring] / (
        (f_perm_scale * f_C[sector][ring]) + f_min_perm
    )
    f_PtestB = f_PLastSum[i_left_seg] + HALF_ROOT_TWO * (
        f_rad_width[ring] + f_rad_width[ring]
    ) / ((f_perm_scale * f_C[sector][ring]) + f_min_perm)
    f_PtestC = f_PLastSum[i_right_seg] + HALF_ROOT_TWO * (
        f_rad_width[ring] + f_rad_width[ring]
    ) / ((f_perm_scale * f_C[sector][ring]) + f_min_perm)

    if f_PtestA < f_PtestB:
        f_Pbest = f_PtestA
    else:
        f_Pbest = f_PtestB

    if f_Pbest < f_PtestC:
        f_PCurrSum[sector] = f_Pbest
    else:
        f_PCurrSum[sector] = f_PtestC

    return f_PCurrSum


@njit(cache=True)
def zone_condition_sum_and_cell_count(
    n_sectors: int,
    n_rings: int,
    sector_indicies: np.ndarray,
    ring_indicies: np.ndarray,
    condition: np.ndarray,
) -> "tuple[np.ndarray, np.ndarray]":
    condition_by_zone = np.zeros((n_sectors, n_rings))
    cell_count_by_zone = np.zeros((n_sectors, n_rings))

    n_rows, n_cols = sector_indicies.shape
    for r in range(n_rows):
        for c in range(n_cols):
            sector = sector_indicies[r, c]
            ring = ring_indicies[r, c]
            if sector < 0 or ring < 0:
                continue
            condition_by_zone[sector, ring] += condition[r, c]
            cell_count_by_zone[sector, ring] += 1.0

    return condition_by_zone, cell_count_by_zone


@njit(cache=True)
def connectivity_of_cell_core(
    sector_index_subset: np.ndarray,
    ring_index_subset: np.ndarray,
    n_sectors: int,
    n_rings: int,
    condition_subset: np.ndarray,
    f_lambda: np.ndarray,
    gen_mode_flag: int,
    number_of_species_gens: int,
    f_rad_width: np.ndarray,
    f_perm_scale: float,
    f_min_perm: float,
) -> float:
    f_H, f_zone_cells = zone_condition_sum_and_cell_count(
        n_sectors,
        n_rings,
        sector_index_subset,
        ring_index_subset,
        condition_subset,
    )

    f_Hlmax = f_zone_cells.copy()
    f_C = np.zeros((n_sectors, n_rings))
    f_Clmax = np.ones((n_sectors, n_rings))

    for sector in range(n_sectors):
        for ring in range(n_rings):
            cells = f_zone_cells[sector, ring]
            if cells > 0:
                f_C[sector, ring] = f_H[sector, ring] / cells
            else:
                f_C[sector, ring] = 0.0

    n_lambda = f_lambda.shape[0]
    d_Wlmaxsum = np.zeros(n_lambda)
    d_Wsum = np.zeros(n_lambda)
    f_PLastSum = np.zeros(n_sectors)
    f_PCurrSum = np.zeros(n_sectors)

    for sector in range(n_sectors):
        for i in range(n_sectors):
            f_PLastSum[i] = 0.0
            f_PCurrSum[i] = 0.0
        f_Psum = 0.0
        f_Plmaxsum = 0.0

        for ring in range(n_rings):
            if f_Clmax[sector, ring] <= 0:
                break

            for i in range(n_sectors):
                f_PCurrSum[i] = 0.0

            f_PCurrSum = path_length(
                f_PCurrSum,
                f_PLastSum,
                f_rad_width,
                f_perm_scale,
                f_min_perm,
                f_C,
                sector,
                ring,
            )
            f_Psum = f_PCurrSum[sector]
            f_PLastSum, f_PCurrSum = f_PCurrSum, f_PLastSum

            f_Plmaxsum += f_rad_width[ring] / (
                (f_perm_scale * f_Clmax[sector, ring]) + f_min_perm
            )

            for i_dist in range(n_lambda):
                lambda_val = f_lambda[i_dist]
                d_dist = f_Psum / lambda_val
                d_dist2 = f_Plmaxsum / lambda_val

                if gen_mode_flag == 0:
                    kernel_function1 = math.exp(-d_dist)
                    kernel_function2 = math.exp(-d_dist2)
                elif gen_mode_flag == 1:
                    kernel_function1 = _regularized_gamma_upper_integer(
                        number_of_species_gens, d_dist
                    )
                    kernel_function2 = _regularized_gamma_upper_integer(
                        number_of_species_gens, d_dist2
                    )
                else:
                    denom = math.sqrt(2.0 * number_of_species_gens) * lambda_val
                    arg1 = (f_Psum - number_of_species_gens * lambda_val) / denom
                    arg2 = (f_Plmaxsum - number_of_species_gens * lambda_val) / denom
                    kernel_function1 = 0.5 * (1.0 - math.erf(arg1))
                    kernel_function2 = 0.5 * (1.0 - math.erf(arg2))

                d_Wsum[i_dist] += kernel_function1 * f_H[sector, ring]
                d_Wlmaxsum[i_dist] += kernel_function2 * f_Hlmax[sector, ring]

    d_Conn = 0.0
    for i_dist in range(n_lambda):
        if d_Wlmaxsum[i_dist] > 0:
            d_Conn += d_Wsum[i_dist]

    if n_lambda > 0:
        d_Conn /= n_lambda

    return d_Conn


def connectivity_of_cell(
    row: int,
    col: int,
    sector_index_by_cell: np.ndarray,
    ring_index_by_cell: np.ndarray,
    n_sectors: int,
    n_rings: int,
    h: np.ndarray,
    f_lambda: np.ndarray,
    gen_mode: str,
    number_of_species_gens: int,
    dartboard_radius: int,
    dartboard_d: int,
    f_rad_width: np.ndarray,
    f_perm_scale: float,
    f_min_perm: float,
) -> float:
    """
    Thin wrapper around the numba-accelerated per-cell solver that extracts the
    dartboard slice to pass into the compiled core implementation.
    """

    n_rows, n_cols = h.shape
    row_min = max(0, row - dartboard_radius)
    row_min_diff = min(0, row - dartboard_radius)
    row_max = min(n_rows, row + dartboard_radius + 1)
    row_max_diff = max(n_rows, row + dartboard_radius + 1) - n_rows
    col_min = max(0, col - dartboard_radius)
    col_min_diff = min(0, col - dartboard_radius)
    col_max = min(n_cols, col + dartboard_radius + 1)
    col_max_diff = max(n_cols, col + dartboard_radius + 1) - n_cols

    sector_slice = sector_index_by_cell[
        -row_min_diff : dartboard_d - row_max_diff,
        -col_min_diff : dartboard_d - col_max_diff,
    ]
    ring_slice = ring_index_by_cell[
        -row_min_diff : dartboard_d - row_max_diff,
        -col_min_diff : dartboard_d - col_max_diff,
    ]
    condition_slice = h[row_min:row_max, col_min:col_max]

    mode_flag = _gen_mode_flag(gen_mode, number_of_species_gens)

    return connectivity_of_cell_core(
        sector_slice,
        ring_slice,
        n_sectors,
        n_rings,
        condition_slice,
        f_lambda,
        mode_flag,
        number_of_species_gens,
        f_rad_width,
        f_perm_scale,
        f_min_perm,
    )


@njit(parallel=True)
def connectivity_grid_kernel(
    condition: np.ndarray,
    land_cells: np.ndarray,
    sector_index_by_cell: np.ndarray,
    ring_index_by_cell: np.ndarray,
    n_sectors: int,
    n_rings: int,
    dartboard_radius: int,
    dartboard_d: int,
    f_lambda: np.ndarray,
    gen_mode_flag: int,
    number_of_species_gens: int,
    f_rad_width: np.ndarray,
    f_perm_scale: float,
    f_min_perm: float,
) -> np.ndarray:
    n_points = land_cells.shape[0]
    results = np.zeros(n_points)
    n_rows, n_cols = condition.shape

    for idx in prange(n_points):
        row = land_cells[idx, 0]
        col = land_cells[idx, 1]

        row_min = row - dartboard_radius
        row_min_diff = 0
        if row_min < 0:
            row_min_diff = row_min
            row_min = 0

        row_max = row + dartboard_radius + 1
        row_max_diff = 0
        if row_max > n_rows:
            row_max_diff = row_max - n_rows
            row_max = n_rows

        col_min = col - dartboard_radius
        col_min_diff = 0
        if col_min < 0:
            col_min_diff = col_min
            col_min = 0

        col_max = col + dartboard_radius + 1
        col_max_diff = 0
        if col_max > n_cols:
            col_max_diff = col_max - n_cols
            col_max = n_cols

        sector_slice = sector_index_by_cell[
            -row_min_diff : dartboard_d - row_max_diff,
            -col_min_diff : dartboard_d - col_max_diff,
        ]
        ring_slice = ring_index_by_cell[
            -row_min_diff : dartboard_d - row_max_diff,
            -col_min_diff : dartboard_d - col_max_diff,
        ]
        condition_slice = condition[row_min:row_max, col_min:col_max]

        results[idx] = connectivity_of_cell_core(
            sector_slice,
            ring_slice,
            n_sectors,
            n_rings,
            condition_slice,
            f_lambda,
            gen_mode_flag,
            number_of_species_gens,
            f_rad_width,
            f_perm_scale,
            f_min_perm,
        )

    return results


def connectivity_of_grid(
    condition: np.ndarray,
    n_processes: int,
    land_array: np.ndarray,
    lambda_parameter: float,
    gen_mode: str,
    number_of_gens: int,
) -> np.ndarray:
    """
    Find connectivity for every cell of `condition`.

    Args:
        condition: 2D scalar condition matrix
        n_processes: Size of process pool to use (parallelise across pixels)
        land_array: 2D boolean array indicating land cells
        lambda_parameter: Average dispersal distance
        gen_mode: 'one-generation' or 'multi-generation'
        number_of_gens: Number of generations

    Returns:
        Connectivity matrix
    """

    connectivity = np.zeros(condition.shape)

    n_sectors: int = 8

    # set the different radii for the next bins according to a dart pattern
    # avail_ring_radii = np.array([1.5, 3.5, 5.5, 8.5, 12.5, 18, 25.7])
    avail_ring_radii = np.array(
        [
            1.5,
            3.5,
            5.5,
            8.5,
            12.5,
            18,
            25.7,
            36.7,
            52.1,
            73.8,
            104.4,
            147.6,
            200,
            300,
            500,
        ]
    )  # ; // number of cells

    ring_radii = np.array([])
    # TODO: could this be calculated from ring_radii, it's (mostly) just a diff
    # although... why isn't the first value 1.5?
    # avail_rad_width = np.array([1, 2, 2, 3, 4, 5.5, 7.7])
    # avail_rad_width = np.array([1, 2, 2, 2.5, 4, 5.5, 7.7, 11, 15.4, 21.7, 30.6, 43.2, 52.4, 100, 200]) #; // number of cells #!!!!!!!!!!!!!!!!!!!
    avail_rad_width = np.array(
        [1, 2, 2, 3, 4, 5.5, 7.7, 11, 15.4, 21.7, 30.6, 43.2, 52.4, 100, 200]
    )  # ; // number of cells #I changed entry 3 from 2.5 to 3 so that it matches the ring radii
    rad_width = np.array([])
    try:
        dartboard_r_min = -1 * lambda_parameter * np.log(0.01)
    except TypeError:
        breakpoint()
    for r in range(len(avail_ring_radii)):
        if avail_ring_radii[r] >= dartboard_r_min:
            dartboard_out_ring_r = avail_ring_radii[r]
            dartboard_radius = math.ceil(dartboard_out_ring_r)
            ring_radii = np.append(ring_radii, avail_ring_radii[r])
            rad_width = np.append(rad_width, avail_rad_width[r])
            break
        else:
            ring_radii = np.append(ring_radii, avail_ring_radii[r])
            rad_width = np.append(rad_width, avail_rad_width[r])
        if dartboard_r_min > avail_ring_radii[-1]:
            dartboard_out_ring_r = avail_ring_radii[-1]
            dartboard_radius = math.ceil(dartboard_out_ring_r)

    # TODO: 7 is still hardcoded in many places, should use this variable instead
    n_rings: int = len(ring_radii)

    dartboard_d = (
        2 * dartboard_radius
    ) + 1  # diameter of the dartboard is the dartboard radius +1 cell in the centre

    sector_index_by_cell = -np.ones((dartboard_d, dartboard_d), dtype=int)
    ring_index_by_cell = -np.ones((dartboard_d, dartboard_d), dtype=int)

    for row in range(dartboard_d):
        for col in range(dartboard_d):

            # sector and ring indicies
            sector, ring = find_ring_and_sector(row, col, ring_radii, dartboard_radius)

            if ring == len(ring_radii) + 1:
                break

            sector_index_by_cell[row, col] = sector
            ring_index_by_cell[row, col] = ring

    # minimum permeability
    f_min_perm = 0.5  # this corresponds to 1-m in the overleaf doc
    f_perm_scale: float = (
        0.5  # this corresponds to m in the overleaf doc, # should be 1 - f_min_perm
    )

    # f_lambda = np.array([0.02, 0.2, 2, 20])    # this corresponds to 1/alpha in the overleaf doc
    f_lambda = np.array([lambda_parameter])  # value we use for now instead

    land_cells = np.argwhere(land_array)

    if land_cells.size == 0:
        return connectivity

    land_cells = np.ascontiguousarray(land_cells, dtype=np.int64)
    condition = np.ascontiguousarray(condition)
    sector_index_by_cell = np.ascontiguousarray(sector_index_by_cell)
    ring_index_by_cell = np.ascontiguousarray(ring_index_by_cell)
    rad_width = np.ascontiguousarray(rad_width)
    f_lambda = np.ascontiguousarray(f_lambda)

    mode_flag = _gen_mode_flag(gen_mode, number_of_gens)

    values = connectivity_grid_kernel(
        condition,
        land_cells,
        sector_index_by_cell,
        ring_index_by_cell,
        n_sectors,
        n_rings,
        dartboard_radius,
        dartboard_d,
        f_lambda,
        mode_flag,
        number_of_gens,
        rad_width,
        f_perm_scale,
        f_min_perm,
    )

    connectivity[land_cells[:, 0], land_cells[:, 1]] = values

    return connectivity


def draw_matrix(
    arr: np.ndarray, title: str, filepath: str, cmap: str = "viridis"
) -> None:
    """
    Use imshow to draw matrix (2D array), coloured by value. Save to disk.

    Args:
        arr: Matrix to draw
        title: Title of plot
        filepath: Path to save plot to, including extension
        cmap: Matplotlib colourmap to shade values by
    """

    f, ax = plt.subplots()
    divider = make_axes_locatable(ax)
    img = ax.imshow(arr, cmap=cmap)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(img, cax=cax, orientation="vertical")
    ax.set_title(title)
    plt.savefig(filepath)
    plt.close(f)
    return


def landscape_connectivity(
    condition,
    n_processes,
    land_array,
    lambda_parameter,
    gen_mode,
    number_of_gens,
    draw_plots: bool = True,
):
    connectivity: np.ndarray = connectivity_of_grid(
        condition, n_processes, land_array, lambda_parameter, gen_mode, number_of_gens
    )
    result: np.ndarray = condition * connectivity
    landscape_connectivity = math.sqrt(np.sum(result))

    if draw_plots:
        draw_matrix(connectivity, "Connectivity", "connectivity.png", "plasma")
        draw_matrix(
            result,
            "Condition * Connectivity",
            "condition_x_connectivity.png",
            "cividis",
        )

    return landscape_connectivity


if __name__ == "__main__":

    # set seed for reproducibility
    np.random.seed(2)

    # how many processes to parallelise with, default to number of cpus available
    n_processes: int = os.cpu_count()

    # create matrix of randomly sampled [0, 1) values
    # value represents the habitat condition of each pixel
    i, j = 2600, 1400
    condition = np.random.rand(i, j)

    draw_matrix(condition, "Condition", "condition.png")

    start = time.time()

    connectivity: np.ndarray = connectivity_of_grid(
        condition,
        n_processes,
        land_array=np.ones_like(condition),
        lambda_parameter=5.0,
        gen_mode="one_generation",
        number_of_gens=1,
    )

    # final result is connectivity multiplied by condition
    # (so each cell's final value takes into account its own condition too)
    result: np.ndarray = condition * connectivity

    duration = time.time() - start

    print(f"Time elapsed: {duration:.1f} seconds")

    draw_matrix(connectivity, "Connectivity", "connectivity.png", "plasma")
    draw_matrix(
        result, "Condition * Connectivity", "condition_x_connectivity.png", "cividis"
    )

    landscape_connectivity: float = math.sqrt(np.sum(result))
    print(landscape_connectivity)
    """
    try:
        # with seed=2, (30, 40)
        # assert np.isclose(np.sum(result), 42963.06921981)

        # with seed=2, (80, 120)
        assert np.isclose(landscape_connectivity, 484199.381269)
    except AssertionError as error:
        print(f"Unexpected result: {landscape_connectivity:.6f}")
        raise error
    """
    print("Success!")
