from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import connectivity  # noqa: E402


@pytest.fixture(scope="module")
def sample_condition() -> np.ndarray:
    return np.array(
        [
            [0.2, 0.5, 0.5, 0.5, 0.2],
            [0.5, 0.8, 0.9, 0.8, 0.5],
            [0.5, 0.9, 1.0, 0.9, 0.5],
            [0.5, 0.8, 0.9, 0.8, 0.5],
            [0.2, 0.5, 0.5, 0.5, 0.2],
        ],
        dtype=np.float64,
    )


@pytest.fixture(scope="module")
def sample_land(sample_condition: np.ndarray) -> np.ndarray:
    return np.ones_like(sample_condition, dtype=bool)


def test_gen_mode_flag_variants() -> None:
    assert connectivity._gen_mode_flag("one_generation", 1) == 0
    assert connectivity._gen_mode_flag("multi-generation", 3) == 1
    assert connectivity._gen_mode_flag("multi-generation", 20) == 2


def test_gen_mode_flag_invalid_value() -> None:
    with pytest.raises(ValueError):
        connectivity._gen_mode_flag("unsupported", 1)


def test_identify_sector_cardinal_quadrants() -> None:
    assert connectivity.identify_sector(-5, 0) == 0  # north
    assert connectivity.identify_sector(5, 0) == 4  # south
    assert connectivity.identify_sector(0, 5) == 2  # east
    assert connectivity.identify_sector(3, -3) == 5  # south-west


def test_connectivity_of_grid_returns_zero_when_no_land() -> None:
    condition = np.ones((3, 3), dtype=np.float64)
    land = np.zeros_like(condition, dtype=bool)
    result = connectivity.connectivity_of_grid(
        condition,
        n_processes=1,
        land_array=land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
    )
    assert np.count_nonzero(result) == 0


def test_connectivity_of_grid_matches_known_values(
    sample_condition: np.ndarray, sample_land: np.ndarray
) -> None:
    result = connectivity.connectivity_of_grid(
        sample_condition,
        n_processes=1,
        land_array=sample_land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
    )

    expected = np.array(
        [
            [0.000013587423109176614, 0.00003853635337806715, 0.00005048106870085347, 0.00003853635337806715, 0.000013587423109176614],
            [0.000037972694320753494, 0.00010319481817372062, 0.00013009414790730698, 0.00010319481817372062, 0.000037972694320753494],
            [0.00005048106870085347, 0.0001423701247790868, 0.00018699247158477777, 0.00014237012477908678, 0.00005048106870085347],
            [0.00003853635337806715, 0.00010709607597972253, 0.0001413701932608621, 0.00010709607597972253, 0.00003853635337806715],
            [0.000013023764051862964, 0.00003284595814234525, 0.00003641595439935369, 0.000032845958142345244, 0.000013023764051862964],
        ],
        dtype=np.float64,
    )

    np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-9)


def test_landscape_connectivity_matches_expected(
    sample_condition: np.ndarray, sample_land: np.ndarray
) -> None:
    value = connectivity.landscape_connectivity(
        sample_condition,
        n_processes=1,
        land_array=sample_land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
        draw_plots=False,
    )

    assert value == pytest.approx(0.03572824440437287, rel=1e-6)
