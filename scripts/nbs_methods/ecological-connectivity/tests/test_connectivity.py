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


@pytest.fixture(scope="module")
def default_params(sample_land: np.ndarray) -> dict:
    return {
        "n_processes": 1,
        "land_array": sample_land,
        "lambda_parameter": 0.1,
        "gen_mode": "one_generation",
        "number_of_gens": 1,
    }


#def test_gen_mode_flag_variants() -> None:
#    assert connectivity._gen_mode_flag("one_generation", 1) == 0
#    assert connectivity._gen_mode_flag("multi-generation", 3) == 1
#    assert connectivity._gen_mode_flag("multi-generation", 20) == 2
#
#
#def test_gen_mode_flag_invalid_value() -> None:
#    with pytest.raises(ValueError):
#        connectivity._gen_mode_flag("unsupported", 1)


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
    )

    assert value == pytest.approx(0.03572824440437287, rel=1e-6)


def test_landscape_connectivity_all_ones(default_params: dict) -> None:
    ones = np.ones((5, 5), dtype=np.float64)
    value = connectivity.landscape_connectivity(ones, **default_params)
    assert value == pytest.approx(0.08759331098811107, rel=1e-6)


def test_landscape_connectivity_all_zeros(default_params: dict) -> None:
    zeros = np.zeros((5, 5), dtype=np.float64)
    value = connectivity.landscape_connectivity(zeros, **default_params)
    assert value == pytest.approx(0.0, abs=1e-12)


def test_landscape_connectivity_half_landscape(default_params: dict) -> None:
    arr = np.zeros((5, 5), dtype=np.float64)
    arr[:, : arr.shape[1] // 2] = 1.0
    value = connectivity.landscape_connectivity(arr, **default_params)
    assert value == pytest.approx(0.048588026793122734, rel=1e-6)


def test_landscape_connectivity_checkerboard(default_params: dict) -> None:
    rows, cols = np.indices((5, 5))
    arr = (rows + cols) % 2
    arr = arr.astype(np.float64)
    value = connectivity.landscape_connectivity(arr, **default_params)
    assert value == pytest.approx(0.039494222107133735, rel=1e-6)


def test_landscape_connectivity_diagonal(default_params: dict) -> None:
    arr = np.eye(5, dtype=np.float64)
    value = connectivity.landscape_connectivity(arr, **default_params)
    assert value == pytest.approx(0.020373457120273282, rel=1e-6)


def test_connectivity_grid_all_ones(sample_land: np.ndarray) -> None:
    arr = np.ones((5, 5), dtype=np.float64)
    result = connectivity.connectivity_of_grid(
        arr,
        n_processes=1,
        land_array=sample_land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
    )
    expected = np.array(
        [
            [0.00018159971905, 0.000272399579, 0.000272399579, 0.000272399579, 0.00018159971905],
            [0.000272399579, 0.0004085993678623636, 0.0004085993678623636, 0.0004085993678623636, 0.000272399579],
            [0.000272399579, 0.0004085993678623636, 0.0004085993678623636, 0.0004085993678623636, 0.000272399579],
            [0.000272399579, 0.0004085993678623636, 0.0004085993678623636, 0.0004085993678623636, 0.000272399579],
            [0.00018159971905, 0.000272399579, 0.000272399579, 0.000272399579, 0.00018159971905],
        ]
    )
    np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-9)


def test_connectivity_grid_all_zeros(sample_land: np.ndarray) -> None:
    arr = np.zeros((5, 5), dtype=np.float64)
    result = connectivity.connectivity_of_grid(
        arr,
        n_processes=1,
        land_array=sample_land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
    )
    assert np.count_nonzero(result) == 0


def test_connectivity_grid_half_landscape(sample_land: np.ndarray) -> None:
    arr = np.zeros((5, 5), dtype=np.float64)
    arr[:, : arr.shape[1] // 2] = 1.0
    result = connectivity.connectivity_of_grid(
        arr,
        n_processes=1,
        land_array=sample_land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
    )
    expected = np.array(
        [
            [1.815997190499e-04, 1.815997190499e-04, 9.079985952497e-05, 0.0, 0.0],
            [2.723995785749e-04, 2.723995785749e-04, 1.361997892875e-04, 0.0, 0.0],
            [2.723995785749e-04, 2.723995785749e-04, 1.361997892875e-04, 0.0, 0.0],
            [2.723995785749e-04, 2.723995785749e-04, 1.361997892875e-04, 0.0, 0.0],
            [1.815997190499e-04, 1.815997190499e-04, 9.079985952497e-05, 0.0, 0.0],
        ]
    )
    np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-9)


def test_connectivity_grid_checkerboard(sample_land: np.ndarray) -> None:
    rows, cols = np.indices((5, 5))
    arr = ((rows + cols) % 2).astype(np.float64)
    result = connectivity.connectivity_of_grid(
        arr,
        n_processes=1,
        land_array=sample_land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
    )
    expected = np.array(
        [
            [9.079985952497e-05, 1.361997892875e-04, 1.361997892875e-04, 1.361997892875e-04, 9.079985952497e-05],
            [9.241945631728e-05, 1.378193860798e-04, 1.832193158423e-04, 1.378193860798e-04, 9.241945631728e-05],
            [9.241945631728e-05, 1.832193158423e-04, 1.378193860798e-04, 1.832193158423e-04, 9.241945631728e-05],
            [9.241945631728e-05, 1.378193860798e-04, 1.832193158423e-04, 1.378193860798e-04, 9.241945631728e-05],
            [4.701952655480e-05, 9.241945631728e-05, 9.241945631728e-05, 9.241945631728e-05, 4.701952655480e-05],
        ]
    )
    np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-9)


def test_connectivity_grid_diagonal(sample_land: np.ndarray) -> None:
    arr = np.eye(5, dtype=np.float64)
    result = connectivity.connectivity_of_grid(
        arr,
        n_processes=1,
        land_array=sample_land,
        lambda_parameter=0.1,
        gen_mode="one_generation",
        number_of_gens=1,
    )
    expected = np.array(
        [
            [9.079985952497e-05, 9.079985952497e-05, 4.539992976248e-05, 0.0, 0.0],
            [4.701952655480e-05, 9.241945631728e-05, 9.079985952497e-05, 4.539992976248e-05, 0.0],
            [4.539992976248e-05, 4.701952655480e-05, 9.241945631728e-05, 9.079985952497e-05, 4.539992976248e-05],
            [0.0, 4.539992976248e-05, 4.701952655480e-05, 9.241945631728e-05, 9.079985952497e-05],
            [0.0, 0.0, 4.539992976248e-05, 4.701952655480e-05, 4.701952655480e-05],
        ]
    )
    np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-9)
