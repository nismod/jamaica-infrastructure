from typing import Any

import pandas as pd
from scipy import integrate
import numpy as np
from tqdm import tqdm

tqdm.pandas()


def is_sole_value(series: pd.Series, value: Any) -> bool:
    """Check the non-null entries of a pandas Series all match `value`."""
    return set(series.dropna().unique()) == {value}


def numeric_only_dataframe(df: pd.DataFrame) -> bool:
    """Does a DataFrame contain only numeric columns?"""
    return df.dtypes.map(pd.api.types.is_numeric_dtype).all()


def get_asset(network_csv_path: str, asset_gpkg: str, asset_layer: str) -> pd.Series:
    """
    Get the path of an asset by its gpkg and layer strings.
    """
    df = pd.read_csv(network_csv_path)
    row = df[(df['asset_gpkg'] == asset_gpkg) & (df['asset_layer'] == asset_layer)]
    if len(row) == 0:
        raise ValueError(f"No asset found for gpkg={asset_gpkg} and layer={asset_layer}")
    elif len(row) > 1:
        raise ValueError(f"Multiple assets found for gpkg={asset_gpkg} and layer={asset_layer}")
    return row.squeeze()


def risks(
    dataframe,
    index_columns,
    probabilities,
    expected_risk_column,
    flood_protection_period=0,
    flood_protection_name=None,
):
    """
    Organise the dataframe to pivot with respect to index columns
    Find the expected risks
    """
    if flood_protection_name is None and flood_protection_period == 0:
        # When there is no flood protection at all
        expected_risk_column = f"{expected_risk_column}_undefended"
        probability_columns = [str(p) for p in probabilities]

    elif flood_protection_period > 0:
        if flood_protection_name is None:
            expected_risk_column = (
                f"{expected_risk_column}_{flood_protection_period}_year_protection"
            )
        else:
            expected_risk_column = f"{expected_risk_column}_{flood_protection_name}"

        probabilities = [
            pr for pr in probabilities if pr <= 1.0 / flood_protection_period
        ]
        probability_columns = [str(p) for p in probabilities]
    else:
        # When there is no flood protection at all
        expected_risk_column = f"{expected_risk_column}_{flood_protection_name}"
        probability_columns = [str(p) for p in probabilities]

    dataframe.columns = dataframe.columns.astype(str)
    dataframe[expected_risk_column] = list(
        integrate.trapezoid(
            dataframe[probability_columns].to_numpy(),
            np.array([probabilities * len(dataframe.index)]).reshape(
                dataframe[probability_columns].shape
            ),
        )
    )

    return dataframe[index_columns + [expected_risk_column]]
