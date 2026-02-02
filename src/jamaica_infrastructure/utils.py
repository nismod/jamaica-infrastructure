from typing import Any

import numpy as np
import pandas as pd


def is_sole_value(series: pd.Series, value: Any) -> bool:
    """
    Check the non-null entries of a pandas Series all match `value`. If `value`
    is null, check the Series is entirely null.
    """
    if isinstance(value, float) and np.isnan(value):
        return all(series.isna())
    else:
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
