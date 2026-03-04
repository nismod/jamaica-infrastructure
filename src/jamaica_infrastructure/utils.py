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
    row = df[(df["asset_gpkg"] == asset_gpkg) & (df["asset_layer"] == asset_layer)]
    if len(row) == 0:
        raise ValueError(
            f"No asset found for gpkg={asset_gpkg} and layer={asset_layer}"
        )
    elif len(row) > 1:
        raise ValueError(
            f"Multiple assets found for gpkg={asset_gpkg} and layer={asset_layer}"
        )
    return row.squeeze()


def parse_jic2005(c: str) -> tuple[str]:
    """Parse existing subsector codes to standard JIC 2005 two- or three-digit codes"""

    # special case 1
    if c == "500-2/3":
        return "50", "502", "503"

    # special case 2
    if c == "401/410":
        return "401", "41"

    # special case 3
    if c == "132.0":
        return ("132",)

    if c == "RES":
        return ()

    if c == "":
        return ()

    # length-two codes should have zero prefix (i.e. 20 means 020)
    if len(c) == 2:
        return (f"0{c}",)

    # zero-terminated codes should be two-digit (i.e. 650 means 65)
    if len(c) == 3:
        if c[2] == "0":
            return (c[:2],)
        else:
            return (c,)

    # six-digit codes should be split (i.e. 180190 means 180, 190 => 18, 19)
    if len(c) == 6:
        codes = []
        for code in parse_jic2005(c[:3]):
            codes.append(code)
        for code in parse_jic2005(c[3:]):
            codes.append(code)
        return tuple(codes)

    # hyphenate codes indicate a range (i.e. 011-2 means 011, 012)
    if "-" in c and "/" not in c:
        start, end_digit = c.split("-")
        prefix = start[:2]
        suffix_int = int(start[2])
        end_int = int(end_digit)

        if end_int == suffix_int:
            return parse_jic2005(start)

        if end_int < suffix_int:
            start = f"{prefix}{end_int}"
            tmp = suffix_int
            suffix_int = end_int
            end_int = tmp

        cs = []
        for c in parse_jic2005(start):
            cs.append(c)
        while True:
            suffix_int = suffix_int + 1
            cs.append(f"{prefix}{suffix_int}")
            if suffix_int == end_int:
                break
        return tuple(cs)

    raise ValueError(f"Unhandled JIC 2005 subsector code '{c}'")
