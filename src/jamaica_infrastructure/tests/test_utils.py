import numpy as np
import pandas as pd

from jamaica_infrastructure.utils import (
    numeric_only_dataframe,
    is_sole_value,
    parse_jic2005,
)


class TestNumericOnlyDataframe:
    def test_numeric_only_dataframe(self):
        df = pd.DataFrame({"a": [-0, 2, 3], "b": [4, 5, 6], "c": [1.2, 3.4, np.nan]})
        assert numeric_only_dataframe(df)

        df = pd.DataFrame({"a": [1, 2, 3], "b": ["foo", "bar", "baz"]})
        assert not numeric_only_dataframe(df)

        df = pd.DataFrame({"a": [1, "2", 3]})
        assert not numeric_only_dataframe(df)


class TestIsSoleValue:
    def test_is_sole_value(self):
        assert is_sole_value(pd.Series([1, 1, 1]), 1)
        assert is_sole_value(pd.Series(["a", "a"]), "a")
        assert is_sole_value(pd.Series([np.nan, 1]), 1)
        assert is_sole_value(pd.Series([np.nan, 1]), 1)
        assert is_sole_value(pd.Series([np.nan, np.nan]), np.nan)

        assert not is_sole_value(pd.Series([1, 2, 2]), 2)
        assert not is_sole_value(pd.Series([1, 2, 2]), 2)
        assert not is_sole_value(pd.Series([1, 2, 2]), 4)
        assert not is_sole_value(pd.Series([np.nan, 1]), 2)
        assert not is_sole_value(pd.Series([np.nan, np.nan]), 1)
        assert not is_sole_value(pd.Series([np.nan, "a", 1]), "a")
        assert not is_sole_value(pd.Series([np.nan, "a", 1]), np.nan)


class TestParseJIC2005:
    def test_special_cases(self):
        code = "500-2/3"
        expected = "50", "502", "503"
        assert parse_jic2005(code) == expected

        code = "401/410"
        expected = "401", "41"
        assert parse_jic2005(code) == expected

        code = "132.0"
        expected = ("132",)
        assert parse_jic2005(code) == expected

        code = "RES"
        expected = ()
        assert parse_jic2005(code) == expected

    def test_length_two_codes_should_have_zero_prefix(self):
        code = "20"
        expected = ("020",)
        assert parse_jic2005(code) == expected

    def test_zero_terminated_codes_should_be_two_digit(self):
        code = "650"
        expected = ("65",)
        assert parse_jic2005(code) == expected

    def test_six_digit_codes_should_be_split(self):
        code = "180190"
        expected = (
            "18",
            "19",
        )
        assert parse_jic2005(code) == expected

    def test_hyphenate_codes_indicate_a_range(self):
        code = "011-2"
        expected = (
            "011",
            "012",
        )
        assert parse_jic2005(code) == expected

        code = "012-1"
        expected = (
            "011",
            "012",
        )
        assert parse_jic2005(code) == expected
