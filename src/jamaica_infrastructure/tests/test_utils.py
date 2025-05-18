import numpy as np
import pandas as pd

from jamaica_infrastructure.utils import numeric_only_dataframe, is_sole_value


class Test_numeric_only_dataframe:
    def test_numeric_only_dataframe(self):
        df = pd.DataFrame({"a": [-0, 2, 3], "b": [4, 5, 6], "c": [1.2, 3.4, np.nan]})
        assert numeric_only_dataframe(df)

        df = pd.DataFrame({"a": [1, 2, 3], "b": ["foo", "bar", "baz"]})
        assert not numeric_only_dataframe(df)

        df = pd.DataFrame({"a": [1, "2", 3]})
        assert not numeric_only_dataframe(df)


class Test_is_sole_value:
    def test_is_sole_value(self):
        assert is_sole_value(pd.Series([1, 1, 1]), 1)
        assert is_sole_value(pd.Series(["a", "a"]), "a")
        assert is_sole_value(pd.Series([np.nan, 1]), 1)
        assert not is_sole_value(pd.Series([1, 2, 2]), 2)
        assert not is_sole_value(pd.Series([1, 2, 2]), 2)
        assert not is_sole_value(pd.Series([1, 2, 2]), 4)
        assert not is_sole_value(pd.Series([np.nan, 1]), 2)
        assert not is_sole_value(pd.Series([np.nan, np.nan]), 1)
