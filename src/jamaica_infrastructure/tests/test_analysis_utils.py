import numpy as np
import pandas as pd

from jamaica_infrastructure.analysis.utils import numeric_only_dataframe


class Test_numeric_only_dataframe:
    def test_numeric_only_dataframe(self):
        df = pd.DataFrame({"a": [-0, 2, 3], "b": [4, 5, 6], "c": [1.2, 3.4, np.nan]})
        assert numeric_only_dataframe(df)

        df = pd.DataFrame({"a": [1, 2, 3], "b": ["foo", "bar", "baz"]})
        assert not numeric_only_dataframe(df)

        df = pd.DataFrame({"a": [1, "2", 3]})
        assert not numeric_only_dataframe(df)