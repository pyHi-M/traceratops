import sys
import types

import numpy as np

from traceratops.plot_matrix_comparison import calculates_ensemble_matrices

sys.modules.setdefault(
    "seaborn", types.SimpleNamespace(set=lambda *args, **kwargs: None)
)


def test_calculates_ensemble_matrices_returns_2d_median_matrix():
    matrix = np.array(
        [
            [[np.nan, np.nan], [1.0, 3.0]],
            [[1.0, 3.0], [np.nan, np.nan]],
        ]
    )

    ensemble_matrices = calculates_ensemble_matrices(
        [matrix], mode="median", max_distance=np.inf
    )

    assert len(ensemble_matrices) == 1
    assert ensemble_matrices[0].shape == (2, 2)
    np.testing.assert_allclose(
        ensemble_matrices[0],
        [[np.nan, 2.0], [2.0, np.nan]],
        equal_nan=True,
    )
