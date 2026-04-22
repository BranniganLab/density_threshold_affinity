import matplotlib
import pytest

@pytest.fixture(autouse=True, scope="session")
def _force_agg_backend():
    matplotlib.use("Agg", force=True)
    return True
