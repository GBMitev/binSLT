import pytest
from binslt.utils.parse_files import get_states_columns

def test_get_states_columns():
    expected_output = ["NN","E","gns","J","tau","e/f","Manifold","v","Lambda","Sigma","Omega"]
    result_output = get_states_columns()

    assert expected_output == result_output
