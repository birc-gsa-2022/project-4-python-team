# This directory will be checked with pytest. It will examine
# all files that start with test_*.py and run all functions with
# names that start with test_
from fm import bwt, compress


def test_1984():
    assert 2 + 2 == 4


def test_correct_bwt_construction():
    x = "mississippi$"
    cor_l = "ipssm$pissi"

    assert bwt(x) == cor_l


def test_compression():
    x = "mississippi$"
    cor_com = "1i1p2s1m1$1p1i2s1i"

    assert compress(x) == cor_com
