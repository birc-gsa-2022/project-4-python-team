# This directory will be checked with pytest. It will examine
# all files that start with test_*.py and run all functions with
# names that start with test_
from fm import bwt, compress, fm_index
import cProfile


def test_1984():
    assert 2 + 2 == 4


def test_correct_bwt_construction():
    x = "mississippi"
    cor_l = "ipssm$pissii"

    assert bwt(x) == cor_l


def test_compression():
    x = "mississippi"
    cor_com = "1i1p2s1m1$1p1i2s2i"

    assert compress(x) == cor_com


def profile():
    cProfile.run('fm_index("mississippi")')


test_1984()
test_correct_bwt_construction()
test_compression()
profile()
