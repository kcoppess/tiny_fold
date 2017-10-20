import pytest
import parameters as p
import log_partition as log_z
import partition as z
import g_matrix as gm


def partition_value_test(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return z.linear(g_base_pair, p.g_loop, p.energies[3], len(seq))

def log_partition_value_test(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return log_z.linear(g_base_pair, p.g_loop, p.energies[3], len(seq))

def almost_equal(x, y, epsilon ):
      return abs(x - y) < epsilon 

def test_partition():
    assert(partition_value_test("AAAAA") == 1)
    assert(almost_equal(partition_value_test("AAAAU"), 1.000012479163, 0.0000000001))
    assert(almost_equal(partition_value_test("GCCCC"), 1.000185778783, 0.0000000001))
    assert(almost_equal(partition_value_test("ACCCCGU"), 1.0021716139156, 0.0000000001))
    assert(almost_equal(partition_value_test("AACCCCGUU"), 1.02330274599, 0.0000000001))
