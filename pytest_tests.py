import pytest
import parameters as p
import log_partition as log_z
import partition as z
import g_matrix as gm
import base_pair_probabilities as bpp
import numpy as np

seqs = [ "AAAAA",
         "AAAAU",
         "GCCCC",
         "ACCCCGU",
         "AACCCCGUU"
]
circs = ["AAAAA-",
         "AAAAU-",
         "ACCCUCCC-",
         "AAACCCUUUCCC-"
]

tol = 1e-8

def partition_value_test(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return z.linear(g_base_pair, p.g_loop, p.energies[3], len(seq))

def partition_deriv_test_AU(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[0]) # g_base_pair, g_loop, g_stack, N, g

def partition_deriv_test_GU(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[1]) # g_base_pair, g_loop, g_stack, N, g

def partition_deriv_test_GC(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[2]) # g_base_pair, g_loop, g_stack, N, g

def partition_deriv_test_stack(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[3]) # g_base_pair, g_loop, g_stack, N, g

def log_partition_value_test(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return log_z.linear(g_base_pair, p.g_loop, p.energies[3], len(seq))

def circ_partition_value_test(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq)-1)
    return z.circular(g_base_pair, p.g_loop, p.energies[3], len(seq)-1)
def circ_partition_deriv_test_AU(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq)-1)
    return z.circular_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq)-1, p.energies[0])
def circ_partition_deriv_test_GU(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq)-1)
    return z.circular_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq)-1, p.energies[1])
def circ_partition_deriv_test_GC(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq)-1)
    return z.circular_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq)-1, p.energies[2])
def circ_partition_deriv_test_stack(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq)-1)
    return z.circular_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq)-1, p.energies[3])

def bpp_value_test(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return bpp.mccaskill_linear(g_base_pair, p.g_loop, p.energies[3], len(seq))
def bpp_deriv_test_AU(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return bpp.mccaskill_linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[0])
def bpp_deriv_test_GU(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return bpp.mccaskill_linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[1])
def bpp_deriv_test_GC(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return bpp.mccaskill_linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[2])
def bpp_deriv_test_stack(seq):
    g_base_pair = gm.generator(seq, p.energies[0], p.energies[1], p.energies[2], len(seq))
    return bpp.mccaskill_linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[3])

def almost_equal(x, y, epsilon ):
    print x-y
    return np.linalg.norm(x - y) < epsilon 

def test_partition():
    for seq, val in zip(seqs, [1,1.000012479163,1.000185778783,1.0021716139156,1.02330274599]):
        assert(almost_equal(partition_value_test(seq), val, tol))

def test_log_partition():
    assert(True)

def test_derivs_AU():
    for seq, val in zip(seqs, [0,-2.10624588E-05,0,-0.00335171279,-0.0746191942895]):
        assert(almost_equal(partition_deriv_test_AU(seq), val, tol))

def test_derivs_GU():
    for seq, val in zip(seqs, [0,0,0,0,0]):
        assert(almost_equal(partition_deriv_test_GU(seq), val, tol))

def test_derivs_GC():
    for seq, val in zip(seqs, [0,0,-0.000313559325,-0.00364420966,-0.039022635596]):
        assert(almost_equal(partition_deriv_test_GC(seq), val, tol))

def test_derivs_stack():
    for seq, val in zip(seqs, [0,0,0,-0.003330650339,-0.074311205720273069]):
        assert(almost_equal(partition_deriv_test_stack(seq), val, tol))

def test_circ_partition():
    for seq, val in zip(circs, [1, 1, 1.0000023076971, 1.0003791934327]):
        assert(almost_equal(circ_partition_value_test(seq), val, tol))
def test_circ_derivs_AU():
    for seq, val in zip(circs, [0, 0, -3.8949547054E-06, -0.0016844202181]):
        assert(almost_equal(circ_partition_deriv_test_AU(seq), val, tol))
def test_circ_derivs_GU():
    for seq, val in zip(circs, [0,0,0,0]):
        assert(almost_equal(circ_partition_deriv_test_GU(seq), val, tol))
def test_circ_derivs_GC():
    for seq, val in zip(circs, [0,0,0,0]):
        assert(almost_equal(circ_partition_deriv_test_GC(seq), val, tol))
def test_circ_derivs_stack():
    for seq, val in zip(circs, [0,0,0,-0.001044413431]):
        assert(almost_equal(circ_partition_deriv_test_stack(seq), val, tol))

def test_bpp():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    aaaau[0,4] = 1.2479007383961551e-05
    gcccc = np.zeros((5,5))
    gcccc[0,4] = 0.00018574427553608532
    accccgu = np.zeros((7,7))
    accccgu[0,6] = 0.0019815320102396245
    accccgu[1,5] = 0.0021544561056377051
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 1.2197252009936157e-05
    aaccccguu[0,8] = 0.020625537448534667
    aaccccguu[1,7] = 0.022553953705591098
    aaccccguu[1,8] = 1.2197252009936157e-05
    aaccccguu[2,6] = 0.022593773148001303
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_value_test(seq), val, tol))
'''
def test_bpp_derivs_AU():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    aaaau[0,4] = 1.2479007383961551e-05
    gcccc = np.zeros((5,5))
    gcccc[0,4] = 0.00018574427553608532
    accccgu = np.zeros((7,7))
    accccgu[0,6] = 0.0019815320102396245
    accccgu[1,5] = 0.0021544561056377051
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 1.2197252009936157e-05
    aaccccguu[0,8] = 0.020625537448534667
    aaccccguu[1,7] = 0.022553953705591098
    aaccccguu[1,8] = 1.2197252009936157e-05
    aaccccguu[2,6] = 0.022593773148001303
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_deriv_test_AU(seq), val, tol))
def test_bpp():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    aaaau[0,4] = 1.2479007383961551e-05
    gcccc = np.zeros((5,5))
    gcccc[0,4] = 0.00018574427553608532
    accccgu = np.zeros((7,7))
    accccgu[0,6] = 0.0019815320102396245
    accccgu[1,5] = 0.0021544561056377051
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 1.2197252009936157e-05
    aaccccguu[0,8] = 0.020625537448534667
    aaccccguu[1,7] = 0.022553953705591098
    aaccccguu[1,8] = 1.2197252009936157e-05
    aaccccguu[2,6] = 0.022593773148001303
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_value_test(seq), val, tol))
def test_bpp():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    aaaau[0,4] = 1.2479007383961551e-05
    gcccc = np.zeros((5,5))
    gcccc[0,4] = 0.00018574427553608532
    accccgu = np.zeros((7,7))
    accccgu[0,6] = 0.0019815320102396245
    accccgu[1,5] = 0.0021544561056377051
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 1.2197252009936157e-05
    aaccccguu[0,8] = 0.020625537448534667
    aaccccguu[1,7] = 0.022553953705591098
    aaccccguu[1,8] = 1.2197252009936157e-05
    aaccccguu[2,6] = 0.022593773148001303
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_value_test(seq), val, tol))
def test_bpp():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    aaaau[0,4] = 1.2479007383961551e-05
    gcccc = np.zeros((5,5))
    gcccc[0,4] = 0.00018574427553608532
    accccgu = np.zeros((7,7))
    accccgu[0,6] = 0.0019815320102396245
    accccgu[1,5] = 0.0021544561056377051
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 1.2197252009936157e-05
    aaccccguu[0,8] = 0.020625537448534667
    aaccccguu[1,7] = 0.022553953705591098
    aaccccguu[1,8] = 1.2197252009936157e-05
    aaccccguu[2,6] = 0.022593773148001303
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_value_test(seq), val, tol))
'''
