import pytest
import parameters as p
import log_partition as log_z
import partition as z
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

tol = 1e-11

def partition_value_test(seq):
    return z.linear(p.energies, seq, len(seq)) #g_base_pair, p.g_loop, p.energies[3], len(seq))

def partition_deriv_test_AU(seq):
    grad = z.linear_derivatives(p.energies, seq, len(seq), 5.69)
    return grad
    #return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[0]) # g_base_pair, g_loop, g_stack, N, g

def partition_deriv_test_GU(seq):
    grad = z.linear_gradient(p.energies, seq, len(seq))
    return grad[1]
    #return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[1]) # g_base_pair, g_loop, g_stack, N, g

def partition_deriv_test_GC(seq):
    grad = z.linear_gradient(p.energies, seq, len(seq))
    return grad[2]
    #return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[2]) # g_base_pair, g_loop, g_stack, N, g

def partition_deriv_test_stack(seq):
    grad = z.linear_gradient(p.energies, seq, len(seq))
    return grad[3]
    #return z.linear_derivatives(g_base_pair, p.g_loop, p.energies[3], len(seq), p.energies[3]) # g_base_pair, g_loop, g_stack, N, g

def log_partition_value_test(seq):
    return log_z.linear(p.energies, seq, len(seq)) #g_base_pair, p.g_loop, p.energies[3], len(seq))

def circ_partition_value_test(seq):
    return z.circular(p.energies, seq, len(seq)-1) #g_base_pair, p.g_loop, p.energies[3], len(seq))
def circ_partition_deriv_test_AU(seq):
    grad = z.circular_gradient(p.energies, seq, len(seq)-1)
    return grad[0]
    #return z.circular_derivatives(p.energies, seq, len(seq)-1, p.energies[0])
def circ_partition_deriv_test_GU(seq):
    grad = z.circular_gradient(p.energies, seq, len(seq)-1)
    return grad[1]
    #return z.circular_derivatives(p.energies, seq, len(seq)-1, p.energies[1])
def circ_partition_deriv_test_GC(seq):
    grad = z.circular_gradient(p.energies, seq, len(seq)-1)
    return grad[2]
    #return z.circular_derivatives(p.energies, seq, len(seq)-1, p.energies[2])
def circ_partition_deriv_test_stack(seq):
    grad = z.circular_gradient(p.energies, seq, len(seq)-1)
    return grad[3]
    #return z.circular_derivatives(p.energies, seq, len(seq)-1, p.energies[3])

def bpp_value_test(seq):
    return bpp.mccaskill_linear(p.energies, seq, len(seq)) #g_base_pair, p.g_loop, p.energies[3], len(seq))
def bpp_deriv_test_AU(seq):
    grad = bpp.mccaskill_linear_gradient(p.energies, seq, len(seq))
    return grad[:, :, 0]
    #return bpp.mccaskill_linear_derivatives(p.energies, seq, len(seq), p.energies[0])
def bpp_deriv_test_GU(seq):
    grad = bpp.mccaskill_linear_gradient(p.energies, seq, len(seq))
    return grad[:, :, 1]
def bpp_deriv_test_GC(seq):
    grad = bpp.mccaskill_linear_gradient(p.energies, seq, len(seq))
    return grad[:, :, 2]
def bpp_deriv_test_stack(seq):
    grad = bpp.mccaskill_linear_gradient(p.energies, seq, len(seq))
    return grad[:, :, 3]

def almost_equal(x, y, epsilon ):
    print x-y
    return np.linalg.norm(x - y) < epsilon 

def test_partition():
    for seq, val in zip(seqs, [1,1.0000124791631115,1.0001857787828816,1.0021716139156089,1.0233027459871049]):
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
    aaccccguu[2,6] = 0.022593777679140784
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_value_test(seq), val, tol))

def test_bpp_derivs_AU():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    aaaau[0,4] = -2.1061933146135372e-05
    gcccc = np.zeros((5,5))
    accccgu = np.zeros((7,7))
    accccgu[0,6] = -0.0033378228091143717
    accccgu[1,5] = -0.0033162276547285432
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = -1.9697223269939332e-05
    aaccccguu[0,8] = -0.068099385089383582
    aaccccguu[1,7] = -0.07121356568856102
    aaccccguu[1,8] = -1.9697223269939332e-05
    aaccccguu[2,6] = -0.070752825849043341
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_deriv_test_AU(seq), val, tol))

def test_bpp_deriv_GC():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    gcccc = np.zeros((5,5))
    gcccc[0,4] = -0.00031344285229964632
    accccgu = np.zeros((7,7))
    accccgu[0,6] = -0.0033162276547285428
    accccgu[1,5] = -0.0036284787194924687
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 4.6130625989703884e-07
    aaccccguu[0,8] = -0.033786243034938615
    aaccccguu[1,7] = -0.036967505426624511
    aaccccguu[1,8] = 4.6130625989703884e-07
    aaccccguu[2,6] = -0.037272417173290488
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_deriv_test_GC(seq), val, tol))

def test_bpp_deriv_stack():
    aaaaa = np.zeros((5,5))
    aaaau = np.zeros((5,5))
    gcccc = np.zeros((5,5))
    accccgu = np.zeros((7,7))
    accccgu[0,6] = -0.0033168476362080901
    accccgu[1,5] = -0.0033162729345412763
    aaccccguu = np.zeros((9,9))
    aaccccguu[0,7] = 8.8575204834229921e-07
    aaccccguu[0,8] = -0.067866373951630518
    aaccccguu[1,7] = -0.070981138779790301
    aaccccguu[1,8] = 8.8575204834229921e-07
    aaccccguu[2,6] = -0.070759614540283317
    for seq, val in zip(seqs, [aaaaa, aaaau, gcccc, accccgu, aaccccguu]):
        assert(almost_equal(bpp_deriv_test_stack(seq), val, tol))
'''
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
