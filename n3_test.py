import parameters as p
import numpy as np
import partition as z
import partition_N3 as zn3
import matplotlib.pyplot as plt
import g_matrix as gm
import time
import scipy.optimize as so

def func(param, sequence):
    N = len(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = False # FIXME
        N = N - 1 # ignoring the '-' part of sequence


    # free energy parameters for loop closure and base pair stacking
    g_loop = 1.0 # kcal/mol
    g_stack = param[3] # kcal/mol

    g_base_pair = gm.generator(sequence, param[0], param[1], param[2], N)
    
    if is_circular:
        return 0.
    else:
        return zn3.linear(g_base_pair, g_loop, g_stack, N)


part = []
part_n3 = []
part_time = []
partn3_time = []
num = 0

part_au = []
part_n3_au = []
part_gu = []
part_n3_gu = []
part_gc = []
part_n3_gc = []
part_stack = []
part_n3_stack = []

for sequence in p.training_sequences:
    M = len(sequence)

    is_circular = False
    if sequence[-1] == '-':
        is_circular = False # FIXME
        M = M - 1

    g_base_pair = gm.generator(sequence, 5.69, 6., 4.09, M)
    
    print sequence
    
    grad = so.approx_fprime([5.69, 6., 4.09, -7.09], func, 1e-9, sequence)
    part_au.append(grad[0])
    part_gu.append(grad[1])
    part_gc.append(grad[2])
    part_stack.append(grad[3])

    part_n3_au.append(zn3.linear_derivatives(g_base_pair, 1., -7.09, M, 5.69))
    part_n3_gu.append(zn3.linear_derivatives(g_base_pair, 1., -7.09, M, 6.))
    part_n3_gc.append(zn3.linear_derivatives(g_base_pair, 1., -7.09, M, 4.09))
    part_n3_stack.append(zn3.linear_derivatives(g_base_pair, 1., -7.09, M, -7.09))
    #start = time.clock()
    #part.append(z.linear(g_base_pair, 1.0, -7.09, M))
    #end = time.clock()
    #part_time.append(end-start)
    #start = time.clock()
    #part_n3.append(zn3.linear(g_base_pair, 1.0, -7.09, M))
    #end = time.clock()
    #partn3_time.append(end-start)

x = np.linspace(min([min(part_au), min(part_gu),min(part_gc),min(part_stack)])-.5, max([max(part_au), max(part_gu),max(part_gc),max(part_stack)]) +.5, 100)
print np.linalg.norm(np.array(part_n3_au + part_n3_gu + part_n3_gc + part_n3_stack) - np.array(part_au + part_gu + part_gc + part_stack))/np.sqrt(len(part_au + part_gu + part_gc + part_stack))
print 'AU '+str(np.linalg.norm(np.array(part_n3_au) - np.array(part_au)) / np.sqrt(len(part_au)))
print 'GU '+str(np.linalg.norm(np.array(part_n3_gu) - np.array(part_gu)) / np.sqrt(len(part_gu)))
print 'GC '+str(np.linalg.norm(np.array(part_n3_gc) - np.array(part_gc)) / np.sqrt(len(part_gc)))
print 'Stack '+str(np.linalg.norm(np.array(part_n3_stack) - np.array(part_stack)) / np.sqrt(len(part_stack)))

plt.scatter(part_au, part_n3_au, color='blue', label='$g_{AU}$')
plt.scatter(part_gu, part_n3_gu, color='red', label='$g_{GU}$')
plt.scatter(part_gc, part_n3_gc, color='green', label='$g_{GC}$')
plt.scatter(part_stack, part_n3_stack, color='magenta', label='$g_{stack}$')
plt.legend(loc=2)
plt.plot(x, x, color='black')
plt.xlim([min(x), max(x)])
plt.ylim([min(x), max(x)])
plt.xlabel('Numerical')
plt.ylabel('$N^3$')
plt.title('$N^3$ Partition Gradient Check (linear)')
plt.show()

#print 'N4 '+str(np.mean(np.array(part_time)))
#print 'N3 '+str(np.mean(np.array(partn3_time)))

#x = np.linspace(min(part)-.5, max(part)+.5, 100)
#print np.linalg.norm(np.array(part) - np.array(part_n3))/np.sqrt(len(part))


#plt.scatter(part, part_n3)
#plt.legend()
#plt.plot(x, x, color='black')
#plt.xlim([min(x), max(x)])
#plt.ylim([min(x), max(x)])
#plt.xlabel('$N^4$')
#plt.ylabel('$N^3$')
#plt.title('Partition Two Ways')
#plt.show()
