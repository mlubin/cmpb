from cmpb import *
import numpy as np
from scipy.sparse import coo_matrix
from ctypes import c_void_p, byref, POINTER, c_double, c_int64, \
c_char, c_char_p
from operator import add as op_add

c_int64_p = POINTER(c_int64)
c_double_p = POINTER(c_double)

# ----------------- problem data --------------------
nvar = 3
nconstr = 2
nnz = 5

c = np.float64(np.array([-3, -2, -4]))
b = np.float64(np.array([3, 2]))


I = [0,0,1,0,1]
J = [0,1,1,2,2]
V = [1.0,1.0,1.0,1.0,1.0]


numconstrcones = 1;
constrconetypes = [ MPBZEROCONE ]
constrconeindices1 = [ 0, 1 ]
constrconeindices =  [ constrconeindices1 ]
constrconelengths = [ 2 ]

numvarcones = 1
varconetypes = [ MPBNONNEGCONE ]
varconeindices1 = [ 0, 1, 2 ]
varconeindices = [ varconeindices1 ]
varconelengths = [ 3 ]
# ---------------------------------------------------



def low_level():
    # required: setup the julia context
    lib.mpb_initialize();

    solver = c_void_p()
    model = c_void_p()

    lib.mpb_new_solver("ECOS","ECOSSolver()", byref(solver))
    lib.mpb_new_model(solver, byref(model))

    constr_idx_arr = np.zeros(numconstrcones, dtype=c_int64_p)
    for i, idx_list in enumerate(constrconeindices):
        constr_idx_arr[i] = ndarray_pointer(np.array(idx_list).astype(np.int64))
    constrconeindices_ptr = ndarray_pointer(constr_idx_arr)

    var_idx_arr = np.zeros(numvarcones, dtype=c_int64_p)
    for i, idx_list in enumerate(varconeindices):
        var_idx_arr[i] = ndarray_pointer(np.array(idx_list).astype(np.int64))
    varconeindices_ptr = ndarray_pointer(var_idx_arr)


    err = lib.mpb_loadproblem(model, nvar, nconstr,
    	ndarray_pointer(c),
    	ndarray_pointer(np.array(I).astype(np.int64)),
    	ndarray_pointer(np.array(J).astype(np.int64)),
    	ndarray_pointer(np.array(V).astype(np.float64)),
    	nnz, ndarray_pointer(b),
        numconstrcones,
        ndarray_pointer(np.array(constrconetypes).astype(np.int64)),
        ndarray_pointer(constr_idx_arr),
        ndarray_pointer(np.array(constrconelengths).astype(np.int64)),
        numvarcones,
        ndarray_pointer(np.array(varconetypes).astype(np.int64)),
        ndarray_pointer(var_idx_arr),
        ndarray_pointer(np.array(varconelengths).astype(np.int64)));

    assert err != 0

    err = lib.mpb_optimize(model);

    assert err != 0

    status = np.ndarray(20, c_char)
    ret = lib.mpb_status(model, status.ctypes.data_as(c_char_p), 20);
    print "STATUS: ", reduce(op_add, status);

    objval = np.zeros(1, dtype=np.float64);
    lib.mpb_getobjval(model, ndarray_pointer(objval))
    assert abs(objval[0] - (-11)) < 1e-3

    sol = np.zeros(nvar, dtype=np.float64)
    lib.mpb_getsolution(model, ndarray_pointer(sol))

    assert abs(sol[0]-1.0) < 1e-3
    assert abs(sol[1]-0.0) < 1e-3
    assert abs(sol[2]-2.0) < 1e-3

    lib.mpb_free_model(model)
    lib.mpb_free_solver(solver)

    lib.mpb_atexit(0)
    return 0

def high_level():
    constrcones = MPBCone(constrconetypes, constrconeindices)
    varcones = MPBCone(varconetypes, varconeindices)
    A = coo_matrix((V, (I, J)), shape=(nconstr, nvar))

    problem = MPBModel("ECOS", "ECOSSolver()", c, A, b,
    constrcones, varcones)
    problem.optimize()
    print problem.status()
    assert abs(problem.get("objval") - (-11)) < 1e-3

    sol = problem.getsolution()
    assert abs(sol[0]-1.0) < 1e-3
    assert abs(sol[1]-0.0) < 1e-3
    assert abs(sol[2]-2.0) < 1e-3

if __name__ == "__main__":
    low_level()
    high_level()






