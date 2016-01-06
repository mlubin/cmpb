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


I = np.array([0,0,1,0,1]).astype(np.int64)
J = np.array([0,1,1,2,2]).astype(np.int64)
V = np.array([1.0,1.0,1.0,1.0,1.0]).astype(np.float64)


numconstrcones = 1;
constrconetypes = np.array([ MPBZEROCONE ]).astype(np.int64)
constrconeindices1 = np.array([ 0, 1 ]).astype(np.int64)
constrconeindices =  [ constrconeindices1 ]
constrconelengths = np.array([ 2 ]).astype(np.int64)

numvarcones = 1
varconetypes = np.array([ MPBNONNEGCONE ]).astype(np.int64)
varconeindices1 = np.array([ 0, 1, 2 ]).astype(np.int64)
varconeindices = [ varconeindices1 ]
varconelengths = np.array([ 3 ]).astype(np.int64)
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
        tmp = idx_list.ctypes.data
        constr_idx_arr[i] = tmp
        print tmp
        print constr_idx_arr.shape
        # import pdb; pdb.set_trace()
        # constr_idx_arr[i] = ndarray_pointer(np.array(idx_list).astype(np.int64))
    constrconeindices_ptr = constr_idx_arr.ctypes.data_as(c_int64_pp)

    var_idx_arr = np.zeros(numvarcones, dtype=c_int64_p)
    for i, idx_list in enumerate(varconeindices):
        var_idx_arr[i] = idx_list.ctypes.data
    varconeindices_ptr = var_idx_arr.ctypes.data_as(c_int64_pp)

    # I_arr = np.array(I).astype(np.int64)
    # J_arr = np.array(J).astype(np.int64)
    # V_arr = np.array(V).astype(np.float64)
    err = lib.mpb_loadproblem(model, nvar, nconstr,
    	ndarray_pointer(c),
    	ndarray_pointer(I),
    	ndarray_pointer(J),
    	ndarray_pointer(V),
    	nnz, ndarray_pointer(b),
        numconstrcones,
        ndarray_pointer(constrconetypes),
        constrconeindices_ptr,
        ndarray_pointer(constrconelengths),
        numvarcones,
        ndarray_pointer(varconetypes),
        varconeindices_ptr,
        ndarray_pointer(varconelengths));

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
    constrcones = MPBCones(constrconetypes, constrconeindices)
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






