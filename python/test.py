from cmpb import *
import numpy as np
from scipy.sparse import coo_matrix
from ctypes import c_void_p, byref, POINTER, c_double, c_int64, \
c_char, c_char_p
from operator import add as op_add

c_int64_p = POINTER(c_int64)
c_double_p = POINTER(c_double)

def low_level():
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
    constrconeindices = np.array([ 0, 1 ]).astype(np.int64)
    # constrconeindices =  [ constrconeindices1 ]
    constrconelengths = np.array([ 2 ]).astype(np.int64)

    numvarcones = 1
    varconetypes = np.array([ MPBNONNEGCONE ]).astype(np.int64)
    varconeindices = np.array([ 0, 1, 2 ]).astype(np.int64)
    # varconeindices = [ varconeindices1 ]
    varconelengths = np.array([ 3 ]).astype(np.int64)
    # ---------------------------------------------------

    # required: setup the julia context
    lib.mpb_initialize();

    solver = c_void_p()
    model = c_void_p()

    lib.mpb_new_solver("ECOS","ECOSSolver()", byref(solver))
    lib.mpb_new_model(solver, byref(model))

    err = lib.mpb_loadproblem(model, nvar, nconstr,
    	ndarray_pointer(c),
    	ndarray_pointer(I),
    	ndarray_pointer(J),
    	ndarray_pointer(V),
    	nnz, ndarray_pointer(b),
        numconstrcones,
        ndarray_pointer(constrconetypes),
        ndarray_pointer(constrconeindices),
        ndarray_pointer(constrconelengths),
        numvarcones,
        ndarray_pointer(varconetypes),
        ndarray_pointer(varconeindices),
        ndarray_pointer(varconelengths));

    assert err == 0

    err = lib.mpb_optimize(model);

    assert err == 0

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
    # ----------------- problem data --------------------
    nvar = 3
    nconstr = 2
    nnz = 5

    c = np.array([-3, -2, -4])
    b = np.array([3, 2])


    I = [0,0,1,0,1]
    J = [0,1,1,2,2]
    V = [1.0,1.0,1.0,1.0,1.0]


    numconstrcones = 1;
    constrconetypes = [ MPBZEROCONE ]
    constrconeindices = [ 0, 1 ]
    # constrconeindices =  [ constrconeindices1 ]
    constrconelengths = [ 2 ]

    numvarcones = 1
    varconetypes = [ MPBNONNEGCONE ]
    varconeindices = [ 0, 1, 2 ]
    # varconeindices = [ varconeindices1 ]
    varconelengths = [ 3 ]
    # ---------------------------------------------------
    constrcones = MPBCones(constrconetypes, constrconelengths, constrconeindices)
    varcones = MPBCones(varconetypes, varconelengths, varconeindices)
    A = coo_matrix((V, (I, J)), shape=(nconstr, nvar))

    problem = MPBModel("ECOS", "ECOSSolver()", c, A, b,
    constrcones, varcones)
    problem.optimize()
    print problem.status()
    assert abs(problem.getproperty("objval") - (-11)) < 1e-3

    sol = problem.getsolution()
    assert abs(sol[0]-1.0) < 1e-3
    assert abs(sol[1]-0.0) < 1e-3
    assert abs(sol[2]-2.0) < 1e-3
    del problem

def int_constr():
    # ----------------- problem data --------------------
    c = np.array([0, -2, -1]);
    I = [0,1,2,3];
    J = [0,0,1,2];
    V = [1.,-1.,-1.,-1.];
    b = np.array([1,0,0,0]);
    nvar = 3;
    nconstr = 4;
    nnz = 4;
    numconstrcones = 2;
    constrconetypes = [ MPBZEROCONE, MPBSOC ];
    constrconeindices = [ 0, 1, 2, 3 ];
    constrconelengths = [ 1, 3 ];
    numvarcones = 1;
    varconetypes = [ MPBFREECONE ]
    varconeindices = [ 0, 1, 2 ]
    # varconeindices = [ varconeindices1 ]
    varconelengths = [ 3 ]
    vartypes = [ MPBCONTVAR, MPBBINVAR, MPBBINVAR ]

    # ---------------------------------------------------
    constrcones = MPBCones(constrconetypes, constrconelengths, constrconeindices)
    varcones = MPBCones(varconetypes, varconelengths, varconeindices)
    A = coo_matrix((V, (I, J)), shape=(nconstr, nvar))

    problem = MPBModel("Pajarito, GLPKMathProgInterface, ECOS","PajaritoSolver(verbose=1,mip_solver=GLPKSolverMIP(),cont_solver=ECOSSolver(verbose=0))", c, A, b,
    constrcones, varcones, vartypes)
    problem.optimize()
    print problem.status()
    assert abs(problem.getproperty("objval") - (-2)) < 1e-3

    sol = problem.getsolution()
    assert abs(sol[0]-1.0) < 1e-3
    assert abs(sol[1]-1.0) < 1e-3
    assert abs(sol[2]-0.0) < 1e-3

if __name__ == "__main__":
    low_level()
    high_level()
    int_constr()
