from ctypes import CDLL, POINTER, c_char, c_char_p, c_int, c_int64, \
    c_double, c_void_p, byref
from numpy import int64, float64, ndarray, array, void, hstack
from os import uname, path
from site import getsitepackages
from scipy.sparse import coo_matrix
from operator import add as op_add
from yaml import safe_load 

# int64_t *, double * pointers
c_int64_p = POINTER(c_int64)
c_int64_pp = POINTER(c_int64)
c_double_p = POINTER(c_double)

# np.ndarray --> c pointer
# valid types: char *, int64_t *, double *, OR exception)
def ndarray_pointer(array):
    if not isinstance(array, ndarray):
        TypeError("input must be a numpy ndarray")
    if not array.dtype in (int64, float64, c_char):
        AttributeError("input array must have int64 or float64 elements")
    if array.dtype == int64:
        return array.ctypes.data_as(c_int64_p)
    elif array.dtype == c_char:
        return array.ctypes.data_as(c_char_p)
    else:
        return array.ctypes.data_as(c_double_p)


# ------------- #
# MPB constants #
# ------------- #
MPBFREECONE = 0
MPBZEROCONE = 1
MPBNONNEGCONE = 2
MPBNONPOSCONE = 3
MPBSOC = 4
MPBSOCROTATED = 5
MPBSDPCONE = 6
MPBEXPPRIMAL = 7
MPBEXPDUAL = 8

STATUS_BUFFER_LENGTH = 100

# ----------------- #
# find and load lib #
# ----------------- #
CMPB_HOME = path.abspath(path.join(path.dirname(__file__), ".."))
LIBNAME = "libcmpb"
EXT = ".dylib" if uname()[0] == "Darwin" else ".so" # no Windows support
localbuild = CMPB_HOME
sitepath = getsitepackages()[0]

if path.exists(path.join(sitepath, LIBNAME + EXT)):
    libpath = path.join(sitepath, LIBNAME + EXT)
else:
    libpath = path.join(localbuild, LIBNAME + EXT)

try:
    lib = CDLL(libpath)
except:
    print("libcmpb not found at {}".format(libpath))
    raise


# ------------------ #
# MPB public methods #
# ------------------ #

# define arguments
lib.mpb_initialize.argtypes = []
lib.mpb_numvar.argtypes = [c_void_p, c_int64_p]
lib.mpb_numconstr.argtypes = [c_void_p, c_int64_p]
lib.mpb_getobjval.argtyps = [c_void_p, c_double_p]
lib.mpb_getobjbound.argtypes = [c_void_p, c_double_p]
lib.mpb_getobjgap.argtypes = [c_void_p, c_double_p]
lib.mpb_getsolvetime.argtypes = [c_void_p, c_double_p]
lib.mpb_new_solver.argtypes = [c_char_p, c_char_p, POINTER(c_void_p)]
lib.mpb_free_solver.argtypes = [c_void_p]
lib.mpb_new_model.argtypes = [c_void_p, POINTER(c_void_p)]
lib.mpb_free_model.argtypes = [c_void_p]
lib.mpb_atexit.argtypes = [c_int]
lib.mpb_loadproblem.argtypes = [c_void_p, c_int64, c_int64, c_double_p,
    c_int64_p, c_int64_p, c_double_p, c_int64, c_double_p,
    c_int64, c_int64_p, c_int64_p, c_int64_p,
    c_int64, c_int64_p, c_int64_p, c_int64_p]
lib.mpb_getsolution.argtypes = [c_void_p, c_double_p]
lib.mpb_getdual.argtypes = [c_void_p, c_double_p]
lib.mpb_status.argtypes = [c_void_p, c_char_p, c_int64]

# define return types
lib.mpb_initialize.restype = c_int
lib.mpb_numvar.restype = c_int
lib.mpb_numconstr.restype = c_int
lib.mpb_getobjval.restype = c_int
lib.mpb_getobjbound.restype = c_int
lib.mpb_getobjgap.restype = c_int
lib.mpb_getsolvetime.restype = c_int
lib.mpb_new_solver.restype = c_int
lib.mpb_free_solver.restype = c_int
lib.mpb_new_model.restype = c_int
lib.mpb_free_model.restype = c_int
lib.mpb_atexit.restype = None
lib.mpb_loadproblem.restype = c_int
lib.mpb_getsolution.restype = c_int
lib.mpb_getdual.restype = c_int
lib.mpb_optimize.restype = c_int
lib.mpb_status.restype = c_int

# --------------- #
# Python bindings #
# --------------- #

solver_file = open(path.join(CMPB_HOME, "solvers.yml"), 'r')
MPB_solverlist = safe_load(solver_file)
solver_file.close()

def MPB_CHECKERR(err):
    if err != 0:
        msg = "Error occurred in call to libcmpb"
        Warning(msg)
        # RuntimeError(msg)

'''
create/exit MPB environment
'''
def MPB_initialize():
    MPB_CHECKERR( lib.mpb_initialize() )

def MPB_exit():
    lib.mpb_atexit(0)


'''
MPBCones constructor

@param types: list of cone types

@param index_lists: list of lists of indices in each cone
'''
class MPBCones(object):
    def __init__(self, types, lengths, indices):
        self.num = len(types)
        self.types = array(types).astype(int64)
        self.indices = array(indices).astype(int64)
        self.lengths = array(lengths).astype(int64)

        self.type_ptr = ndarray_pointer(self.types)
        self.index_ptr = ndarray_pointer(self.indices)
        self.length_ptr = ndarray_pointer(self.lengths)

'''
wrapper for MathProgBase solver
'''
class MPBSolver(object):
    def __init__(self, packagename, solvername):
        self.ptr = c_void_p()
        MPB_CHECKERR( lib.mpb_new_solver(
            packagename, solvername, byref(self.ptr)) )

    # automatically release solver when no
    # references to MPBSolver object remain
    def __del__(self):
        MPB_CHECKERR( lib.mpb_free_solver(self.ptr) )

'''
wrapper for MathProgBase model
'''
class MPBModel(object):
    '''
    MathProgBase model constructor

    @param packagename:
        string, ---- (TODO: description)

    @param solvername:
        string, ---- (TODO: description)

    @param c:
        problem data, real vector \in R^n

    @param A:
        problem data, real matrix \in R^{m x n}.
        expected in scipy.sparse.coo_matrix format

    @param b:
        problem data, real matrix \in R^m

    @param constrcones:
        description of constraint cones as MPBCones object

    @param varcones:
        description of variable cones as MPBCones object
    '''
    def __init__(self, packagename, solvername,
        c, A, b, constrcones, varcones):

        if not isinstance(A, coo_matrix):
            TypeError("input A must be a scipy.sparse coo_matrix")
        if not isinstance(constrcones, MPBCones):
            TypeError(
                "input constrcones must be an object of type MPBCones")
        if not isinstance(varcones, MPBCones):
            TypeError(
                "input constrcones must be an object of type MPBCones")
        if not sum(constrcones.lengths) == A.shape[0]:
            ValueError("inputs constrcones and A incompatibly sized")
        if not sum(varcones.lengths) == A.shape[1]:
            ValueError("inputs varcones and A incompatibly sized")


        # initialize MathProgBase environment
        MPB_initialize()

        # intialize MathProgBase solver
        self.solver = MPBSolver(packagename, solvername)

        # initialize MathProgBase model
        self.ptr = c_void_p()
        MPB_CHECKERR( lib.mpb_new_model(self.solver.ptr, self.ptr) )

        # load problem data into MathProgBase model
        self.numvar = A.shape[1]
        self.numconstr = A.shape[0]

        row_arr = A.row.astype(int64)
        col_arr = A.col.astype(int64)
        data_arr = A.data.astype(float64)
        b_arr = b.astype(float64)
        c_arr = c.astype(float64)

        MPB_CHECKERR( lib.mpb_loadproblem(self.ptr,
            self.numvar, self.numconstr, ndarray_pointer(c_arr),
            ndarray_pointer(row_arr), ndarray_pointer(col_arr),
            ndarray_pointer(data_arr), A.nnz, ndarray_pointer(b_arr),
            constrcones.num, constrcones.type_ptr,
            constrcones.index_ptr, constrcones.length_ptr,
            varcones.num, varcones.type_ptr,
            varcones.index_ptr, varcones.length_ptr) )

        # create arrays for solution and dual
        self.solution = ndarray(self.numvar, dtype=float64)
        self.dual = ndarray(self.numconstr, dtype=float64)


    '''
    gets int/float properties from model
    '''
    def getproperty(self, property_name):
        if property_name == "numvar":
            call = lib.mpb_numvar
            dtype = int64
        elif property_name == "numconstr":
            call = lib.mpb_numconstr
            dtype = int64
        elif property_name == "objval":
            call = lib.mpb_getobjval
            dtype = float64
        elif property_name == "objbound":
            call = lib.mpb_getobjbound
            dtype = float64
        elif property_name == "pbjgap":
            call = lib.mpb_getobjgap
            dtype = float64
        elif property_name == "solvetime":
            call = lib.mpb_getsolvetime
            dtype = float64
        else:
            print "invalid property key"

        prop = ndarray(1, dtype=dtype)
        MPB_CHECKERR( call(self.ptr, ndarray_pointer(prop)) )
        return prop[0]

    def getsolution(self):
        MPB_CHECKERR( lib.mpb_getsolution(self.ptr,
            ndarray_pointer(self.solution)) )
        return self.solution

    def getdual(self):
        MPB_CHECKERR( lib.mpb_getsolution(self.ptr,
            ndarray_pointer(self.dual)) )
        return self.dual

    def optimize(self):
        MPB_CHECKERR( lib.mpb_optimize(self.ptr) )

    def status(self):
        len_buffer = STATUS_BUFFER_LENGTH
        status_buffer = ndarray(len_buffer, dtype = c_char)
        MPB_CHECKERR( lib.mpb_status(self.ptr,
            ndarray_pointer(status_buffer), len_buffer) )
        return reduce(op_add, status_buffer)


    # automatically free model (and solver) when no
    # references to MPBModel object remain & exit MPB environment
    def __del__(self):
        del self.solver
        MPB_CHECKERR( lib.mpb_free_model(self.ptr) )
        MPB_exit()

