#include <julia.h>
#include <string.h>

// internal helper functions

jl_function_t* mpb_get_function(const char *name) {
    jl_value_t *mpb = jl_eval_string("MathProgBase");
    jl_function_t *f = jl_get_function((jl_module_t*)mpb, name);
    assert(jl_is_function(f));
    return f;
}

// MUST GC_PUSH the returned value!
jl_value_t* mpb_ptr_to_intvec(int64_t *values, int64_t len) {

    jl_value_t *int_vector_type = jl_apply_array_type(jl_int64_type, 1);

    jl_value_t *arr = (jl_value_t*)jl_alloc_array_1d(int_vector_type, len);
    int64_t *data = (int64_t*)jl_array_data(arr);

    int64_t i;
    for (i = 0; i < len; i++) {
        data[i] = values[i];
    }

    return arr;
}

// MUST GC_PUSH the returned value!
jl_value_t* mpb_ptr_to_floatvec(double *values, int64_t len) {

    jl_value_t *float_vector_type = jl_apply_array_type(jl_float64_type, 1);

    jl_value_t *arr = (jl_value_t*)jl_alloc_array_1d(float_vector_type, len);
    double *data = (double*)jl_array_data(arr);

    int64_t i;
    for (i = 0; i < len; i++) {
        data[i] = values[i];
    }

    return arr;
}

void mpb_fill_floatvec(jl_value_t *arr, double *output, int64_t len) {

    jl_value_t *float_vector_type = jl_apply_array_type(jl_float64_type, 1);
    assert(jl_typeis(sol, float_vector_type));
    double *data = (double*)jl_array_data(sol);
    assert(jl_array_len(arr) >= len);

    int i;
    double *data = (double*)jl_array_data(sol);
    for (i = 0; i < numvar; i++) {
        output[i] = data[i];
    }
}

// MUST GC_PUSH the returned value!
jl_value_t* mpb_triplet_to_sparse(int64_t *I, int64_t *J, double *V, int64_t nnz, int64_t m, int64_t n) {

    // in the same order as the arguments to sparse()
    jl_value_t **objects;
    JL_GC_PUSHARGS(objects, 5);

    objects[1] = mpb_ptr_to_intvec(I,nnz);
    objects[2] = mpb_ptr_to_intvec(J,nnz);
    objects[3] = mpb_ptr_to_floatvec(V,nnz);
    objects[4] = jl_box_int64(m);
    objects[5] = jl_box_int64(n);

    jl_function_t *sparse_f = jl_get_function(jl_base_module, "sparse");
    jl_value_t *spmat = jl_call(sparse_f, objects, 5);
    assert(!jl_exception_occurred());

    JL_GC_POP();

    return spmat;
}

#define MPBFREECONE 0
#define MPBZEROCONE 1
#define MPBNONNEGCONE 2
#define MPBNONPOSCONE 3
#define MPBSOC 4
#define MPBSOCROTATED 5
#define MPBSDPCONE 6
#define MPBEXPPRIMAL 7
#define MPBEXPDUAL 8

jl_value_t* int_to_symbol(int64_t val) {
    switch (val) {
        case MPBFREECONE:
            return (jl_value_t*)jl_symbol("Free");
        case MPBZEROCONE:
            return (jl_value_t*)jl_symbol("Zero");
        case MPBNONNEGCONE:
            return (jl_value_t*)jl_symbol("NonNeg");
        case MPBNONPOSCONE:
            return (jl_value_t*)jl_symbol("NonPos");
        case MPBSOC:
            return (jl_value_t*)jl_symbol("SOC");
        case MPBSOCROTATED:
            return (jl_value_t*)jl_symbol("SOCRotated");
        case MPBSDPCONE:
            return (jl_value_t*)jl_symbol("SDP");
        case MPBEXPPRIMAL:
            return (jl_value_t*)jl_symbol("ExpPrimal");
        case MPBEXPDUAL:
            return (jl_value_t*)jl_symbol("ExpDual");
        default:
            assert(0 && "Unrecognized cone value");
            return (jl_value_t*)jl_symbol("Error");
    }
}

// MUST GC_PUSH the returned value
// Generates a vector with a list of cones, to be used as input for loadproblem!
jl_value_t* mpb_conevector(int64_t numcones, int64_t *conetypes, int64_t **coneindices, int64_t *conelengths) {

    jl_value_t *conevector = NULL, *indexvector = NULL, *conesymbol = NULL, *conetuple = NULL;
    JL_GC_PUSH4(&conevector, &indexvector, &conesymbol, &conetuple);
    conevector = jl_eval_string("Array(Tuple{Symbol,Any},0)");
    jl_function_t *push_f = jl_get_function(jl_base_module, "push!");
    int (i = 0; i < numcones; i++) {
        conesymbol = int_to_symbol(conetypes[i]);
        indexvector = mpb_ptr_to_intvec(coneindices[i], conelengths[i]);
        conetuple = jl_new_struct(jl_eval_string("Tuple{Symbol,Vector{Int64}}"), conesymbol, indexvector);
        assert(!jl_exception_occurred());
        jl_call2(push_f, conevector, conetuple); 
        assert(!jl_exception_occurred());
    }

    JL_GC_POP();

    return conevector;
}


// public functions

int mpb_initialize() {

    if (!jl_is_initialized()) {
        jl_init(JULIA_INIT_DIR);
        jl_eval_string("import MathProgBase");
        assert(!jl_exception_occurred());
        jl_value_t *mpb = jl_eval_string("MathProgBase");
        assert(jl_is_module(mpb));

        // TODO: set up atexit callback
    }

    return 0; // success
}

#define MPBGetIntProperty(prop) \
int \
mpb_ ## prop (void *model, int64_t *output) \
{ \
    jl_value_t *val = jl_call1(mpb_get_function(prop), (jl_value_t*)model); \
    assert(!jl_exception_occurred()); \
    assert(jl_is_int64(val)); \
    *output = jl_unbox_int64(val); \
    return 0; \
}

#define MPBGetFloatProperty(prop) \
int \
mpb_ ## prop (void *model, double *output) \
{ \
    jl_value_t *val = jl_call1(mpb_get_function(prop), (jl_value_t*)model); \
    assert(!jl_exception_occurred()); \
    assert(jl_is_float64(val)); \
    *output = jl_unbox_float64(val); \
    return 0; \
}

MPBGetIntProperty(numvar);
MPBGetIntProperty(numconstr);

MPBGetFloatProperty(getobjval);
MPBGetFloatProperty(getobjbound);
MPBGetFloatProperty(getobjgap);
MPBGetFloatProperty(getsolvetime);

int64_t numcones, int64_t *conetypes, int64_t **coneindices, int64_t *conelengths

int mpb_loadproblem(void *model, // model pointer
                    int64_t numvar, // number of variables
                    int64_t numconstr, // number of rows in A
                    double *c, // objective coefficient vector
                    int64_t *A_I, // row indices of A in triplet format
                    int64_t *A_J, // column indices of A in triplet format
                    double *A_V, // nonzero values of A in triplet format
                    int64_t A_nnz, // lengths of A_I, A_J, and A_V
                    double *b, // right-hand side vector
                    int64_t numconstrcones, // number of constraint cones
                    int64_t *constrconetypes, // types of each constraint cone
                    int64_t **constrconeindices, // vector of indices for each constraint cone
                    int64_t *constrconelengths, // number of indices in each constraint cone
                    int64_t numvarcones, // number of variable cones
                    int64_t *varconetypes, // types of each variable cone
                    int64_t **varconeindices, // vector of indices for each variable cone
                    int64_t *varconelengths // number of indices in each variable cone
                    ) {

    jl_value_t *cvec = NULL, *Amat = NULL, *bvec = NULL, *constr_cones = NULL,
               *var_cones = NULL;
    JL_GC_PUSH5(&cvec, &Amat, &bvec, &constr_cones, &var_cones);

    cvec         = mpb_ptr_to_floatvec(c, numvar);
    Amat         = mpb_triplet_to_sparse(A_I, A_J, A_V, A_nnz, numconstr, numvar);
    bvec         = mpb_ptr_to_floatvec(b, numconstr);
    constr_cones = mpb_conevector(numconstrcones, constrconetypes,
                               constrconeindices, constrconelengths);
    var_cones    = mpb_conevector(numvarcones, varconetypes,
                               varconeindices, varconelengths);

    jl_function_t *loadproblem_f = mpb_get_function("loadproblem!");
    jl_value_t *loadargs[] = { (jl_value_t*)model, cvec, Amat, bvec, constr_cones, var_cones };

    jl_call(loadproblem_f, loadargs, 6);
    assert(!jl_exception_occurred());

    JL_GC_POP();

    return 0;
}

// length of output should be at least numvar
int mpb_getsolution(void *model, double *output) {

    jl_value_t *sol = jl_call1(mpb_get_function("getsolution"), (jl_value_t*)model);
    assert(!jl_exception_occurred())
    JL_GC_PUSH(&sol);

    int64_t numvar;
    int ret = mpb_numvar(model, &numvar);
    assert(ret == 0);

    mpb_fill_floatvec(sol, output, numvar);

    JL_GC_POP();

    return 0;
}

// length of output should be at least numconstr
int mpb_getdual(void *model, double *output) {

    jl_value_t *sol = jl_call1(mpb_get_function("getdual"), (jl_value_t*)model);
    assert(!jl_exception_occurred())
    JL_GC_PUSH(&sol);

    int64_t numconstr;
    int ret = mpb_numconstr(model, &numconstr);
    assert(ret == 0);

    mpb_fill_floatvec(sol, output, numconstr);

    JL_GC_POP();

    return 0;
}

int mpb_optimize(void *model) {

    jl_call1(mpb_get_function("optimize!"), (jl_value_t*)model);
    assert(!jl_exception_occurred());

    return 0;
}

// len is length of output vector
int mpb_status(void *model, char *output, int64_t len) {

    jl_value_t *status = jl_call1(mpb_get_function("status"), (jl_value_t*)model);
    assert(!jl_exception_occurred());
    assert(jl_is_symbol(status));
    const char *name = ((jl_sym_t*)status)->name; // TODO: replace with jl_symbol_name

    assert(strlen(name) <= len);

    strncpy(output, name, len);

    return 0;
}
