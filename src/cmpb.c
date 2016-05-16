#include <julia.h>
#include <string.h>
#include <stdio.h>
#include <cmpb.h>


void jl_(jl_value_t*); // for debugging

// internal helper functions

jl_function_t* mpb_get_function(const char *name) {
    jl_value_t *mpb = jl_eval_string("MathProgBase");
    jl_function_t *f = jl_get_function((jl_module_t*)mpb, name);
    assert(jl_is_function(f));
    return f;
}

// MUST GC_PUSH the returned value!
// Adds 1 to all terms!
jl_value_t* mpb_ptr_to_intvec(const int64_t *values, int64_t len) {

    jl_value_t *int_vector_type = jl_apply_array_type(jl_int64_type, 1);

    jl_value_t *arr = (jl_value_t*)jl_alloc_array_1d(int_vector_type, len);
    int64_t *data = (int64_t*)jl_array_data(arr);

    int64_t i;
    for (i = 0; i < len; i++) {
        data[i] = values[i]+1;
    }

    return arr;
}

// MUST GC_PUSH the returned value!
jl_value_t* mpb_ptr_to_floatvec(const double *values, int64_t len) {

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
    assert(jl_typeis(arr, float_vector_type));
    assert(jl_array_len(arr) >= len);

    int i;
    double *data = (double*)jl_array_data(arr);
    for (i = 0; i < len; i++) {
        output[i] = data[i];
    }
}

// MUST GC_PUSH the returned value!
jl_value_t* mpb_triplet_to_sparse(const int64_t *I, const int64_t *J, const double *V, int64_t nnz, int64_t m, int64_t n) {

    // in the same order as the arguments to sparse()
    jl_value_t **objects;
    JL_GC_PUSHARGS(objects, 5);

    objects[0] = mpb_ptr_to_intvec(I,nnz);
    objects[1] = mpb_ptr_to_intvec(J,nnz);
    objects[2] = mpb_ptr_to_floatvec(V,nnz);
    objects[3] = jl_box_int64(m);
    objects[4] = jl_box_int64(n);

    jl_function_t *sparse_f = jl_get_function(jl_base_module, "sparse");
    jl_value_t *spmat = jl_call(sparse_f, objects, 5);
    assert(!jl_exception_occurred());

    JL_GC_POP();

    return spmat;
}


jl_value_t* int_to_conesymbol(int64_t val) {
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

jl_value_t* int_to_categorysymbol(int64_t val) {
    switch (val) {
        case MPBCONTVAR:
            return (jl_value_t*)jl_symbol("Cont");
        case MPBINTVAR:
            return (jl_value_t*)jl_symbol("Int");
        case MPBBINVAR:
            return (jl_value_t*)jl_symbol("Bin");
        case MPBSEMICONTVAR:
            return (jl_value_t*)jl_symbol("SemiCont");
        case MPBSEMIINTVAR:
            return (jl_value_t*)jl_symbol("SemiInt");
        default:
            assert(0 && "Unrecognized variable category value");
            return (jl_value_t*)jl_symbol("Error");
    }
}

// MUST GC_PUSH the returned value
// Generates a vector with a list of cones, to be used as input for loadproblem!
jl_value_t* mpb_conevector(int64_t numcones, const int64_t *conetypes, const int64_t *coneindices, const int64_t *conelengths) {

    jl_value_t *conevector = NULL, *indexvector = NULL, *conesymbol = NULL, *conetuple = NULL;
    JL_GC_PUSH4(&conevector, &indexvector, &conesymbol, &conetuple);
    conevector = jl_eval_string("Array(Tuple{Symbol,Any},0)");
    jl_function_t *push_f = jl_get_function(jl_base_module, "push!");
    int i;
    int offset = 0;
    for (i = 0; i < numcones; i++) {
        conesymbol = int_to_conesymbol(conetypes[i]);
        indexvector = mpb_ptr_to_intvec(coneindices + offset, conelengths[i]);
        conetuple = jl_new_struct(jl_eval_string("Tuple{Symbol,Vector{Int64}}"), conesymbol, indexvector);
        assert(!jl_exception_occurred());
        jl_call2(push_f, conevector, conetuple);
        assert(!jl_exception_occurred());
        offset = offset + conelengths[i];
    }

    JL_GC_POP();

    return conevector;
}

// MUST GC_PUSH the returned value
jl_value_t* mpb_objectid(jl_value_t *object) {

    jl_function_t *objectid_f = jl_get_function(jl_base_module, "object_id");
    jl_value_t *id = jl_call1(objectid_f, object);
    return id;
}

void mpb_register_object(jl_value_t *object) {

    jl_value_t *object_dict = jl_eval_string("__mpb_object_dict");
    assert(jl_typeis(object_dict, jl_eval_string("Dict{Any,Any}")));
    jl_function_t *setindex_f = jl_get_function(jl_base_module, "setindex!");
    jl_call3(setindex_f, object_dict, mpb_objectid(object), object);
}

void mpb_release_object(jl_value_t *object) {

    jl_value_t *object_dict = jl_eval_string("__mpb_object_dict");
    assert(jl_typeis(object_dict, jl_eval_string("Dict{Any,Any}")));
    jl_function_t *pop_f = jl_get_function(jl_base_module, "pop!");
    jl_call2(pop_f, object_dict,mpb_objectid(object));
}

// public functions

int mpb_initialize() {

    if (!jl_is_initialized()) {
        jl_init(JULIA_INIT_DIR);
        jl_eval_string("import MathProgBase");
        assert(!jl_exception_occurred());
        jl_value_t *mpb = jl_eval_string("MathProgBase");
        assert(jl_is_module(mpb));
        // global dictionary of objects we've given out to users
        // this prevents the GC from freeing them
        jl_eval_string("__mpb_object_dict = Dict{Any,Any}()");
    }

    return 0; // success
}

int mpb_new_solver(const char *packagename, const char *solvername, void **output) {
    char tmp[200];
    snprintf(tmp, 200, "using %s", packagename);
    jl_eval_string(tmp);
    assert(!jl_exception_occurred());

    jl_value_t *solver = jl_eval_string(solvername);
    JL_GC_PUSH1(&solver);
    mpb_register_object(solver);
    JL_GC_POP();

    *output = solver;

    return 0;
}


int mpb_supportscone(void *solver, int64_t conetype) {

    jl_value_t *cones = NULL, *sym = NULL;
    JL_GC_PUSH2(&solver, &sym);
    cones = jl_call1(mpb_get_function("supportedcones"), solver);
    assert(!jl_exception_occurred());
    sym = int_to_conesymbol(conetype);

    jl_function_t *in_f = jl_get_function(jl_base_module, "in");
    jl_value_t * isin = jl_call2(in_f, sym, cones);
    assert(!jl_exception_occurred());
    assert(jl_is_bool(isin));
    int8_t isin_b = jl_unbox_bool(isin);

    JL_GC_POP();

    return isin_b == 1;
}

// not very important to free the solver, it usually doesn't
// take up any memory
int mpb_free_solver(void *solver) {

    mpb_release_object(solver);
    return 0;
}

int mpb_new_model(void *solver, void **output) {

    jl_value_t * model = jl_call1(mpb_get_function("ConicModel"), solver);
    JL_GC_PUSH1(&model);
    mpb_register_object(model);
    JL_GC_POP();

    *output = model;

    return 0;
}

int mpb_free_model(void *model) {

    mpb_release_object(model);
    return 0;
}

// julia cleanup:
// ensures a clean shutdown, solvers properly released
// use exitcode = 0 for normal case
void mpb_atexit(int exitcode) {
    jl_atexit_hook(exitcode);
}

#define MPBGetIntProperty(prop) \
int \
mpb_ ## prop (void *model, int64_t *output) \
{ \
    jl_value_t *val = jl_call1(mpb_get_function(#prop), (jl_value_t*)model); \
    assert(!jl_exception_occurred()); \
    assert(jl_is_int64(val)); \
    *output = jl_unbox_int64(val); \
    return 0; \
}

#define MPBGetFloatProperty(prop) \
int \
mpb_ ## prop (void *model, double *output) \
{ \
    jl_value_t *val = jl_call1(mpb_get_function(#prop), (jl_value_t*)model); \
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


int mpb_loadproblem(void *model, // model pointer
                    int64_t numvar, // number of variables
                    int64_t numconstr, // number of rows in A
                    const double *c, // objective coefficient vector
                    const int64_t *A_I, // row indices of A in triplet format
                    const int64_t *A_J, // column indices of A in triplet format
                    const double *A_V, // nonzero values of A in triplet format
                    int64_t A_nnz, // lengths of A_I, A_J, and A_V
                    const double *b, // right-hand side vector
                    int64_t numconstrcones, // number of constraint cones
                    const int64_t *constrconetypes, // types of each constraint cone
                    const int64_t *constrconeindices, // vector of indices for each constraint cone
                    const int64_t *constrconelengths, // number of indices in each constraint cone
                    int64_t numvarcones, // number of variable cones
                    const int64_t *varconetypes, // types of each variable cone
                    const int64_t *varconeindices, // vector of indices for each variable cone
                    const int64_t *varconelengths // number of indices in each variable cone
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
    assert(!jl_exception_occurred());
    JL_GC_PUSH1(&sol);

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
    assert(!jl_exception_occurred());
    JL_GC_PUSH1(&sol);

    int64_t numconstr;
    int ret = mpb_numconstr(model, &numconstr);
    assert(ret == 0);

    mpb_fill_floatvec(sol, output, numconstr);

    JL_GC_POP();

    return 0;
}

// assuming vartypes is of length len
int mpb_setvartype(void *model, const int64_t *vartypes, int64_t len) {
    int64_t numvar;
    int ret = mpb_numvar(model, &numvar);
    assert(ret == 0);
    assert(len <= numvar);

    jl_value_t *catvector = NULL, *catsymbol;
    JL_GC_PUSH1(&catvector);
    jl_function_t *push_f = jl_get_function(jl_base_module, "push!");
    catvector = jl_eval_string("Array(Symbol,0)");

    int i;
    for (i = 0; i < len; i++) {
        catsymbol = int_to_categorysymbol(vartypes[i]);
        jl_call2(push_f, catvector, catsymbol);
        assert(!jl_exception_occurred());
    }

    jl_call2(mpb_get_function("setvartype!"), (jl_value_t*)model, catvector);
    assert(!jl_exception_occurred());

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

// STUB: replace this with Julia call to check package
//  (and load it to ensure precompiled?)
int mpb_checkpackage(const char *packagename){
    // int ret = JULIACALL( packagename in keys(Pkg.installed )
    // if ret
    //     JULIACALL( import packagename )
    // return ret
    return 0;
}


