#include <julia.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

const double c[] = {-3, -2, -4};
const int64_t I[] = {1,1,2,1,2};
const int64_t J[] = {1,2,2,3,3};
const double V[] = {1.0,1.0,1.0,1.0,1.0};
const double b[] = {3,2};
const int nvar = 3;
const int nconstr = 2;
const int nnz = 5;

void jl_(jl_value_t*);

int main(int argc, char *argv[])
{
    /* required: setup the julia context */
    jl_init(JULIA_INIT_DIR);

    jl_eval_string("println(sqrt(2.0))");
    jl_eval_string("import MathProgBase");
    assert(!jl_exception_occurred());
    jl_value_t *mpb = jl_eval_string("MathProgBase");
    assert(jl_is_module(mpb));

    jl_function_t *ConicModel = jl_get_function((jl_module_t*)mpb, "ConicModel");
    assert(jl_is_function(ConicModel));
    
    // constr_cones and var_cones are Vector{Tuple{Symbol,Any}}
    jl_value_t** objects;
    JL_GC_PUSHARGS(objects, 20);
    // var_cones, constr_cones, zerocone, nonnegcone
   
    objects[0] = jl_eval_string("Array(Tuple{Symbol,Any},0)");
    objects[1] = jl_eval_string("Array(Tuple{Symbol,Any},0)");
    //JL_GC_PUSH1(&constr_cones);

    objects[2] = (jl_value_t*)jl_symbol("Zero");
    objects[3] = (jl_value_t*)jl_symbol("NonNeg");
    //JL_GC_PUSH1(&zerocone)
    //JL_GC_PUSH1(&nonnegcone);

    objects[4] = jl_new_struct(jl_eval_string("Tuple{Symbol,Symbol}"),objects[2],objects[3]);
    //JL_GC_PUSH1(&tup);
    //jl_show(JL_STDOUT,tup);
    jl_(objects[4]);

    // arrays
    jl_value_t* float_vector_type = jl_apply_array_type(jl_float64_type, 1);
    jl_value_t* int_vector_type = jl_apply_array_type(jl_int64_type, 1);
    // c -- objective vector
    objects[5] = (jl_value_t*)jl_alloc_array_1d(float_vector_type, nvar);
    double *data = (double*)jl_array_data(objects[5]);
    for (int i = 0; i < nvar; i++)
        data[i] = c[i];
    
    // b -- constraint rhs
    objects[6] = (jl_value_t*)jl_alloc_array_1d(float_vector_type, nconstr);
    data = (double*)jl_array_data(objects[6]);
    for (int i = 0; i < nconstr; i++)
        data[i] = b[i];

    // I,J,V -- sparse matrix in triplet form
    objects[7] = (jl_value_t*)jl_alloc_array_1d(int_vector_type, nnz); 
    objects[8] = (jl_value_t*)jl_alloc_array_1d(int_vector_type, nnz); 
    objects[9] = (jl_value_t*)jl_alloc_array_1d(float_vector_type, nnz);
    int64_t *dataI = (int64_t*)jl_array_data(objects[7]);
    int64_t *dataJ = (int64_t*)jl_array_data(objects[8]);
    double *dataV = (double*)jl_array_data(objects[9]);

    for (int i = 0; i < nnz; i++) {
        dataI[i] = I[i];
        dataJ[i] = J[i];
        dataV[i] = V[i];
    }
    
    jl_function_t *sparse_f = jl_get_function(jl_base_module, "sparse");
    
    objects[10] = jl_box_int64(nconstr);
    objects[11] = jl_box_int64(nvar);
    objects[12] = jl_call(sparse_f, objects+7,5);

    jl_(objects[7]);
    jl_(objects[8]);
    jl_(objects[9]);
    jl_(objects[12]);
    

    // cones. take the easy way for now
    objects[13] = jl_eval_string("[(:Zero,1:2)]");
    objects[14] = jl_eval_string("[(:NonNeg,1:3)]");

    jl_eval_string("import ECOS");
    assert(!jl_exception_occurred());

    objects[15] = jl_eval_string("MathProgBase.ConicModel(ECOS.ECOSSolver())");
    assert(!jl_exception_occurred());

    jl_function_t *loadproblem_f = jl_get_function((jl_module_t*)mpb, "loadproblem!");
    jl_value_t *loadargs[] = { objects[15], objects[5], objects[12], objects[6], objects[13], objects[14] };
    jl_call(loadproblem_f, loadargs,6);
    jl_function_t *optimize_f = jl_get_function((jl_module_t*)mpb, "optimize!");
    jl_call1(optimize_f, objects[15]);
    
    jl_function_t *status_f = jl_get_function((jl_module_t*)mpb, "status");
    jl_value_t* status = jl_call1(status_f, objects[15]);
    assert(jl_is_symbol(status));
    jl_(status);


    jl_function_t *getobjval_f = jl_get_function((jl_module_t*)mpb, "getobjval");
    jl_value_t *objval = jl_call1(getobjval_f, objects[15]);
    assert(jl_is_float64(objval));
    double objval_d = jl_unbox_float64(objval);
    assert(fabs(objval_d - (-11)) < 1e-3);

    jl_function_t *getsolution_f = jl_get_function((jl_module_t*)mpb, "getsolution");
    jl_value_t *sol = jl_call1(getsolution_f, objects[15]);
    assert(jl_typeis(sol, float_vector_type));
    double *solval = (double*)jl_array_data(sol);
    assert(fabs(solval[0]-1.0) < 1e-3);
    assert(fabs(solval[1]-0.0) < 1e-3);
    assert(fabs(solval[2]-2.0) < 1e-3);

    JL_GC_POP();


    //jl_eval_string("import SCS");

    /* strongly recommended: notify julia that the
         program is about to terminate. this allows
         julia time to cleanup pending write requests
         and run all finalizers
    */
    jl_atexit_hook(0);
    return 0;
}
