#include <cmpb.h>
#include <julia.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

const double c[] = {-3, -2, -4};
const int64_t I[] = {0,0,1,0,1};
const int64_t J[] = {0,1,1,2,2};
const double V[] = {1.0,1.0,1.0,1.0,1.0};
const double b[] = {3,2};
const int nvar = 3;
const int nconstr = 2;
const int nnz = 5;
const int numconstrcones = 1;
const int64_t constrconetypes[] = { MPBZEROCONE };
const int64_t constrconeindices1[] = { 0, 1 };
const int64_t* constrconeindices[] = { constrconeindices1 };
const int64_t constrconelengths[] = { 2 };
const int numvarcones = 1;
const int64_t varconetypes[] = { MPBNONNEGCONE };
const int64_t varconeindices1[] = { 0, 1, 2 };
const int64_t* varconeindices[] = { varconeindices1 };
const int64_t varconelengths[] = { 3 };

void jl_(jl_value_t*);

int main(int argc, char *argv[])
{
    int i;
    /* required: setup the julia context */
    mpb_initialize();

    jl_eval_string("println(sqrt(2.0))");
    assert(!jl_exception_occurred());
    
    jl_eval_string("import ECOS");
    assert(!jl_exception_occurred());

    jl_value_t *model = NULL;
    JL_GC_PUSH1(&model);
    model = jl_eval_string("MathProgBase.ConicModel(ECOS.ECOSSolver())");
    assert(!jl_exception_occurred());

    int ret = mpb_loadproblem(model, nvar, nconstr, c, I, J, V, nnz, b,
                              numconstrcones, constrconetypes, constrconeindices, constrconelengths,
                              numvarcones, varconetypes, varconeindices, varconelengths);

    ret = mpb_optimize(model);
    
    char status[20];
    ret = mpb_status(model, status, 20);
    printf("STATUS: %s\n",status);

    double objval;
    mpb_getobjval(model, &objval);
    assert(fabs(objval - (-11)) < 1e-3);

    double *sol = (double*) malloc(nvar*sizeof(double));
    mpb_getsolution(model, sol);
    assert(fabs(sol[0]-1.0) < 1e-3);
    assert(fabs(sol[1]-0.0) < 1e-3);
    assert(fabs(sol[2]-2.0) < 1e-3);

    JL_GC_POP();


    /* strongly recommended: notify julia that the
         program is about to terminate. this allows
         julia time to cleanup pending write requests
         and run all finalizers
    */
    jl_atexit_hook(0);
    return 0;
}
