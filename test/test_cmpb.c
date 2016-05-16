#include <cmpb.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

void testcont1(void *solver)
{
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
    const int64_t constrconeindices[] = { 0, 1 };
    const int64_t constrconelengths[] = { 2 };
    const int numvarcones = 1;
    const int64_t varconetypes[] = { MPBNONNEGCONE };
    const int64_t varconeindices[] = { 0, 1, 2 };
    const int64_t varconelengths[] = { 3 };

    void *model;
    assert(mpb_supportscone(solver, MPBNONNEGCONE));
    assert(mpb_supportscone(solver, MPBSOC));
    assert(mpb_supportscone(solver, MPBEXPPRIMAL));
    assert(!mpb_supportscone(solver, MPBSDPCONE));

    mpb_new_model(solver, &model);

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

    mpb_free_model(model);
}

/*
    Problem SOCINT1
    min 0x - 2y - 1z
     st  x            == 1
         x >= ||(y,z)||
         (y,z) binary
*/
void testint1(void *solver)
{
    const double c[] = {0, -2, -1};
    const int64_t I[] = {0,};
    const int64_t J[] = {0};
    const double V[] = {1};
    const double b[] = {1};
    const int nvar = 3;
    const int nconstr = 1;
    const int nnz = 1;
    const int numconstrcones = 1;
    const int64_t constrconetypes[] = { MPBZEROCONE };
    const int64_t constrconeindices[] = { 0 };
    const int64_t constrconelengths[] = { 1 };
    const int numvarcones = 1;
    const int64_t varconetypes[] = { MPBSOC };
    const int64_t varconeindices[] = { 0, 1, 2 };
    const int64_t varconelengths[] = { 3 };
    const int64_t vartypes[] = { MPBCONTVAR, MPBBINVAR, MPBBINVAR };

    void *model;

    mpb_new_model(solver, &model);

    int ret = mpb_loadproblem(model, nvar, nconstr, c, I, J, V, nnz, b,
                              numconstrcones, constrconetypes, constrconeindices, constrconelengths,
                              numvarcones, varconetypes, varconeindices, varconelengths);

    ret = mpb_setvartype(model, vartypes, nvar);
    ret = mpb_optimize(model);

    char status[20];
    ret = mpb_status(model, status, 20);
    printf("STATUS: %s\n",status);

    double objval;
    mpb_getobjval(model, &objval);
    assert(fabs(objval - (-2)) < 1e-3);

    double *sol = (double*) malloc(nvar*sizeof(double));
    mpb_getsolution(model, sol);
    assert(fabs(sol[0]-1.0) < 1e-3);
    assert(fabs(sol[1]-1.0) < 1e-3);
    assert(fabs(sol[2]-0.0) < 1e-3);

    mpb_free_model(model);
}

/*
    Problem SOCINT2
    min 0x - 2y - 1z
     st  x            == 1
         x >= ||(y,z)|| # as a constraint
         (y,z) binary
*/
void testint2(void *solver)
{
    const double c[] = {0, -2, -1};
    const int64_t I[] = {0,1,2,3};
    const int64_t J[] = {0,0,1,2};
    const double V[] = {1,-1,-1,-1};
    const double b[] = {1,0,0,0};
    const int nvar = 3;
    const int nconstr = 4;
    const int nnz = 4;
    const int numconstrcones = 2;
    const int64_t constrconetypes[] = { MPBZEROCONE, MPBSOC };
    const int64_t constrconeindices[] = { 0, 1, 2, 3 };
    const int64_t constrconelengths[] = { 1, 3 };
    const int numvarcones = 1;
    const int64_t varconetypes[] = { MPBFREECONE };
    const int64_t varconeindices[] = { 0, 1, 2 };
    const int64_t varconelengths[] = { 3 };
    const int64_t vartypes[] = { MPBCONTVAR, MPBBINVAR, MPBBINVAR };

    void *model;

    mpb_new_model(solver, &model);

    int ret = mpb_loadproblem(model, nvar, nconstr, c, I, J, V, nnz, b,
                              numconstrcones, constrconetypes, constrconeindices, constrconelengths,
                              numvarcones, varconetypes, varconeindices, varconelengths);

    ret = mpb_setvartype(model, vartypes, nvar);
    ret = mpb_optimize(model);

    char status[20];
    ret = mpb_status(model, status, 20);
    printf("STATUS: %s\n",status);

    double objval;
    mpb_getobjval(model, &objval);
    assert(fabs(objval - (-2)) < 1e-3);

    double *sol = (double*) malloc(nvar*sizeof(double));
    mpb_getsolution(model, sol);
    assert(fabs(sol[0]-1.0) < 1e-3);
    assert(fabs(sol[1]-1.0) < 1e-3);
    assert(fabs(sol[2]-0.0) < 1e-3);

    mpb_free_model(model);
}


int main(int argc, char *argv[])
{
    // Set up the julia context
    mpb_initialize();

    void *solver;

    mpb_new_solver("ECOS","ECOSSolver()", &solver);

    testcont1(solver);

    mpb_free_solver(solver);

    //mpb_new_solver("Gurobi","GurobiSolver()", &solver);
    mpb_new_solver("Pajarito, GLPKMathProgInterface, ECOS","PajaritoSolver(verbose=1,mip_solver=GLPKSolverMIP(),cont_solver=ECOSSolver(verbose=0))", &solver);

    // Pajarito doesn't pass int1 test because it has discrete variables in cones
    //testint1(solver);
    testint2(solver);

    mpb_free_solver(solver);

    // julia cleanup 
    mpb_atexit(0);
    return 0;
}
