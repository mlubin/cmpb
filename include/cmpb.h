#ifndef ___CMPB_H__
#define ___CMPB_H__


#include <stdint.h>

#define MPBFREECONE 0
#define MPBZEROCONE 1
#define MPBNONNEGCONE 2
#define MPBNONPOSCONE 3
#define MPBSOC 4
#define MPBSOCROTATED 5
#define MPBSDPCONE 6
#define MPBEXPPRIMAL 7
#define MPBEXPDUAL 8

int mpb_initialize(void);
void mpb_atexit(int exitcode);

int mpb_numvar(void *model,int64_t *output);
int mpb_numconstr(void *model,int64_t *output);

int mpb_getobjval(void *model,double *output);
int mpb_getobjbound(void *model,double *output);
int mpb_getobjgap(void *model,double *output);
int mpb_getsolvetime(void *model,double *output);

int mpb_new_solver(const char *packagename, const char *solvername, void **output);
int mpb_free_solver(void *solver);
int mpb_new_model(void *solver, void **output);
int mpb_free_model(void *model);
void mpb_atexit(int exitcode);


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
                    const int64_t **constrconeindices, // vector of indices for each constraint cone
                    const int64_t *constrconelengths, // number of indices in each constraint cone
                    int64_t numvarcones, // number of variable cones
                    const int64_t *varconetypes, // types of each variable cone
                    const int64_t **varconeindices, // vector of indices for each variable cone
                    const int64_t *varconelengths // number of indices in each variable cone
                    );

int mpb_getsolution(void *model, double *output);
int mpb_getdual(void *model, double *output);
int mpb_optimize(void *model);
int mpb_status(void *model, char *output, int64_t len);


#endif
