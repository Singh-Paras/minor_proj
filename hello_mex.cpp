// hello_mex.cpp
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[]) {
    mexPrintf("Hello from MEX!\n");
}