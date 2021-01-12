/**
    LinewiseSolver_mexWrapper.cpp
    Purpose: Call from Matlab to apply linewise affine linear partitioning

    @author Lukas Kiefer
    @version 0.9
*/

#include "linewiseAffineMS.h"
#include "mex.h"

void mexFunction(int nlhs,  mxArray *plhs[], int nrhs,
        const mxArray *prhs[])
{
    // Shortcuts
	#define U_OUT       plhs[0]
	#define A_OUT       plhs[1]
	#define B_OUT       plhs[2]
    
	#define F_IN        prhs[0]
	#define A_DATA_IN	prhs[1]
	#define B_DATA_IN	prhs[2]
	#define DIR_IN      prhs[3]
	#define GAMMA_IN	prhs[4]
	#define ETA_IN		prhs[5]
	#define C_LINEAR_IN prhs[6]
	#define S_LINEAR_IN	prhs[7]
	#define C_CONST_IN	prhs[8]
	#define S_CONST_IN	prhs[9]
    #define NR_THREADS_IN	prhs[10]
    
    // Determine data dimensions
    const mwSize  nr_dims  = mxGetNumberOfDimensions(F_IN);
    const mwSize *dataDims = mxGetDimensions(F_IN);
    
    const long unsigned int m           = dataDims[0];
    const long unsigned int n           = dataDims[1];
    // Check for multi- or single channel image
    int n_ch;
    if (nr_dims < 3)
        n_ch = 1;
    else
        n_ch = dataDims[2];
    
    const long unsigned int nr_channels = n_ch;
    // Model parameters
    const double gamma_s	= mxGetScalar(GAMMA_IN);
    const double eta_s      = mxGetScalar(ETA_IN);
    // Number of threads for openMP (Multicore)
    const int nr_threads = mxGetScalar(NR_THREADS_IN);
    // Get Pointers
    double* data_raw   = mxGetPr(F_IN);
    double* dir_raw     = mxGetPr(DIR_IN);
    double* C_mixed_raw = mxGetPr(C_LINEAR_IN);
    double* S_mixed_raw = mxGetPr(S_LINEAR_IN);
    double* C_const_raw = mxGetPr(C_CONST_IN);
    double* S_const_raw = mxGetPr(S_CONST_IN);
    double* a_data_raw = mxGetPr(A_DATA_IN);
    double* b_data_raw = mxGetPr(B_DATA_IN);
            
    
    
    if(mxGetM(DIR_IN) != 2)
        mexErrMsgTxt("Each direction vector must have length 2");
    
    // Create output
    const mwSize output_dims[3] = {m,n,nr_channels};
    
    U_OUT = mxCreateNumericArray(nr_dims,output_dims,mxDOUBLE_CLASS,mxREAL);
    A_OUT = mxCreateNumericArray(nr_dims,output_dims,mxDOUBLE_CLASS,mxREAL);
    B_OUT = mxCreateNumericArray(nr_dims,output_dims,mxDOUBLE_CLASS,mxREAL);
    
    // Create mem pointers
    double* u_raw = mxGetPr(U_OUT);
    double* a_out_raw = mxGetPr(A_OUT);
    double* b_out_raw = mxGetPr(B_OUT);
   
    // Call frame function
    ArmadilloConverter(data_raw,a_data_raw,b_data_raw,
                                  C_mixed_raw,S_mixed_raw,C_const_raw,S_const_raw,
                                  u_raw,a_out_raw,b_out_raw,
                                  dir_raw,gamma_s,eta_s,m,n,nr_channels,nr_threads);
    
    return;
}
 