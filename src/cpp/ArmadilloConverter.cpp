/**
    ArmadilloConverter.cpp
    Purpose: Converter from pointers to armadillo objects 
             called by the mexFunction LinewiseSolver_mexWrapper from MATLAB

    @author Lukas Kiefer
    @version 1.0
*/

#include "linewiseAffineMS.h"

void ArmadilloConverter(double* u_data_raw, double* a_data_raw,double* b_data_raw,
                        double* C_linear_raw,double* S_linear_raw,
                        double* C_const_raw, double* S_const_raw,
                        double* u_out_raw, double* a_out_raw, double* b_out_raw,
                        double* dir_raw, const double gamma_s,const double eta_s,
                        const int m, const int n, const int nr_channels,const int nr_threads)
{
    omp_set_num_threads(nr_threads);
    // Max possible size of an image scanline
    const int max_stripe_length = max(m,n);
    // Create Armadillo objects from raw pointers on MATLAB objects
    // The used constructors make sure that the memory corresponding to the pointers is used
    // Image data
    cube u_data = cube(u_data_raw,m,n,nr_channels,false,true);
    // Slope data
    cube a_data = cube(a_data_raw,m,n,nr_channels,false,true);
    cube b_data = cube(b_data_raw,m,n,nr_channels,false,true);
    // Givens rotation coefficients
    mat C_linear = mat(C_linear_raw,2*max_stripe_length,2,false,true);
    mat S_linear = mat(S_linear_raw,2*max_stripe_length,2,false,true);
    mat C_const = mat(C_const_raw,max_stripe_length,1,false,true);
    mat S_const = mat(S_const_raw,max_stripe_length,1,false,true);
    // Direction of lines
    vec dir = vec(dir_raw,2,false,true);
    // Output
    // Pathwise regularized image
    cube u_out = cube(u_out_raw,m,n,nr_channels,false,true);
    // Corresponding regularized slopes
    cube a_out = cube(a_out_raw,m,n,nr_channels,false,true);
    cube b_out = cube(b_out_raw,m,n,nr_channels,false,true);

    // Call linewise function
    LinewisePartitioning(u_out,a_out,b_out,m,n,nr_channels,u_data,a_data,b_data,dir,gamma_s,eta_s,C_linear,S_linear,C_const,S_const);
}