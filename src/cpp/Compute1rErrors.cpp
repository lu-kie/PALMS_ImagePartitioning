/**
    Compute1rErrors.cpp
    Purpose: Computes the approximation errors for univariate Jet estimation
             with Givens rotations for discrete intervals [1,r] for all r

    @author Lukas Kiefer
    @version 1.0
*/

#include "linewiseAffineMS.h"
#include "Interval.h"

vec Compute1rErrors(mat u_data, mat a_data,mat b_data, const int nr_channels,
                    const mat &C_linear, const mat &S_linear, const mat &C_const, const mat &S_const)
{   
    // Define local variables
    int n = u_data.n_cols;
    vec  eps1R = zeros(n);
    // Create an interval to compute the approximation errors on all intervals [1,r]
    Interval err = Interval(1,1,nr_channels,u_data.col(0),a_data.col(0),b_data.col(0));
    for(int r =1; r < n; r++){
        // Compute the approximation error for interval [1,r] using the error update function of the Interval class
        err.addBottomDataPoint(nr_channels, C_linear,S_linear,C_const, S_const,u_data.col(r),a_data.col(r),b_data.col(r));
        eps1R(r) = err.getEps();
    }
    
    return eps1R;
}
