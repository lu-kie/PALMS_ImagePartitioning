/**
    FindBest1DPartition.cpp
    Purpose: Computes the best (piecewise affine-linear) partition for univariate Jet data
             (u_data,a_data,b_data) and jump penalty gamma by dynamic programming

    @author Lukas Kiefer
    @version 1.0
*/
#include "Interval.h"
#include "linewiseAffineMS.h"

void FindBest1DPartition(mat u_data,mat a_data,mat b_data,const int n, const int nr_channels, double &gamma,
                         double eta, const vec &eps_1r, mat &C_linear, mat &S_linear,
                         mat &C_const, mat &S_const, ivec &L)
{
    // Allocate vector with optimal functional values for each r=1,...,n
    vec B = zeros(n,1);
    L(0) = 0; 
    // List containing possible segments, i.e., discrete intervals
    std::list<Interval> segments;
    segments.push_front(Interval(2,2,nr_channels,eta*u_data.col(1),a_data.col(1),b_data.col(1)));
    // Local aux variables
    double b;
    vec udata_new(nr_channels);
    vec adata_new(nr_channels);
    vec bdata_new(nr_channels);
    
    for(int r=2; r<=n; r++) {
        // Init with approximation error of single-segment partition, i.e. l = 1:
        B(r-1) = eps_1r(r-1);
        L(r-1) = 0;

        // Loop (backwards in l) through candidates for (best) last changepoint 
        for(list<Interval>::iterator it = segments.begin(); it != segments.end(); ++it) {
            Interval &curr_interval = *it;
            while (curr_interval.getR() < r){
                // Get data of new index r
                udata_new = eta*u_data.col(curr_interval.getR());
                adata_new = a_data.col(curr_interval.getR());
                bdata_new = b_data.col(curr_interval.getR());
                // Extend current interval by new data and update its approximation error with Givens rotations
                curr_interval.addBottomDataPoint(nr_channels, C_linear,S_linear,C_const, S_const,udata_new,adata_new,bdata_new);
            }
            // Check if current interval has better energy
            b = B(curr_interval.getL() - 2) + gamma + curr_interval.getEps();
            if (b <= B(r-1)) {
                B(r-1) = b;
                L(r-1) = curr_interval.getL()-1;
            }
            // Pruning-strategy (omit unnecessary computations of approximation errors)
            if (curr_interval.getEps()+gamma > B(r-1)){
                break;
            }
            
        }
        // Add interval with left r bound to the list of candidates
        if (r<=n-1) {
            // add l=r+1
            segments.push_front(Interval(r+1,r+1,nr_channels,eta*u_data.col(r),a_data.col(r),b_data.col(r)));
        }

    }
}
