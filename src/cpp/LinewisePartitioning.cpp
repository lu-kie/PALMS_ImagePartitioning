/**
    LinewisePartitioning.cpp
    Purpose: Applies univariate affine Jet-partitioning to all lines of the 
             input image for the input direction

    @author Lukas Kiefer
    @version 1.0
*/

#include "linewiseAffineMS.h"

void LinewisePartitioning(cube &u_out,cube &a_out,cube &b_out,
                          const int m,const int n,const int nr_channels,
                          cube &u_data,cube &a_data,cube &b_data,vec &dir, double gamma_s,
                          double eta_s,mat &C_linear, mat &S_linear, mat &C_const, mat &S_const)
{  
    
    // Create lists with stripes along direction dir
    vector<Stripe> L_udata, L_adata,L_bdata;
    Extract1Dstripes(u_data,dir,L_udata,m,n,nr_channels);
    Extract1Dstripes(a_data,dir,L_adata,m,n,nr_channels);
    Extract1Dstripes(b_data,dir,L_bdata,m,n,nr_channels);
   
    // Solve univariate partitioning problems along the lines induced by dir
    #pragma omp parallel for schedule(dynamic)
    for(unsigned int iter = 0; iter < L_udata.size(); ++iter) {
        // All lists share the same iterator iter
        Stripe &stripe_udata     = L_udata[iter];
        Stripe &stripe_adata     = L_adata[iter];
        Stripe &stripe_bdata     = L_bdata[iter];
        // Length of current 1D-problem
        int stripe_length = stripe_udata.giveLength();
        // Linear indices of current stripe within 2D domain
        uvec linear_indices = stripe_udata.getIndices();
        // Catch stripes of length 1
        if(stripe_length < 2) {
            u_out.elem(linear_indices) = stripe_udata.getData().t();
            a_out.elem(linear_indices) = stripe_adata.getData().t();
            b_out.elem(linear_indices) = stripe_bdata.getData().t();
            continue;
        }

        // The 1D partition is encoded by the vector L
        ivec L(stripe_length);

        // [1,r]-errors
        vec Eps1R = Compute1rErrors(eta_s*stripe_udata.getData(),stripe_adata.getData(),stripe_bdata.getData(),
                                           nr_channels,C_linear,S_linear,C_const,S_const);
        // Find optimal 1D partition
        FindBest1DPartition(stripe_udata.getData(),stripe_adata.getData(),stripe_bdata.getData(),stripe_length,
                                     nr_channels,gamma_s,eta_s,Eps1R,C_linear,S_linear,C_const,S_const,L);


        // Get solution from partition and assign to stripe
        mat a_curr = zeros(nr_channels, stripe_length);
        mat b_curr = zeros(nr_channels, stripe_length);
        mat u_curr = zeros(nr_channels, stripe_length);
        ReconstructionFromPartition(L,stripe_udata.getData(),stripe_adata.getData(),stripe_bdata.getData(),
                         stripe_length,nr_channels,eta_s,C_linear,S_linear,u_curr,a_curr,b_curr);


        // Assign stripe to 2D outputs
        u_out.elem(linear_indices) = u_curr.t();
        a_out.elem(linear_indices) = a_curr.t();
        b_out.elem(linear_indices) = b_curr.t();
    }
}

