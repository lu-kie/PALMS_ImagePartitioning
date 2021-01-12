/**
    GenerateSystemMatrices.cpp
    Purpose: Generates the system-matrices for affine-linear 1D-jet estimation

    @author Lukas Kiefer
    @version 1.0
*/

#include "linewiseAffineMS.h"

void GenerateSystemMatrices(const int n,double eta, mat &A)
{
    // Create system matrix for linear part in A_q in eq. (27)
    for(unsigned int i=0; i < A.n_rows; i++ ){
        A(i,0) = eta*(i/2 + 1);
        A(i,1) = eta*1;
        i++;
        A(i,0) = 1;
        A(i,1) = 0;
    }
}
