#ifndef LINEWISEAFFINEMS_H
#define LINEWISEAFFINEMS_H

#include <omp.h>
#include <vector>
#include <list>
#include <iterator>
#define ARMA_NO_DEBUG
#include <armadillo>

#include "Stripe.h"

using namespace std;
using namespace arma;

// Converter from pointers to armadillo objects
void ArmadilloConverter(double* u_data_raw, double* a_data_raw,double* b_data_raw,
                        double* C_linear_raw,double* S_linear_raw,
                        double* C_const_raw, double* S_const_raw,
                        double* u_out_raw, double* a_out_raw, double* b_out_raw,
                        double* dir_raw, const double gamma_s,const double eta_s,
                        const int m, const int n, const int nr_channels, const int nr_threads);

// Applies univariate affine Jet-partitioning to all lines of the input image for given direction
void LinewisePartitioning(cube &u_out,cube &a_out,cube &b_out,
                          const int m,const int n,const int nr_channels,
                          cube &u_data,cube &a_data,cube &b_data,vec &dir, double gamma_s,
                          double eta_s,mat &C_linear, mat &S_linear, mat &C_const, mat &S_const);

// Extracts and stores the stripes of the input image
void Extract1Dstripes(const cube &I, const vec &dir,std::vector<Stripe>  &L, const int m, const int n,const int nr_channels);

// Extracts linear indices of a line of the image domain
uvec GetIndexes(int x_lim,int x_cor,int x_dir,int y_lim,int y_cor,int y_dir);

// Generates the "regression"-matrices for alpha,delta in eq. (27) in the affine-linear 1D-jet estimation 
void GenerateSystemMatrices(const int n,double eta, mat &A);

// Computes and stores the approximation errors for intervals [1,r] for all r
vec Compute1rErrors(mat u_data, mat a_data,mat b_data, const int nr_channels,
                    const mat &C_linear, const mat &S_linear, const mat &C_const, const mat &S_const);

// Computes the optimal univariate partitioning for data f and slope data x,y
void FindBest1DPartition(mat u_data,mat a_data,mat b_data,const int n, const int nr_channels, double &gamma,
                         double eta, const vec &eps_1r, mat &C_linear, mat &S_linear,
                         mat &C_const, mat &S_const, ivec &L);

// Computes the corresponding reconstruction for an optimal partition
void ReconstructionFromPartition(const ivec &L, mat u_data, mat a_data, mat b_data, const int n,
                                 const int nr_channels, double eta,
                                 const mat &C_linear, const mat &S_linear,
                                 mat &u_out, mat &a_out, mat &b_out);

#endif  