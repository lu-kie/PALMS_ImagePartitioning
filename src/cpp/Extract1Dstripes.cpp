/**
    Extract1Dstripes.cpp
    Purpose: Extract the stripes of the input image for input direction and 
             stores them as Stripe objects in a vector

    @author Lukas Kiefer
    @version 1.0
*/

#include "linewiseAffineMS.h"
void AssignDataAndIndices(uvec &indices_flat, const int m, const int n, const int nr_channels, const cube &I,std::vector<Stripe>  &L);

void Extract1Dstripes(const cube &I, const vec &dir,std::vector<Stripe>  &L, const int m, const int n,const int nr_channels)
{
    int x_dir = dir(0);
    int y_dir = dir(1);
    if (abs(x_dir) > n || abs(y_dir) > m || x_dir < 0 || ((x_dir == 0) && (y_dir<=0)) ) {
        printf("Error: permitted direction");
        return;
    }
    int x_cor, y_cor; // Starting indexes

    if (y_dir >= 0) {
        // Startpoints: Rows
        for(int row = 0; row < y_dir; row++) {
            for(int j =0; j< n; j++) {
                x_cor = j;
                y_cor = row;
                uvec indices_flat = GetIndexes(n, x_cor, x_dir, m,y_cor,y_dir);
                // Append the data and the linear indices to the vector L
                AssignDataAndIndices(indices_flat,m,n,nr_channels,I,L);
            }
        }
        // Startpoints: Cols
        for(int col = 0; col < x_dir; col++) {
            for(int i= y_dir; i < m ; i++) {
                x_cor = col;
                y_cor = i;
                uvec indices_flat = GetIndexes(n, x_cor, x_dir, m,y_cor,y_dir);
                // Append the data and the linear indices to the vector L
                AssignDataAndIndices(indices_flat,m,n,nr_channels,I,L);
            }
        }
    } else {
        // Startpoints: Rows
        for(int row = m-1; row >= m+y_dir; row--) {
            for(int j =0; j< n; j++) {
                x_cor = j;
                y_cor = row;
                uvec indices_flat = GetIndexes(n, x_cor, x_dir, m,y_cor,y_dir);
                // Append the data and the linear indices to the vector L
                AssignDataAndIndices(indices_flat,m,n,nr_channels,I,L); 
            }
        }
        // Startpoints: Cols
        for(int col = 0; col < x_dir; col++) {
            for(int i= 0; i < m+y_dir; i++) {
                x_cor = col;
                y_cor = i;
                uvec indices_flat = GetIndexes(n, x_cor, x_dir, m,y_cor,y_dir);
                // Append the data and the linear indices to the vector L
                AssignDataAndIndices(indices_flat,m,n,nr_channels,I,L);
            }
        }
    }
}


void AssignDataAndIndices(uvec &indices_flat, const int m, const int n, const int nr_channels, const cube &I,std::vector<Stripe>  &L){
    int t = indices_flat.n_elem;
    uvec indices = uvec(nr_channels * t);
    if(nr_channels == 1)
        indices = indices_flat;
    else {
        for(int q=0; q<nr_channels; q++) {
            indices.rows(q*t , (q+1)*t-1) = indices_flat + (q*(n*m));
        }
    }
    mat data_new = I.elem(indices);
    data_new.reshape(data_new.n_elem/nr_channels,nr_channels);
    inplace_trans(data_new);

    L.push_back(Stripe(data_new,indices));
    return;
}