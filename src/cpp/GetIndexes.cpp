/**
    GetIndexes.cpp
    Purpose: Extracts linear indices of a line of the image domain
             starting index: (x_cor,y_cor)
             end index: (x_lim,y_lim)
             direction: (x_dir,y_dir)

    @author Lukas Kiefer
    @version 1.0
*/

#include "linewiseAffineMS.h"

uvec GetIndexes(int x_lim,int x_cor,int x_dir,int y_lim,int y_cor,int y_dir){
    
    // Determine the number of steps from (x_cor,y_cor) to (x_lim,y_lim) in terms of the direction (x_dir,y_dir)
    int nr_steps;
    
    if (x_dir == 0)
        nr_steps = ((y_lim-y_cor)-1)/y_dir;
    else if (y_dir > 0)
        nr_steps = min( (x_lim-x_cor-1)/x_dir , (y_lim-y_cor-1)/y_dir);
    else if (y_dir < 0)
        nr_steps = min( (x_lim-x_cor-1)/x_dir, -y_cor/y_dir);
    else if (y_dir == 0)
        nr_steps = (x_lim-x_cor-1)/x_dir;
        
    // Get the indexes as subscripts of the image domain
    umat indexes = umat(2,nr_steps+1);

    indexes.row(1) = conv_to<umat>::from(linspace(x_cor,x_cor+nr_steps*x_dir, nr_steps+1)).t();//x_indexes
    indexes.row(0) = conv_to<umat>::from(linspace(y_cor,y_cor+nr_steps*y_dir, nr_steps+1)).t();//y_indexes
    // Convert the subscripts to linear indexes of the image domain
    uvec indices(nr_steps+1,1);
    for(int i = 0; i < nr_steps+1; i++) {
        indices(i) = sub2ind( size(y_lim,x_lim), indexes(0,i),indexes(1,i) );
    }

    return indices;
}
