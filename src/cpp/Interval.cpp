#include "Interval.h"

// Getter
int Interval::getL(){
    return l;
}
int Interval::getR(){
    return r;
}
double Interval::getEps(){
    return eps;
}
mat Interval::getUdata(){
    return u_data;
}
mat Interval::getAdata(){
    return a_data;
}
mat Interval::getBdata(){
    return b_data;
}

// Setter
void Interval::setEps(double e){
    eps = e;
}
void Interval::setL(int left){
    l = left;
}
void Interval::setR(int right){
    r = right;
}
void Interval::setUdata(mat y){
    u_data = b_data;
}
void Interval::setAdata(mat z){
    a_data = z;
}
void Interval::setBdata(mat w){
    b_data = w;
}
// Destructor
Interval::~Interval(){
}
// Constructor
Interval::Interval(int left,int right,const int nr_channels,
        mat g, mat z, mat w){
            
    l = left;
    r = right;
    eps = 0.0;
    u_data = g;
    a_data = z;
    b_data = w;
}


// Give interval / data length
int Interval::giveLength() const
{
    return r-l+1;
}
// Add data point to the bottom
void Interval::addBottomDataPoint(const int nr_channels,
        const mat &C_linear, const mat &S_linear,
        const mat &C_const, const mat &S_const,
        vec udata_new, vec adata_new, vec bdata_new){
    
    double c,s, hj_old, fr_old, xr_old, yj_old,yr_old;
    int h = giveLength(); //h is the current interval length

    // u_data,a_data
    // Handle special case that OLD interval length is 1
    if (h == 1) {
        c = C_linear(1,0);
        s = S_linear(1,0);
        vec f1_old = u_data.col(0);
        vec x1_old = a_data.col(0);
        u_data.col(0) = c*f1_old + s*x1_old;
        a_data.col(0) = -s*f1_old + c*x1_old;
    }
    // Eliminate new row 1
    for(int j=0; j < 2; j++) {
        c = C_linear(2*h,j);
        s = S_linear(2*h,j);
        if (j==0) {
            for(int q = 0; q < nr_channels; q++) { // Faster than vectorized
                hj_old = u_data(q,0);
                fr_old = udata_new(q);

                u_data(q,0) = c*hj_old + s*fr_old;
                udata_new(q) = -s*hj_old + c*fr_old;
            }
        }
        if (j==1) {
            for(int q = 0; q < nr_channels; q++) {
                hj_old = a_data(q,0);
                fr_old = udata_new(q);

                a_data(q,0) = c*hj_old + s*fr_old;
                udata_new(q) = -s*hj_old + c*fr_old;
            }
        }
    }
    // Eliminate new row 2
    for(int j=0; j < 2; j++) {
        c = C_linear(2*h+1,j);
        s = S_linear(2*h+1,j);
        if (j==0) {
            for(int q = 0; q < nr_channels; q++) {
                hj_old = u_data(q,0);
                xr_old = adata_new(q);

                u_data(q,0) = c*hj_old + s*xr_old;
                adata_new(q) = -s*hj_old + c*xr_old;
            }
        }
        if(j==1) {
            for(int q = 0; q < nr_channels; q++) {
                hj_old = a_data(q,0);
                xr_old = adata_new(q);

                a_data(q,0) = c*hj_old + s*xr_old;
                adata_new(q) = -s*hj_old + c*xr_old;
            }
        }

    }

    // b_data
    // Look up Givens rotation coefficients
    c = C_const(h,0);
    s = S_const(h,0);
    // Eliminate new row 3
    for(int q = 0; q < nr_channels; q++) {
        yj_old = b_data(q,0);
        yr_old = bdata_new(q,0);
        b_data(q,0) = c*yj_old + s*yr_old;
        bdata_new(q,0) = -s*yj_old + c*yr_old;
    }

    // Update interval error, data, length
    if(h+1 >= 2) { // h is old length
        for(int q=0; q < nr_channels; q++)
            eps+= (udata_new(q)*udata_new(q)) + (adata_new(q)*adata_new(q)) + (bdata_new(q)*bdata_new(q));
    }

    r++;
    if(h < 2) {
        u_data.col(h) = udata_new;
        a_data.col(h) = adata_new;
    }
}

// Update associated data i.e. sparse Givens rotate it (in reconstruction process)
void Interval::givensRotateLinearData(int w, const mat &C_linear, const mat &S_linear){
    
    double c,s;
    // mixed linear constant data
    // Handle first rotation separately
    if(w == 0) {
        c = C_linear(1,0);
        s = S_linear(1,0);
        vec f1_old = u_data.col(0);
        vec x1_old = a_data.col(0);
        u_data.col(0) = c*f1_old + s*x1_old;
        a_data.col(0) = -s*f1_old + c*x1_old;

    } else {
        if(w % 2 == 0) {
            // Do f
            c = C_linear(w,0);
            s = S_linear(w,0);
            vec f1_old = u_data.col(0);
            vec fw_old = u_data.col(w/2);
            u_data.col(0) = c*f1_old + s*fw_old;
            u_data.col(w/2) = -s*f1_old + c*fw_old;

            c = C_linear(w,1);
            s = S_linear(w,1);
            fw_old = u_data.col(w/2);
            vec x1_old = a_data.col(0);
            a_data.col(0) = c*x1_old + s*fw_old;
            u_data.col(w/2) = -s*x1_old + c*fw_old;
        } else {
            // Do x
            c = C_linear(w,0);
            s = S_linear(w,0);
            vec f1_old = u_data.col(0);
            vec xw_old = a_data.col((w-1)/2);
            u_data.col(0) = c*f1_old + s*xw_old;
            a_data.col((w-1)/2) = -s*f1_old + c*xw_old;

            c = C_linear(w,1);
            s = S_linear(w,1);
            xw_old = a_data.col((w-1)/2);
            vec x1_old = a_data.col(0);
            a_data.col(0) = c*x1_old + s*xw_old;
            a_data.col((w-1)/2) = -s*x1_old + c*xw_old;
        }
    }
}