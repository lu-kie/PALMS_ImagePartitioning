/**
    ReconstructionFromPartition.cpp
    Purpose: Computes the corresponding 1D-reconstruction of the input 1D-partition
    
    @author Lukas Kiefer
    @version 1.0
*/

#include "Interval.h"
#include "linewiseAffineMS.h"

bool inline compare_intLengths(const Interval &inter1, const Interval &inter2)
{
    return (inter1.giveLength() < inter2.giveLength());
}

void ReconstructionFromPartition(const ivec &L, mat u_data, mat a_data, mat b_data, const int n,
                                 const int nr_channels, double eta,
                                 const mat &C_linear, const mat &S_linear,
                                 mat &u_out, mat &a_out, mat &b_out){
    // Get system matrices of underlying least squares problem
    mat A = zeros(2*n,2);
    // A holds the columns 1 and 3 from A_q of eq. (27)
    GenerateSystemMatrices(n,eta,A);
    // Create interval objects
    std::list<Interval> Intervals;
    int r = n,l;
    while(true) {
        l = L(r-1)+1;
        // Handle/Catch interval lengths < 2
        if (r-l +1 < 2) {
            u_out.cols(l-1,r-1) = u_data.cols(l-1,r-1);
            a_out.cols(l-1,r-1) = a_data.cols(l-1,r-1);
            vec b_mean(nr_channels);
            for(int ch=0; ch < nr_channels; ch++) {
                vec tmp = b_data(span(ch,ch) , span(l-1,r-1));
                b_mean(ch) = mean(tmp);
            }
            for(int w = l; w <= r; w++ ) {
                b_out.col(w-1) = b_mean;
            }
        }
        else {
            Intervals.push_front(Interval(l,r,nr_channels,eta*u_data.cols(l-1,r-1),a_data.cols(l-1,r-1),b_data.cols(l-1,r-1)));
        }
        if (l==1)
            break;
        r = l-1;
    }
    if(Intervals.empty())
        return;
    // Sort Intervals in length-ascending order
    // (This allows to avoid computing the same row transformations of the system matrix multiple times.)
    Intervals.sort(compare_intLengths);
    // Solve linear equation systems
    mat Aj_old = mat(1,2);
    mat Aw_old = mat(1,2);
    mat p = zeros<mat>(nr_channels,3);

    rowvec I = linspace<vec>(1,n,n).t();
    mat udata_curr,adata_curr;
    double c,s;
    // The master iterator knowing which intervals have to be considered (i.e. which interval lengths)
    list<Interval>::iterator master_it = Intervals.begin();

    // Handle v=0 separately:
    c = C_linear(1,0);
    s = S_linear(1,0);
    double aj_old, aw_old;
    for(int iter = 0; iter < 2; iter++) {
        aj_old =A(0,iter);
        aw_old = A(1,iter);
        A(0,iter) = c*aj_old+s*aw_old;
        A(1,iter) = -s*aj_old+c*aw_old;
    }
    // Update data
    for(list<Interval>::iterator it = master_it; it!= Intervals.end(); ++it) {
        Interval &currentInterval = *it;
        currentInterval.givensRotateLinearData(0,C_linear,S_linear);
    }

    // v > 0:
    for(int v = 1; v < n; v++) {
        // Linear rows
        int w = 2*v;
        for(int j = 0; j < 2; j++) {
            // Read-in Givens rotation coefficients
            c = C_linear(w,j);
            s = S_linear(w,j);
            // Update A
            for(int iter=0; iter < 2; iter++) {
                aj_old = A(j,iter);
                aw_old = A(w,iter);
                A(j,iter) = c*aj_old+s*aw_old;
                A(w,iter) = -s*aj_old+c*aw_old;
            }
        }
        // Constant row
        w++;
        for(int j = 0; j < 2; j++) {
            c = C_linear(w,j);
            s = S_linear(w,j);
            // Update A
            for(int iter=0; iter < 2; iter++) {
                aj_old = A(j,iter);
                aw_old = A(w,iter);
                A(j,iter) = c*aj_old+s*aw_old;
                A(w,iter) = -s*aj_old+c*aw_old;
            }

        }
        // Update data
        for(list<Interval>::iterator it = master_it; it!= Intervals.end(); ++it) {
            Interval &currentInterval = *it;
            // Linear  data
            w--;
            currentInterval.givensRotateLinearData(w,C_linear,S_linear);
            w++;
            currentInterval.givensRotateLinearData(w,C_linear,S_linear);
        }
        //
        while(true) {
            Interval &curr = *master_it;
            l = curr.getL();
            r = curr.getR();
            if(v!=curr.giveLength()-1){
                break; 
            }
            // Back substitution
            p.zeros();
            udata_curr = curr.getUdata();
            adata_curr = curr.getAdata();

            // Get linear coefficients
            for(int ch = 0; ch < nr_channels; ch++) {
                p(ch,2)  = adata_curr(ch,0)/A(1,1); //offset for origin in left interval boarder
                p(ch,1)  = (udata_curr(ch,0) - A(0,1)*p(ch,2))/A(0,0); //slope a
                p(ch,0)  = mean(curr.getBdata().row(ch)); // slope b
            }

            // Fill segment channelwise
            for(int ch=0; ch<nr_channels; ch++){
                // compute the functional values on the interval
                u_out(span(ch,ch),span(curr.getL()-1,curr.getR()-1)) = p(ch,1)*I.cols(0,curr.giveLength()-1)+p(ch,2);
            }
            // Save slopes
            for(int t = l; t <= r; t++) {
                a_out.col(t-1) = p.col(1);
                b_out.col(t-1) = p.col(0);
            }
            //Increase master Iterator (a filled segment won't be considered again)
            ++master_it;
            //Check if master iterator has finished, i.e., if reconstructions on all intervals are finished
            if(master_it == Intervals.end()) {
                return;
            }
        }
    }
    return;
}
