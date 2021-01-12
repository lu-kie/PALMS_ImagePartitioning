#ifndef STRIPE_H
#define STRIPE_H

#define ARMA_NO_DEBUG
#include <armadillo>

using namespace arma;

class Stripe
{
private:
    mat data;
    uvec indices; //for using .elem, i.e. linear ("idx" in matlab) indices
public:
    //Getter
    mat getData();
    uvec getIndices();
    //Setter
    void setData(mat data_new);
    void setIndices(uvec indices_new);
    //Destructor
    ~Stripe();
    //Constructor
    Stripe(mat data_new,uvec indices_new);
    Stripe();
    //Give length
    int giveLength();
};

#endif