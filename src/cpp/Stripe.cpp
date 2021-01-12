#include "Stripe.h"

// Getter
mat Stripe::getData(){
    return data;
}

uvec Stripe::getIndices(){
    return indices;
}
// Give Length
int Stripe::giveLength(){
    return data.n_cols;
}

// Setter
void Stripe::setData(mat data_new){
    data = data_new;
}

void Stripe::setIndices(uvec indices_new){
    indices = indices_new;
}
// Destructor
Stripe::~Stripe(){
    
}
// Constructor
Stripe::Stripe(mat data_new, uvec indices_new){
    data = data_new;
    indices = indices_new;
}

Stripe::Stripe(){}