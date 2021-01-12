#ifndef INTERVAL_H
#define INTERVAL_H
#define ARMA_NO_DEBUG
#include <armadillo>
using namespace arma;
using namespace std;

class Interval
{
private:
    int l; //Left bound of discrete interval
    int r; //Right bound of discrete interval
    double eps; //Approximation error for interval
    mat u_data; // Interval data
    mat a_data; // Slope data
    mat b_data; // Perpendicular slope data
public:
    // Getter:
    int getL();
    int getR();
    double getEps();
    mat getUdata();
    mat getAdata();
    mat getBdata();
    // Setter:
    void setEps(double e);
    void setL(int left);
    void setR(int right);
    void setUdata(mat g);
    void setAdata(mat z);
    void setBdata(mat w);
    // Destructor
    ~Interval();
    // Constructor
    Interval(int left,int right,const int nr_channels,mat g,mat z,mat w);
    // Give interval / data length
    int giveLength() const;
    // Add data point to the bottom of the interval
    void addBottomDataPoint(const int nr_channels,
                            const mat &C_linear,const mat &S_linear,
                            const mat &C_const, const mat &S_const,
                            vec udata_new, vec adata_new, vec bdata_new);
    // Update associated data i.e. sparse Givens rotate it (for reconstruction process)
    void givensRotateLinearData(int w, const mat &C_linear, const mat &S_linear);
    void givensRotateConstData(int v, mat &C_const, mat &S_const);

};
#endif
