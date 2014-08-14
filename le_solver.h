#ifndef _LE_SOLVER_H
#define _LE_SOLVER_H
#include "nld_sys.h"

class LyapunovExpsSolver {
public:
    LyapunovExpsSolver(System *s) {
        nld_sys = s;
        debugFlag = false;
        dim = s->getDim();
        projections.resize(dim);
    }
    std::vector<double> & calcLE(double warmUpTime,
                     double wudt,
                     int numSteps,
                     double stepTime,
                     double dt,
                     state_t &ini,
                     double eps = -1);

    double getKYdim() { return KYdim; }
    void setDbg(bool value) { debugFlag = value; }

private:
    int dim;
    void GramShmidt(ublas::vector_range<state_t> &state, std::vector<double> &norms);
    std::vector<double> projections;
    bool debugFlag;
    void calcKaplanYorkeDimension();
    std::vector<double> LEs;
    double KYdim;
    System *nld_sys;
};

#endif
