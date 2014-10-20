#ifndef _LE_SOLVER_H
#define _LE_SOLVER_H
#include "nld_sys.h"

class LyapunovExpsSolver {
public:
    LyapunovExpsSolver(System *s);
    std::vector<double> & calcLE(double warmUpTime,
                     double wudt,
                     int numSteps,
                     double stepTime,
                     double dt,
                     state_t &ini,
                     double eps = -1);

    double getKYdim() { return KYdim; }
    void setKeepImmLEs(bool value) { keepImmediateLEs = value; }
    void plotFiniteTimeLEs(double windowTime);
private:
    int dim;
    void GramShmidt(ublas::vector_range<state_t> &state, std::vector<double> &norms);
    std::vector<double> projections;
    bool keepImmediateLEs;
    void calcKaplanYorkeDimension();
    std::vector<double> LEs;
    std::vector<std::vector<double> > m_ImmLes;
    double KYdim;
    System *nld_sys;
};

#endif
