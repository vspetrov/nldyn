#ifndef _LE_SOLVER_H
#define _LE_SOLVER_H
#include "nld_sys.h"

class LyapunovExpsSolver {
public:
    LyapunovExpsSolver(System *s) {
        nld_sys = s;
        debugFlag = false;
    }
    std::vector<double> calcLE(double warmUpTime,
                               double wudt,
                               int numSteps,
                               double stepTime,
                               double dt,
                               std::vector<double> &ini,
                               double eps = -1);

    double getKYdim() { return KYdim; }
    void setDbg(bool value) { debugFlag = value; }
private:
    bool debugFlag;
    void calcKaplanYorkeDimension(std::vector<double> &LEs);
    std::vector<double> LEs;
    double KYdim;
    System *nld_sys;
};

#endif
