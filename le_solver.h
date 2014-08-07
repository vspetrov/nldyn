#ifndef _LE_SOLVER_H
#define _LE_SOLVER_H
#include "nld_sys.h"

class LyapunovExpsSolver {
public:
    LyapunovExpsSolver(System *s) { nld_sys = s; }
    std::vector<double> calcLE(double warmUpTime,
                               double wudt,
                               int numSteps,
                               double stepTime,
                               double dt,
                               std::vector<double> &ini);
private:
    System *nld_sys;
};

#endif
