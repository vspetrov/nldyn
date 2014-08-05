#ifndef _LE_SOLVER_H
#define _LE_SOLVER_H
#include "nld_sys.h"

class LyapunovExpsSolver {
public:
    LyapunovExpsSolver(System *s) { system = s; }
    std::vector<double> calcLE(double warmUpTime,
                               int numSteps,
                               double stepTime,
                               double dt,
                               std::vector<double> &ini);
private:
    System *system;
};

#endif
