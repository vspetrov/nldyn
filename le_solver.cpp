#include "le_solver.h"

std::vector<double> calcLE(double warmUpTime,
                           int numSteps,
                           double stepTime,
                           double dt,
                           state_t &ini) {
    system->solve(warmUpTime, dt, ini);
    state_t s0 = system->getState();
}
