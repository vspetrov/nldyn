#ifndef _LORENZ_SYSTEM_H
#define _LORENZ_SYSTEM_H
#include "nld_sys.h"

class Lorenz : public System {
private:
    double sigma;
    double ro;
    double beta;
public:
    Lorenz();
    virtual void rhs(const state_t &state, state_t & out, double time);
    virtual void jac(const state_t &state, matrix_t & out, double time, state_t &dfdt);
};

#endif
