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
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     std::vector<state_t> & out, double time);
};

#endif
