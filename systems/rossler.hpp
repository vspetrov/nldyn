#ifndef _ROSSLER_SYSTEM_H
#define _ROSSLER_SYSTEM_H
#include "nld_sys.h"

class Rossler : public System {
private:
    double a;
    double b;
    double c;
public:
    Rossler();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     std::vector<state_t> & out, double time);
};

#endif
