#ifndef _FHN_SYSTEM_H
#define _FHN_SYSTEM_H
#include "nld_sys.h"

class FHN : public System {
private:
    double a;
    double e;
public:
    FHN();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     std::vector<state_t> & out, double time);
};

#endif
