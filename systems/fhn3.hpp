#ifndef _FHN3_SYSTEM_H
#define _FHN3_SYSTEM_H
#include "nld_sys.h"


class FHN3 : public System {
private:
    static const int num = 3;
    double a[num];
    double e[num];
    double Doo;
    double Deo;
    double Doe;
public:
    FHN3();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     state_t & out, double time);
};

#endif
