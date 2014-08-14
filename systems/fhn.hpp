#ifndef _FHN_SYSTEM_H
#define _FHN_SYSTEM_H
#include "nld_sys.h"

class FHN : public System {
private:
    double a;
    double e;
public:
    FHN();
    virtual void rhs(const state_t &state, state_t & out, double time);
    virtual void jac(const state_t &state,
                     matrix_t & out, double time, state_t &dfdt);
};


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
    virtual void rhs(const state_t &state, state_t & out, double time);
    virtual void jac(const state_t &state,
                     matrix_t & out, double time, state_t &dfdt);
};


#endif
