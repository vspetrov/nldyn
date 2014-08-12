#ifndef _AP_SYSTEM_H
#define _AP_SYSTEM_H
#include "nld_sys.h"

class AP : public System {
private:
    double a;
    double k;
    double ksi;
public:
    AP();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     state_t & out, double time);
};

class AP3 : public System {
private:
    static const int num = 3;
    double a[num];
    double k[num];
    double ksi[num];
    double Doo;
    double Deo;
    double Doe;
public:
    AP3();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     state_t & out, double time);
};

#endif
