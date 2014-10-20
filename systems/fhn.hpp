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
protected:
    int num;
    double *a;
    double *e;
    double Doo;
    double Deo;
    double Doe;

public:
    FHN3();
    ~FHN3() { free(a); free(e); }
    virtual void rhs(const state_t &state, state_t & out, double time);
    virtual void jac(const state_t &state,
                     matrix_t & out, double time, state_t &dfdt);
    void setDoo(double v) { Doo = v; }
    void setDeo(double v) { Deo = v; }
    void setDoe(double v) { Doe = v; }
    void setEpsilon(int i, double _e) { e[i] = _e; }
};


#endif
