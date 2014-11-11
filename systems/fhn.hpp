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
    double Doo;
    double Deo;
    double Doe;
    double *a;
    double *e;

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

class FHN_O2_E_chain : public System {
private:
    int num;
    double Doo;
    double Deo;
    double Doe;
    double *a;
    double *e;
    double Dee;
public:
    FHN_O2_E_chain(int _N);
    ~FHN_O2_E_chain() { delete a; delete e; }
    virtual void rhs(const state_t &state, state_t & out, double time);
    virtual void jac(const state_t &state,
                     matrix_t & out, double time, state_t &dfdt);
    void setDee(double v) { Dee = v; }
    void setDoo(double v) { Doo = v; }
    void setDeo(double v) { Deo = v; }
    void setDoe(double v) { Doe = v; }
    void setEpsilon(int i, double _e) { e[i] = _e; }

};

#endif
