#ifndef _NLD_SYS
#define _NLD_SYS
#include "utils.hpp"
#include "time_series.hpp"



class System {
protected:
    int dim;
    state_t vars;
    state_t _rk4[4];
    TimeSeries ts;
    state_t jacobian;
    bool solveBoth;
public:
    System(int dimension);
    virtual void rhs(const_it_t &state,
                     it_t & out, double time) = 0;

    void rhs_linearized(const_it_t &state,
                        const_it_t & lin_state,
                        it_t & out, double time);

    virtual void jac(const_it_t &state,
                     state_t &out, double time) = 0;
    void rk4_step(state_t &v, double dt, double time);
    void rhs_dt(const state_t &state, state_t & out, double time, double dt);
    void solve(double MaxTime, double dt,
          state_t ini,
          double saveTS);
    void solve(double MaxTime, double dt, state_t ini);
    TimeSeries & getTs();
    int getDim() { return dim; }
    state_t & getState() { return vars; }
    void setSolveBoth(bool value);
    int getDimBoth() { return dim*(dim+1);}
};

#endif
