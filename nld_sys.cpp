#include "nld_sys.h"
#include <assert.h>

System::System(int dimension) {
    dim = dimension;
    vars.resize(dimension);
    for (int i=0; i<4; i++)
        _rk4[i].resize(dimension);
    jacobian.resize(dim);
    for (auto &jac_i : jacobian) jac_i.resize(dim);
    solveBoth = false;
}

void System::setSolveBoth(bool value) {
    int d = dim;
    if (value) {
        d = dim*(dim+1);
    }
    vars.resize(d);
    for (int i=0; i<4; i++)
        _rk4[i].resize(d);

    solveBoth = value;
}

void System::rk4_step(state_t &v, double dt, double time) {
    state_t state(v);
    int len = v.size();
    rhs_dt(state, _rk4[0], time, dt);

    // state = v + _rk4[0]/2.0;
    sum_2vec(v, _rk4[0], 0.5, state);
    rhs_dt(state, _rk4[1], time, dt);

    // state = v + _rk4[1]/2.0;
    sum_2vec(v, _rk4[1], 0.5, state);
    rhs_dt(state, _rk4[2], time, dt);

    // state = v + _rk4[2];
    sum_2vec(v, _rk4[0], state);
    rhs_dt(state, _rk4[3], time, dt);

    for (int i=0; i<len; i++) {
        v[i] += (_rk4[0][i] + _rk4[1][i]*2.0 + _rk4[2][i]*2.0 + _rk4[3][i])/6.0;
    }

}

TimeSeries & System::getTs() {
    return ts;
}

void System::solve(double MaxTime, double dt,
                   state_t ini,
                   double saveTS) {
    ts.clear();
    int ts_point_counter = 0;
    double time = 0;
    vars.swap(ini);
    while (time < MaxTime) {
        this->rk4_step(vars, dt, time);
        if (time > ts_point_counter*saveTS) {
            std::vector<double> p;
            p.push_back(time);
            for (int i=0; i<dim; i++) {
                p.push_back(vars[i]);
            }
            ts.addPoint(p);
            ts_point_counter++;
        }
        time += dt;
    }
}

void System::solve(double MaxTime, double dt,
                   state_t ini) {
    double time = 0;
    vars.swap(ini);
    while (time < MaxTime) {
        this->rk4_step(vars, dt, time);
        time += dt;
    }
}


void System::rhs_dt(const state_t &state, state_t & out, double time, double dt) {
    auto s = state.begin();
    auto o = out.begin();

    if (solveBoth){
        rhs(s,o,time);
        auto s_lin = s+dim;
        auto o_lin = o+dim;
        rhs_linearized(s,s_lin,o_lin,time);
    } else {
        rhs(s,o,time);
    }

    for (auto it=out.begin(); it != out.end(); it++)
        (*it) *= dt;
}

void System::rhs_linearized(const_it_t &state, const_it_t& lin_state,
                            it_t& out, double time) {
    jac(state, jacobian, time);
    for (int j=0; j<dim; j++){
        for (int i=0; i<dim; i++) {
            *out = std::inner_product(lin_state+dim*j,
                                      lin_state+dim*(j+1),
                                      jacobian[i].begin(),0.0);
            out++;
        }
    }
}
