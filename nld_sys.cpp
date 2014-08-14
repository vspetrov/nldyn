#include "nld_sys.h"
#include <assert.h>
#include <boost/bind.hpp>

System::System(int dimension) {
    dim = dimension;
    vars.resize(dimension);
    dfdt.resize(dim);
    for (int i=0; i<4; i++)
        _rk4[i].resize(dimension);
    jacobian = matrix_t(dim, dim);
    ode_explicit_pair = std::make_pair ([this](const state_t & state,
                                               state_t &out, const double time) -> void {
                                            this->rhs(state, out, time);
                                        },[this](const state_t & state,
                                                 matrix_t &out, const double time, state_t &dfdt) -> void {
                                            this->jac(state, out, time, dfdt);
                                        });
    stepper_kind = STEPPER_RK4;
    rhs_b = boost::bind(&System::rhs, this, _1, _2, _3);
    rhs_combined_b = boost::bind(&System::rhs_combined, this, _1, _2, _3);
    rhs_implicit_active = rhs_b;
}

void System::addAnalyzer(Analyzer *analyzer) {
    analyzers.push_back(analyzer);
}


void System::setSolveCombined(bool value) {
    if (value) {
        rhs_implicit_active = rhs_combined_b;
    } else {
        rhs_implicit_active = rhs_b;
    }
    rk4_stepper = ode::runge_kutta4<state_t>();
    rk_fehlberg_stepper = ode::runge_kutta_fehlberg78<state_t>();
    // rosenbrock4_stepper = ode::rosenbrock4<double>();
    rk_ck54_stepper = ode::runge_kutta_cash_karp54_classic<state_t>();
}

void System::solve(double MaxTime, double dt,
                   state_t ini) {

    double time = 0;
    vars.swap(ini);


    while (time < MaxTime) {
        switch(stepper_kind){
        case STEPPER_RK4:
            rk4_stepper.do_step(rhs_implicit_active, vars, time, dt);
            break;
        case STEPPER_RK_FEHLBERG78:
            rk_fehlberg_stepper.do_step(rhs_implicit_active, vars, time, dt);
            break;
        case STEPPER_ROSENBROCK4:
            rosenbrock4_stepper.do_step(ode_explicit_pair, vars, time, dt);
            break;
        case STEPPER_RK_CK54:
            rk_ck54_stepper.do_step(rhs_implicit_active, vars, time, dt);
            break;
        }
        for (auto &a : analyzers) {
            a->addPoint(vars, time);
        }
        time += dt;
    }
}

void System::setStepper(int _stepper) {
    stepper_kind = _stepper;
}


void System::rhs_combined(const state_t &state,
                            state_t & out, double time) {
    assert(out.size() == dim*(dim+1));
    rhs(state, out, time);
    jac(ublas::subrange(state,0,dim), jacobian, time, dfdt);
    for (int i=1; i<out.size() / dim ; i++) {
        auto st = ublas::subrange(state,dim*i,dim*(i+1));
        ublas::subrange(out, dim*i, dim*(i+1)) = ublas::prod(jacobian, st);
    }
}
