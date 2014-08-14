#ifndef _NLD_SYS
#define _NLD_SYS
#include "utils.hpp"
#include "analyzers/analyzer.hpp"

typedef boost::function<void (const state_t&, state_t &, const double)> rhs_fn_t;
typedef boost::function<void (const state_t&, matrix_t &, const double, state_t &dfdt)> jac_fn_t;

class System {
protected:
    int dim;
    state_t vars;
    state_t dfdt;
    state_t _rk4[4];
    std::vector<Analyzer*> analyzers;
    matrix_t jacobian;
    ode::runge_kutta4<state_t> rk4_stepper;
    ode::rosenbrock4<double> rosenbrock4_stepper;
    ode::runge_kutta_fehlberg78<state_t> rk_fehlberg_stepper;
    ode::runge_kutta_cash_karp54_classic<state_t> rk_ck54_stepper;
    int stepper_kind;
    rhs_fn_t rhs_b;
    rhs_fn_t jac_b;
    rhs_fn_t rhs_combined_b;
    rhs_fn_t rhs_implicit_active;
    std::pair<rhs_fn_t,jac_fn_t> ode_explicit_pair;
public:
    enum{
        STEPPER_RK4,
        STEPPER_ROSENBROCK4,
        STEPPER_RK_FEHLBERG78,
        STEPPER_RK_CK54
    };
    void addAnalyzer(Analyzer *analyzer);
    void setStepper(int _stepper);
    System(int dimension);
    virtual void rhs(const state_t & state,
                     state_t & out, const double time) = 0;

    void rhs_combined(const state_t & state,
                         state_t & out, const double time);

    virtual void jac(const state_t & state,
                     matrix_t & out, const double time, state_t &dfdt) = 0;

    void setSolveCombined(bool value);

    void solve(double MaxTime, double dt,
          state_t ini);

    template<class AnalyzerType>
        AnalyzerType* getAnalyzer();
    int getDim() { return dim; }
    state_t & getState() { return vars; }
};

template<class AnalyzerType>
AnalyzerType* System::getAnalyzer() {
    for (auto &a : analyzers) {
        AnalyzerType *ret = dynamic_cast<AnalyzerType*>(a);
        if (ret) return ret;
    }
    return NULL;
}

#endif
