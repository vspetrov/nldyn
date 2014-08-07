#include "nld_sys.h"
#include <stdio.h>
#include <assert.h>

System::System(int dimension) {
    dim = dimension;
    vars.resize(dimension);
    for (int i=0; i<4; i++)
        _rk4[i].resize(dimension);
    jacobian.resize(dim);
    for (auto &jac_i : jacobian) jac_i.resize(dim);
}

void System::rk4_step(double dt, double time, state_t *state_for_jacobian) {
    std::vector<double> state(vars);
    rhs_dt(state, _rk4[0], time, dt, state_for_jacobian);

    state = vars + _rk4[0]/2.0;
    rhs_dt(state, _rk4[1], time, dt, state_for_jacobian);

    state = vars + _rk4[1]/2.0;
    rhs_dt(state, _rk4[2], time, dt, state_for_jacobian);

    state = vars + _rk4[2];
    rhs_dt(state, _rk4[3], time, dt, state_for_jacobian);

    for (int i=0; i<dim; i++) {
        vars[i] += (_rk4[0][i] + _rk4[1][i]*2.0 + _rk4[2][i]*2.0 + _rk4[3][i])/6.0;
    }

}

TimeSeries & System::getTs() {
    return ts;
}

void System::solve(double MaxTime, double dt,
                   std::vector<double> ini,
                   double saveTS) {
    ts.clear();
    int ts_point_counter = 0;
    double time = 0;
    vars.swap(ini);
    while (time < MaxTime) {
        this->rk4_step(dt, time);
        if (saveTS > 0 && time > ts_point_counter*saveTS) {
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

void System::solve_linearized(ts_t &ts, state_t &ini) {

    double dt = ts[1][0]-ts[0][0];
    vars = ini;
    for (auto &ts_point : ts) {
        state_t phase_point(ts_point.begin()+1,ts_point.end());
        double time = ts_point[0];
        rk4_step(dt, time, &phase_point);
    }
}
void System::rhs_dt(std::vector<double> &state, std::vector<double> & out, double time, double dt, state_t *state_for_jacobian) {
    if (state_for_jacobian){
        this->rhs_linearized(*state_for_jacobian,state,out,time);
    } else {
        this->rhs(state,out,time);
    }
    for (auto it=out.begin(); it != out.end(); it++)
        (*it) *= dt;
}

void System::rhs_linearized(std::vector<double> &state, state_t &lin_state, std::vector<double> &dvdt, double time) {
    jac(state, jacobian, time);
    for (int i=0; i<dim; i++)
        dvdt[i] = dot(lin_state, jacobian[i]);
}


Lorenz::Lorenz(int dimension) : System(dimension) {
    this->sigma = 10.0;
    this->beta  = 8./3.;
    this->ro    = 28.0;
}

void Lorenz::rhs(std::vector<double> &state, std::vector<double> &dvdt, double time) {
    dvdt[0] = sigma*(state[1]-state[0]);
    dvdt[1] = state[0]*(ro - state[2]) - state[1];
    dvdt[2] = state[0]*state[1]-beta*state[2];
}


void Lorenz::jac(std::vector<double> &state,
                 std::vector<std::vector<double> >&out, double time) {
    assert(out.size() == dim);

    out[0][0] = -sigma;      out[0][1] = sigma;    out[0][2] = 0;
    out[1][0] = ro-state[2]; out[1][1] = -1;       out[1][2] = -state[0];
    out[2][0] = state[1];    out[2][1] = state[0]; out[2][2] = -beta;
}


Rossler::Rossler() : System(3) {
    a = 0.1;
    b = 0.1;
    c = 18.0;
}

void Rossler::rhs(std::vector<double> &state, std::vector<double> &dvdt, double time) {
    dvdt[0] = -state[1]-state[2];
    dvdt[1] = state[0]+a*state[1];
    dvdt[2] = b+state[2]*(state[0]-c);
}


void Rossler::jac(std::vector<double> &state,
                 std::vector<std::vector<double> >&out, double time) {
    assert(out.size() == dim);

    out[0][0] = 0;      out[0][1] = -1;    out[0][2] = -1;
    out[1][0] = 1; out[1][1] = a;       out[1][2] = 0;
    out[2][0] = state[2];    out[2][1] = 0; out[2][2] = state[0]-c;
}


FHN::FHN() : System(2) {
    a=0;
    e=0.01;
}

void FHN::rhs(std::vector<double> &state, std::vector<double> &dvdt, double time) {
    dvdt[0] = state[0]-state[0]*state[0]*state[0]/3.0-state[1];
    dvdt[1] = e*(state[0]-a);
}


void FHN::jac(std::vector<double> &state,
                 std::vector<std::vector<double> >&out, double time) {
    assert(out.size() == dim);

    out[0][0] = 1-state[0]*state[0];      out[0][1] = -1;
    out[1][0] = e; out[1][1] = 0;
}


TimeSeries::TimeSeries() {
}

void TimeSeries::addPoint(std::vector<double> &p) {
    ts.push_back(p);
}

std::vector<std::vector<double> >& TimeSeries::getAll() {
    return ts;
}
std::vector<double> TimeSeries::row(int i) {
    std::vector<double> row;
    for (auto point : ts) row.push_back(point[i]);
    return row;
}

void TimeSeries::plot(int x, int y) {
    auto xr = row(x);
    auto yr = row(y);
    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp,"plot '-' w l\n");
    assert(xr.size() == yr.size());
    for (int i=0; i<xr.size(); i++) {
        fprintf(gp,"%g %g\n",xr[i],yr[i]);
    }
    // fprintf(gp,"e");
    fclose(gp);
}

TimeSeries::~TimeSeries() {
}
