#include "nld_sys.h"
#include <stdio.h>
#include <assert.h>
#include <string>
#include <sstream>

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


Lorenz::Lorenz() : System(3) {
    this->sigma = 10.0;
    this->beta  = 8./3.;
    this->ro    = 28.0;
}

void Lorenz::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = sigma*(state[1]-state[0]);
    dvdt[1] = state[0]*(ro - state[2]) - state[1];
    dvdt[2] = state[0]*state[1]-beta*state[2];
}


void Lorenz::jac(const_it_t &state,
                 std::vector<state_t> & out, double time) {
    assert(out.size() == dim);

    out[0][0] = -sigma;      out[0][1] = sigma;    out[0][2] = 0;
    out[1][0] = ro-state[2]; out[1][1] = -1;       out[1][2] = -state[0];
    out[2][0] = state[1];    out[2][1] = state[0]; out[2][2] = -beta;
}


Rossler::Rossler() : System(3) {
    a = 0.2;
    b = 0.2;
    c = 5.7;
}

void Rossler::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = -state[1]-state[2];
    dvdt[1] = state[0]+a*state[1];
    dvdt[2] = b+state[2]*(state[0]-c);
}


void Rossler::jac(const_it_t &state,
                  std::vector<state_t> & out, double time) {
    assert(out.size() == dim);

    out[0][0] = 0;      out[0][1] = -1;    out[0][2] = -1;
    out[1][0] = 1; out[1][1] = a;       out[1][2] = 0;
    out[2][0] = state[2];    out[2][1] = 0; out[2][2] = state[0]-c;
}


FHN::FHN() : System(2) {
    a=0;
    e=0.01;
}

void FHN::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = state[0]-state[0]*state[0]*state[0]/3.0-state[1];
    dvdt[1] = e*(state[0]-a);
}


void FHN::jac(const_it_t &state,
              std::vector<state_t> & out, double time) {
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
    fprintf(gp,"e");
    fclose(gp);
}

void TimeSeries::plotRows(std::vector<int> &idx) {
    std::vector<std::string> files;
    std::string gp_cmd="plot ";
    for (auto id : idx) {
        auto r = row(id);
        std::ostringstream oss;
        oss << id;
        std::string filename = ".nldyn.ts.row"+oss.str() + ".dat";
        files.push_back(filename);
        FILE *f = fopen(filename.c_str(),"w");
        for (auto v : r)
            fprintf(f,"%g\n",v);
        fclose(f);
    }
    for (int i=0; i<files.size() -1 ; i++)
        gp_cmd += "'"+files[i]+"' w l, ";
    gp_cmd += "'"+files.back()+"' w l\n";

    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp,"%s",gp_cmd.c_str());
    fprintf(gp,"e");
    fclose(gp);

    for (auto f : files) remove(f.c_str());
}

void TimeSeries::plotRows() {
    std::vector<int> idx;
    for (int i=0; i<ts[0].size(); i++)
        idx.push_back(i);
    plotRows(idx);
}
TimeSeries::~TimeSeries() {
}
