#include "ap.hpp"
#include <assert.h>

AP::AP() : System(2) {
    a=0.05;
    k=8;
    ksi = -0.5;
    vars[0] = 0.1;
    vars[1] = 0;
}

const static double _eps_min = 0.1;
const static double _eps_max = 1.0;
const static double _eps_rate = 3.0;
const static double _eps_th   = 0.05;

static inline double _eps(double v) {
    return _eps_min+(_eps_max-_eps_min)/(1.0+exp(-(_eps_th-v)*_eps_rate));
}

static inline double _eps_d(double v) {
    double E = exp(-(_eps_th-v)*_eps_rate);
    -(_eps_max-_eps_min)/((1.0+E)*(1.0+E))*E*_eps_rate;
}

void AP::rhs(const state_t &state, state_t & dvdt, double time) {
    dvdt[0] = k*state[0]*(1.0-state[0])*(state[0]-a)-state[0]*state[1];
    dvdt[1] = _eps(state[0])*(k*state[0]-state[1] + ksi);
}


void AP::jac(const state_t &state,
             matrix_t & out, double time, state_t &dfdt) {
    out(0,0) = -3.0*k*state[0]*state[0]+2.0*k*(1.0+a)*state[0]-k*a-state[1]; out(0,1) = -state[0];
    double eps = _eps(state[0]);
    out(1,0) = eps*k+(k*state[0]-state[1]+ksi)*_eps_d(state[0]);                        out(1,1) = -eps;
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
}


//------------------------------ AP 3

AP3::AP3() : System(6) {
    a[0]=a[1]=a[2]=0.05;
    k[0]=k[1]=k[2]=8.0;
    ksi[0] = -0.5;
    ksi[1] = -0.7;
    ksi[2] = 0;
    for (int i=0; i<num; i++) {
        vars[i*2] = 0.1;
        vars[i*2+1] = 0;
    }
    Doo = 0.005;
    Doe = 0.05;
    Deo = 0.01;
}

void AP3::rhs(const state_t &state, state_t & dvdt, double time) {
    double eps[3];
    for (int i=0; i< num; i++)
        eps[i] = (state[i*2] < 0.05) ? 1.0 : 0.1;

    dvdt[0] = k[0]*state[0]*(1.0-state[0])*(state[0]-a[0])-state[0]*state[1] +
        Doo*(state[2]-state[0]) + Deo*(state[4]-state[0]);
    dvdt[1] = eps[0]*(k[0]*state[0]-state[1] + ksi[0]);

    dvdt[2] = k[1]*state[2]*(1.0-state[2])*(state[2]-a[1])-state[2]*state[3] +
        Doo*(state[0]-state[2]) + Deo*(state[4]-state[2]);
    dvdt[3] = eps[1]*(k[1]*state[2]-state[3] + ksi[1]);

    dvdt[4] = k[2]*state[4]*(1.0-state[4])*(state[4]-a[2])-state[4]*state[5] +
        Doe*(state[0]+state[2]-2.0*state[4]);
    dvdt[5] = eps[2]*(k[2]*state[4]-state[5] + ksi[2]);


}


void AP3::jac(const state_t &state,
              matrix_t & out, double time, state_t &dfdt) {
    double eps[3];
    for (int i=0; i< num; i++)
        eps[i] = (state[i*2] < 0.05) ? 1.0 : 0.1;


    out.clear();
    std::fill(dfdt.begin(),dfdt.end(),0.0);

    out(0,0) = -3.0*k[0]*state[0]*state[0]+2.0*k[0]*(1.0+a[0])*state[0]-k[0]*a[0]-state[1]-Doo-Deo;
    out(0,1) = -state[0]; out(0,2) = Doo; out(0,4) = Deo;
    out(1,0) = eps[0]*k[0];                        out(1,1) = -eps[0];

    out(2,0) = Doo;
    out(2,2) = -3.0*k[1]*state[2]*state[2]+2.0*k[1]*(1.0+a[1])*state[2]-k[1]*a[1]-state[3]-Doo-Deo;
    out(2,3) = -state[2];  out(2,4) = Deo;
    out(3,2) = eps[1]*k[1];                        out(3,3) = -eps[1];

    out(4,0) = Deo; out(4,2) = Deo;
    out(4,4) = -3.0*k[2]*state[4]*state[4]+2.0*k[2]*(1.0+a[2])*state[4]-k[2]*a[2]-state[5]-2*Doe;
    out(4,5) = -state[4];
    out(5,4) = eps[2]*k[2];                        out(5,5) = -eps[2];

}
