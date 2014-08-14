#include "fhn.hpp"
#include <assert.h>

FHN::FHN() : System(2) {
    a=0;
    e=0.01;
    vars[0] = 0.1;
    vars[1] = 0.1;
}

void FHN::rhs(const state_t &state, state_t & dvdt, double time) {
    dvdt[0] = state[0]-state[0]*state[0]*state[0]/3.0-state[1];
    dvdt[1] = e*(state[0]-a);
}


void FHN::jac(const state_t &state,
              matrix_t & out, double time, state_t &dfdt) {
    out(0,0) = 1-state[0]*state[0];      out(0,1) = -1;
    out(1,0) = e;                        out(1,1) = 0;

    dfdt[0] = 0.0; dfdt[1] = 0.0;
}


FHN3::FHN3() : System(6) {
    a[0]=0.1;
    a[1] = 0.2;
    a[2] = 0.35;
    e[0]=e[1]=e[2]=0.02;
    Doo = 0.001;
    Doe = 0.005;
    Deo = 0.005;
    for (auto &v : vars) v = 0.1;
}

void FHN3::rhs(const state_t &state, state_t & dvdt, double time) {
    dvdt[0] = state[0]-state[0]*state[0]*state[0]/3.0-state[1]
        +Doo*(state[2]-state[0])+Deo*(state[4]-state[0]);
    dvdt[1] = e[0]*(state[0]-state[1]+a[0]);

    dvdt[2] = state[2]-state[2]*state[2]*state[2]/3.0-state[3]
        +Doo*(state[0]-state[2])+Deo*(state[4]-state[0]);
    dvdt[3] = e[1]*(state[2]-state[3]+a[1]);

    dvdt[4] = state[4]-state[4]*state[4]*state[4]/3.0-state[5]
        +Doe*(state[0]+state[2] - 2.0*state[4]);
    dvdt[5] = e[2]*(state[4]-state[5]+a[2]);

}


void FHN3::jac(const state_t &state,
               matrix_t & out, double time, state_t &dfdt) {
    out.clear();

    out(0,0) = 1-state[0]*state[0]-Doo-Deo; out(0,1) = -1; out(0,2) = Doo; out(0,4) = Deo;
    out(1,0) = e[0]; out(1,1) = -e[0];

    out(2,0) = Doo; out(2,2) = 1-state[2]*state[2]-Doo-Deo; out(2,3) = -1; out(2,4) = Deo;
    out(3,2) = e[1]; out(3,3) = -e[1];

    out(4,0) = Doe; out(4,2) = Doe; out(4,4) = 1-state[4]*state[4]-2.0*Doe; out(4,5) = -1;
    out(5,4) = e[2]; out(5,5) = -e[2];

    std::fill(dfdt.begin(),dfdt.end(),0.0);

}
