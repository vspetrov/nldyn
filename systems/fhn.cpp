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
    num = 3;
    a = new double[num];
    e = new double[num];
    a[0]=0.1;
    a[1] = 0.2;
    a[2] = 0.35;
    e[0]=e[1]=e[2]=0.01;
    e[2] = 0.01; //0.01 - no positive, 0.1 - positive LE
    Doo = 0.001; // 0.001 - positive LE, 0.1 - even large
    Doe = 0.05;
    Deo = 0.00;
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


// ------------------------------------------------------------------------------

FHN_O2_E_chain::FHN_O2_E_chain(int _N) : System((2+_N)*2) {
    char *seedv = getenv("NLD_SEED");
    int seed = 111;
    if (seedv) seed = atoi(seedv);

    srand(seed); //14
    num = _N+2 ;
    a = new double[num];
    e = new double[num];
    a[0]=0.2;
    a[1] = 0.3;
    e[0]=e[1]=0.01;
    for (int i=2; i<num; i++){
        // a[i] = 0.35 + rand()/(float)RAND_MAX*0.15;
        // a[i] = 1 - 0.3*(i-2)/(num-2-1) - 0.2*rand()/(float)RAND_MAX;
        a[i] = 0.5;
        e[i] = 0.01;
        // e[i] = 0.005 + rand()/(float)RAND_MAX*0.015;
    }
    a[num/2] = 0.7;
    Doo = 0; // 0.001 - positive LE, 0.1 - even large
    Doe = 0.1;
    Deo = 0.0;
    Dee = 0.075;
    for (auto &v : vars) v = 0.1;
}

void FHN_O2_E_chain::rhs(const state_t &state, state_t & dvdt, double time) {

    dvdt[0] = state[0]-state[0]*state[0]*state[0]/3.0-state[1]
        +Doo*(state[2]-state[0])+Deo*(state[4]-state[0]) ;
        // +Deo*(state[num-1]-state[0]) ;

    dvdt[1] = e[0]*(state[0]-state[1]+a[0]);

    dvdt[2] = state[2]-state[2]*state[2]*state[2]/3.0-state[3]
        +Doo*(state[0]-state[2])+Deo*(state[4]-state[2]) ;
        // +Deo*(state[num-1]-state[2]);

    dvdt[3] = e[1]*(state[2]-state[3]+a[1]);

    dvdt[4] = state[4]-state[4]*state[4]*state[4]/3.0-state[5]
        +Doe*(state[0]+state[2] - 2.0*state[4]) + Dee*(state[6]-state[4]);
    dvdt[5] = e[2]*(state[4]-state[5]+a[2]);

    for (int i=3; i<num; i++) {
        dvdt[2*i] = state[2*i]-state[2*i]*state[2*i]*state[2*i]/3.0-state[2*i+1]
            +Dee*(state[2*(i-1)]+state[2*((i < num - 1) ? i+1 : i)] - 2*state[2*i]);
        dvdt[2*i+1] = e[i]*(state[2*i]-state[2*i+1]+a[i]);
    }

}

void FHN_O2_E_chain::jac(const state_t &state,
               matrix_t & out, double time, state_t &dfdt) {
    out.clear();
    assert(0);
    std::fill(dfdt.begin(),dfdt.end(),0.0);

}
