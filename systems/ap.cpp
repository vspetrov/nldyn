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

void AP::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = k*state[0]*(1.0-state[0])*(state[0]-a)-state[0]*state[1];
    dvdt[1] = _eps(state[0])*(k*state[0]-state[1] + ksi);
}


void AP::jac(const_it_t &state,
              state_t & out, double time) {
    assert(out.size() == dim*dim);
    auto it = out.begin();

    it[0] = -3.0*k*state[0]*state[0]+2.0*k*(1.0+a)*state[0]-k*a-state[1]; it[1] = -state[0];
    it += dim;
    double eps = _eps(state[0]);
    it[0] = eps*k+(k*state[0]-state[1]+ksi)*_eps_d(state[0]);                        it[1] = -eps;
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
    Doo = 0.001;
    Doe = 0.02;
    Deo = 0.01;
}

void AP3::rhs(const_it_t &state, it_t & dvdt, double time) {
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


void AP3::jac(const_it_t &state,
              state_t & out, double time) {
    assert(out.size() == dim*dim);
    double eps[3];
    for (int i=0; i< num; i++)
        eps[i] = (state[i*2] < 0.05) ? 1.0 : 0.1;

    auto it = out.begin();

    it[0] = -3.0*k[0]*state[0]*state[0]+2.0*k[0]*(1.0+a[0])*state[0]-k[0]*a[0]-state[1]-Doo-Deo;
    it[1] = -state[0]; it[2] = Doo; it[4] = Deo;
    it += dim;
    it[0] = eps[0]*k[0];                        it[1] = -eps[0];     it += dim;

    it[0] = Doo;
    it[2] = -3.0*k[1]*state[2]*state[2]+2.0*k[1]*(1.0+a[1])*state[2]-k[1]*a[1]-state[3]-Doo-Deo;
    it[3] = -state[2];  it[4] = Deo;     it += dim;
    it[2] = eps[1]*k[1];                        it[3] = -eps[1];     it += dim;

    it[0] = Deo; it[2] = Deo;
    it[4] = -3.0*k[2]*state[4]*state[4]+2.0*k[2]*(1.0+a[2])*state[4]-k[2]*a[2]-state[5]-2*Doe;
    it[5] = -state[4];      it += dim;
    it[4] = eps[2]*k[2];                        it[5] = -eps[2];

}
