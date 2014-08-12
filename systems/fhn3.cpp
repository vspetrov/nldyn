#include "fhn3.hpp"
#include <assert.h>

FHN3::FHN3() : System(6) {
    a[0]=0.1;
    a[1] = 0.2;
    a[2] = 0.35;
    e[0]=e[1]=e[2]=0.2;
    Doo = 0.002;
    Doe = 0.05;
    Deo = 0.02;
    for (auto &v : vars) v = 0.1;
}

void FHN3::rhs(const_it_t &state, it_t & dvdt, double time) {
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


void FHN3::jac(const_it_t &state,
              state_t & out, double time) {
    assert(out.size() == dim*dim);
    std::fill(out.begin(),out.end(),0.0);

    auto it = out.begin();
    it[0] = 1-state[0]*state[0]-Doo-Deo; it[1] = -1; it[2] = Doo; it[4] = Deo; it += dim;
    it[0] = e[0]; it[1] = -e[0]; it += dim;

    it[0] = Doo; it[2] = 1-state[2]*state[2]-Doo-Deo; it[3] = -1; it[4] = Deo; it += dim;
    it[2] = e[1]; it[3] = -e[1]; it += dim;

    it[0] = Doe; it[2] = Doe; it[4] = 1-state[4]*state[4]-2.0*Doe; it[5] = -1; it += dim;
    it[4] = e[2]; it[5] = -e[2];

}
