#include "fhn3.hpp"
#include <assert.h>

FHN3::FHN3() : System(6) {
    a[0]=0;
    a[1] = 0.5;
    a[2] = -1.1;
    e[0]=e[1]=e[2]=0.01;
    Doo = 0.0;
    Doe = 0.5;
    Deo = 0;
    for (auto &v : vars) v = 0.1;
}

void FHN3::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = state[0]-state[0]*state[0]*state[0]/3.0-state[1]
        +Doo*(state[2]-state[0])+Deo*(state[4]-state[0]);
    dvdt[1] = e[0]*(state[0]-a[0]);

    dvdt[2] = state[2]-state[2]*state[2]*state[2]/3.0-state[3]
        +Doo*(state[0]-state[2])+Deo*(state[4]-state[0]);
    dvdt[3] = e[1]*(state[2]-a[1]);

    dvdt[4] = state[4]-state[4]*state[4]*state[4]/3.0-state[5]
        +Doe*(state[0]+state[2] - 2.0*state[4]);
    dvdt[5] = e[2]*(state[4]-a[2]);

}


void FHN3::jac(const_it_t &state,
              std::vector<state_t> & out, double time) {
    assert(out.size() == dim);

}
