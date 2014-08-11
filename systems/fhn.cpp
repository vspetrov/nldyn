#include "fhn.hpp"
#include <assert.h>

FHN::FHN() : System(2) {
    a=0;
    e=0.01;
}

void FHN::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = state[0]-state[0]*state[0]*state[0]/3.0-state[1];
    dvdt[1] = e*(state[0]-a);
}


void FHN::jac(const_it_t &state,
              state_t & out, double time) {
    assert(out.size() == dim*dim);
    auto it = out.begin();

    it[0] = 1-state[0]*state[0];      it[1] = -1; it += dim;
    it[0] = e;                        it[1] = 0;
}
