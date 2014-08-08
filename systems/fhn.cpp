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
              std::vector<state_t> & out, double time) {
    assert(out.size() == dim);

    out[0][0] = 1-state[0]*state[0];      out[0][1] = -1;
    out[1][0] = e; out[1][1] = 0;
}
