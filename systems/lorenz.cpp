#include "lorenz.hpp"
#include <assert.h>

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
