#include "lorenz.hpp"
#include <assert.h>

Lorenz::Lorenz() : System(3) {
    this->sigma = 10.0;
    this->beta  = 8./3.;
    this->ro    = 28.0;
    for (auto &v : vars) v = 1.0;
}

void Lorenz::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = sigma*(state[1]-state[0]);
    dvdt[1] = state[0]*(ro - state[2]) - state[1];
    dvdt[2] = state[0]*state[1]-beta*state[2];
}


void Lorenz::jac(const_it_t &state,
                 state_t & out, double time) {
    assert(out.size() == dim*dim);

    auto it = out.begin();
    it[0] = -sigma;      it[1] = sigma;    it[2] = 0; it += dim;
    it[0] = ro-state[2]; it[1] = -1;       it[2] = -state[0]; it += dim;
    it[0] = state[1];    it[1] = state[0]; it[2] = -beta;
}
