#include "rossler.hpp"
#include <assert.h>


Rossler::Rossler() : System(3) {
    a = 0.2;
    b = 0.2;
    c = 5.7;
}

void Rossler::rhs(const_it_t &state, it_t & dvdt, double time) {
    dvdt[0] = -state[1]-state[2];
    dvdt[1] = state[0]+a*state[1];
    dvdt[2] = b+state[2]*(state[0]-c);
}


void Rossler::jac(const_it_t &state,
                  state_t & out, double time) {
    assert(out.size() == dim*dim);
    auto it = out.begin();
    it[0] = 0;        it[1] = -1;  it[2] = -1; it += dim;
    it[0] = 1;        it[1] = a;   it[2] = 0; it += dim;
    it[0] = state[2]; it[1] = 0;   it[2] = state[0]-c;
}
