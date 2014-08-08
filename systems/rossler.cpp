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
                  std::vector<state_t> & out, double time) {
    assert(out.size() == dim);

    out[0][0] = 0;      out[0][1] = -1;    out[0][2] = -1;
    out[1][0] = 1; out[1][1] = a;       out[1][2] = 0;
    out[2][0] = state[2];    out[2][1] = 0; out[2][2] = state[0]-c;
}
