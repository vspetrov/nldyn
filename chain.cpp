#include <iostream>
#include "systems/fhn.hpp"
#include "systems/lorenz.hpp"
#include "systems/ap.hpp"
#include "le_solver.h"
#include "analyzers/time_series.hpp"
#include "analyzers/mapper.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>

int main(int argc, char **argv) {
    FHN_O2_E_chain fhn_chain(50);
    TimeSeries ts(0.5);
    fhn_chain.addAnalyzer(&ts);


    fhn_chain.solve(5000, 0.01, fhn_chain.getState());
    ts.saveBinary("chain_rst.bin", false, 1);
    // std::vector<int> idx(5);
    // idx[0] = 1;
    // idx[1] = 3;
    // idx[2] = 5;
    // idx[3] = 11;
    // idx[4] = 51;
    // ts.plotRows(idx,0);
    return 0;
}
