#include <iostream>
#include "le_solver.h"


#include "systems/ap.hpp"
#include "systems/lorenz.hpp"
#include "systems/fhn3.hpp"

int main(int argc, char **argv) {

    // System *s = new AP3();
    // System *s = new Lorenz();
    System *s = new FHN3();
    s->solve(3000,0.01,s->getState(),0.1);
    auto ts = s->getTs();
    ts.setNoLegend(true);
    std::vector<int> rows_to_plot;
    for (int i=0; i<s->getDim()/2; i++)
        rows_to_plot.push_back(2*i+1);
    ts.plotRows(rows_to_plot,0);

    LyapunovExpsSolver les(s);
    // les.setDbg(true);
    std::cout << les.calcLE(1000,0.01,10000, 1,0.01,s->getState()) << std::endl;
    std::cout << "KYdim: " << les.getKYdim() << std::endl;
    return 0;
}
