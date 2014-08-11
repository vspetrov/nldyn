#include <iostream>
#include "le_solver.h"


#include "systems/fhn3.hpp"
#include "systems/lorenz.hpp"

int main(int argc, char **argv) {
    // FHN3 fhn3;
    // fhn3.solve(5000,0.01,fhn3.getState(),0.1);
    // auto ts = fhn3.getTs();
    // std::vector<int> rows_to_plot;
    // for (int i=0; i<3; i++)
        // rows_to_plot.push_back(2*i+1);
    // ts.plotRows(rows_to_plot,0);
    Lorenz lrz;
    LyapunovExpsSolver les(&lrz);
    std::cout << les.calcLE(100,0.01,1000,1,0.01,lrz.getState(),0.01) << std::endl;
    std::cout << "KYdim: " << les.getKYdim() << std::endl;

    return 0;
}
