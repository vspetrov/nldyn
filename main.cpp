#include <iostream>
#include "nld_sys.h"
#include "le_solver.h"

int main(int argc, char **argv) {
    // Lorenz lrnz;
    // std::vector<double> ini = {0.1,0.1,0.1};
    // lrnz.solve(100,0.001, ini, 0.01);
    // auto ts = lrnz.getTs();
    // ts.plot(2,3);
    // LyapunovExpsSolver les(&lrnz);
    // std::cout << les.calcLE(1000,0.01,10000,0.1,0.001,ini) << std::endl;

    // std::vector<double> ini = {1,2,3};
    // FHN fhn;
    // // fhn.solve(1000,0.01,ini,0.01);
    // // auto ts = fhn.getTs();
    // // ts.plot(0,1);
    // LyapunovExpsSolver les(&fhn);
    // std::cout << les.calcLE(100,0.01,10000,0.01,0.001,ini) << std::endl;

    std::vector<double> ini = {1,2,3};
    Rossler rslr;
    // // rslr.solve(1000,0.01,ini,0.01);
    // // auto ts = rslr.getTs();
    // // ts.plot(2,1);
    LyapunovExpsSolver les(&rslr);
    std::cout << les.calcLE(500,0.01,50,1,0.01,ini,0.01) << std::endl;
    std::cout << "KYdim: " << les.getKYdim() << std::endl;

    return 0;
}
