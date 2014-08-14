#include <iostream>
#include "systems/fhn.hpp"
#include "systems/lorenz.hpp"
#include "le_solver.h"
int main(int argc, char **argv) {
    System *s;
    FHN3 fhn3;


    s = &fhn3;
    // s->setStepper(System::STEPPER_RK_CK54);
    s->setSaveTS(0.01);
    s->solve(3000,0.01,s->getState());
    auto ts = s->getTs();
    std::vector<int> rows_to_plot;
    rows_to_plot.push_back(1);
    rows_to_plot.push_back(3);
    rows_to_plot.push_back(5);
    ts.setNoLegend(true);
    ts.plotRows(rows_to_plot, 0);

    LyapunovExpsSolver les(s);
    auto LEs = les.calcLE(5000,0.01,50000,1,0.01,s->getState());
    for (auto &v : LEs) std::cout << v << std::endl;

    return 0;
}
