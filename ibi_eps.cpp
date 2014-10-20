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

#include <opencv/cv.h>
int main(int argc, char **argv) {
    System *s;
    std::string rst_dir_prefix = ".";
    if (argc > 1)
        rst_dir_prefix = std::string(argv[1]);

    std::ofstream ofs(rst_dir_prefix + "/ibi_eps.dat");
    std::ofstream ofs_sync(rst_dir_prefix + "/ibi_eps_sync.dat");

    double e_min = 0.005;
    double e_max = 0.02;
    int steps = 100;
    for (double epsilon = e_min; epsilon < e_max+1E-6; epsilon += (e_max - e_min)/(steps-1)){ //0.00005
        FHN3 fhn3;
        fhn3.setDeo(0.01);
        fhn3.setDoe(0.1);
        fhn3.setDoo(0);
        fhn3.setEpsilon(2,epsilon);
        s = &fhn3;
        std::cout << "epsilon = " << epsilon << std::endl;
        std::cout << "Transient skip maps" << std::endl;
        s->solve(20000,0.01,s->getState());
        typedef ConstSectionMapper MT;
        s->addAnalyzer(new MT(4,1.0));
        std::cout << "Calculating IBIs" << std::endl;
        s->solve(500000,0.01,s->getState());
        auto map = s->getAnalyzer<MT>().back()->getMap(4);
        for (int i=0; i<map->size(); i++) {
            double v = (*map)[i].first;
            ofs << epsilon << " " << v << std::endl;
        }
        s->clearAnalyzers();
        fhn3.setDoo(0.1);
        s->solve(20000,0.01,s->getState());
        s->addAnalyzer(new MT(4,1.0));
        std::cout << "Calculating IBIs SYNC" << std::endl;
        s->solve(50000,0.01,s->getState());
        map = s->getAnalyzer<MT>().back()->getMap(4);
        double ibi_av_sync = 0;
        for (int i=0; i<map->size(); i++) {
            ibi_av_sync += (*map)[i].first;
        }
        ibi_av_sync /= map->size();
        ofs_sync << epsilon << " " << ibi_av_sync << std::endl;
    }
    ofs_sync.close();
    ofs.close();
    std::cout << "done" << std::endl;

    return 0;
}
