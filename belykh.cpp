#include <iostream>
#include "systems/phase.hpp"
#include "le_solver.h"
#include "analyzers/time_series.hpp"
#include "analyzers/mapper.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <opencv/cv.h>
#include <sstream>
#include <string>
#include <exception>

#define XSTR(X) STR(X)
#define STR(X) #X
#define ALPHA PI/6

template<typename T>
std::string to_str(T value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

int main(int argc, char **argv) {
    int N = 10;
    double delta = 0.02;
    double alpha = ALPHA;
    double omega0 = 2;

    bool singleRun = false;
    if (argc > 1) singleRun  = true;

    std::string alpha_str(XSTR(ALPHA));
    std::replace(alpha_str.begin(), alpha_str.end(), '/', '_');

    std::string filename=(
        "belykh_{N:"+to_str(N)+"}_{om0:"+to_str(omega0)+
        "}_{alpha:"+alpha_str+"}.dat"
        );



    std::cout << "*********************************************************\n";
    std::cout << "Experiment start:" << "\n\t*** N=" << N
              << " alpha=" << XSTR(ALPHA) << std::endl;
    if (!singleRun)
        std::cout << "Result: " << filename << std::endl;
    std::cout << "*********************************************************\n";

    std::ofstream ofs;
    if (!singleRun)
        ofs.open(filename);

    double K = 0;
    if (singleRun) {
        K=atof(argv[1]);
        delta = atof(argv[2]);
        Kuramoto Kmt(N);
        std::cout << "K = " << K << " delta=" << delta << std::endl;
        Kmt.initOmegasDelta(omega0,delta);
        Kmt.setK(K);
        Kmt.solve(5000,0.01,Kmt.getState());
        Kmt.plotCircle();
        return 0;
    }


    int num_steps = 20;
    double delta_start = 0.0;
    double delta_max = 2;
    double delta_step = (delta_max - delta_start)/(num_steps - 1);

    double t1 = cv::getTickCount();
    ofs << "0 0 0" << std::endl;
    double Kstep = 0.05;
    for (int dc =  1; dc < num_steps; dc++) {
        delta = delta_start + delta_step*dc;
        bool sync = false;
        Kuramoto Kmt(N);
        Kmt.initOmegasDelta(omega0,delta);
        std::cout << "Step " << dc << " out of " << num_steps
                  << "; delta = " << delta
                  <<"; Elapsed: " << (cv::getTickCount() - t1)/cv::getTickFrequency() << std::endl;
        do {
            Kmt.setK(K);
            Kmt.solve(5000,0.01,Kmt.getState());

            Kmt.addAnalyzer(new TimeSeries(0.1, true));
            Kmt.solve(1000,0.01,Kmt.getState());

            // Kmt.plotCircle();sleep(2);
            auto ts  = Kmt.getAnalyzer<TimeSeries>().back()->getAll();
            sync = true;
            for (int i=0; i<ts.size(); i++) {
                auto row = ts[i];
                for (int j=0; j<N; j++) {
                    double real_delta = fabs(row[j]-row[0]);
                    if (fmod(real_delta,2*PI) > alpha) {
                        sync = false;
                        i = ts.size();
                        break;
                    }
                }
            }
        } while (!sync && (K += Kstep));

        double delta_max = 0;
        auto om = Kmt.getOmegas();
        for (int i=1; i<N; i++) {
            if (fabs((*om)[0]-(*om)[i]) > delta_max) delta_max = fabs((*om)[0]-(*om)[i]);
        }
        std::cout  <<"\t --> SYNC: K=" << K << " Ktheor=" <<  delta_max/(cos(3*alpha/2.0)*sin(alpha)) << std::endl;
        ofs << delta_max << " " << K << " " << delta_max/(cos(3*alpha/2.0)*sin(alpha)) << "\n";
    }
    if (!singleRun)
        ofs.close();


    return 0;
}
