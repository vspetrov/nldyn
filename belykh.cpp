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
#define ALPHA PI/12

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
        "belykh_{N:"+to_str(N)+"}_{om0:"+to_str(omega0)+"}_{delta:"+to_str(delta)+
        "}_{alpha:"+alpha_str+"}.dat"
        );



    std::cout << "*********************************************************\n";
    std::cout << "Experiment start:" << "\n\t*** N=" << N
              << " delta=" << delta << " alpha=" << XSTR(ALPHA) << std::endl;
    if (!singleRun)
        std::cout << "Result: " << filename << std::endl;
    std::cout << "*********************************************************\n";

    std::ofstream ofs;
    if (!singleRun)
        ofs.open(filename);

    double K = 0;
    if (singleRun)
        K=atof(argv[1]);

    int numSteps = 60;
    double Kmax = 1;
    double Kstep = Kmax/(numSteps - 1);
    for (int step = 0; singleRun || step < numSteps; step++, K+=Kstep) {
        Kuramoto Kmt(N);
        Kmt.setK(K);
        Kmt.initOmegasDelta(omega0,delta);
        std::cout << "Phase I: Pre-calc..." << std::endl;
        double t1 = cv::getTickCount();
        Kmt.solve(3000,0.01,Kmt.getState());
        std::cout << "\t->Done. Elapsed: " << ((double)cv::getTickCount() - t1) / cv::getTickFrequency() << std::endl;

        double A = 0;
        Kmt.addAnalyzer(new TimeSeries(0.1, true));

        t1 = cv::getTickCount();
        std::cout << "Phase II: solving with TimeSeries..." << std::endl;
        Kmt.solve(7000,0.01,Kmt.getState());
        std::cout << "\t->Done. Elapsed: " << ((double)cv::getTickCount() - t1) / cv::getTickFrequency() << std::endl;
        t1 = cv::getTickCount();

        std::cout << "Phase III: TimeSeries postprocessing..." << std::endl;

        auto ts  = Kmt.getAnalyzer<TimeSeries>().back()->getAll();
        std::vector<bool> isSync(N);
        std::fill(isSync.begin(), isSync.end(), true);


        for (int i=0; i<ts.size(); i++) {
            auto row = ts[i];
            for (int j=0; j<N; j++) {
                double real_delta = fabs(row[j]-row[0]);
                if ((fmod(real_delta,2*PI) > alpha) && isSync[j]) {
                    isSync[j] = false;
                }
                if (isSync[j]) {
                    double angle = fmod(real_delta,2*PI);
                }
            }
        }

        int numSync = 1;
        for (int i=0; i<N; i++) {
            if (i == 0) continue;
            if (isSync[i]) {
                numSync++;
            }
        }


        auto om_ptr = Kmt.getOmegas();
        int numSyncTheor = 1;
        for (int i=1; i<N; i++) {
            double _delta = om_ptr->at(i)-om_ptr->at(0);
            if (K > fabs(_delta)/(cos(3.0*alpha/2.0)*sin(alpha)))
                numSyncTheor++;
        }

        std::cout << "\t->Done. Elapsed: " << ((double)cv::getTickCount() - t1) / cv::getTickFrequency() << std::endl;
        std::cout << "-----\nStep Completed: K=" << K << " NsNum="
                  << numSync << " NsTh=" << numSyncTheor << std::endl;
        std::cout << "---------------------------------------------------------\n";

        if (singleRun) break;
        ofs << K << " " << numSync << " " << numSyncTheor << "\n";


        // Kmt.plotCircle();
    }
    if (!singleRun)
        ofs.close();

    return 0;
}
