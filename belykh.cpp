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

int main(int argc, char **argv) {
    int N = 10;
    Kuramoto Kmt(N);

    double delta = 0.02;
    double K = atof(argv[1]);
    double alpha = PI/6;
    Kmt.setK(K);
    Kmt.initOmegasDelta(2,delta);
    std::cout << "Phase I: Pre-calc..." << std::endl;
    double t1 = cv::getTickCount();
    Kmt.solve(5000,0.01,Kmt.getState());
    std::cout << "Phase I done. Elapsed: " << ((double)cv::getTickCount() - t1) / cv::getTickFrequency() << std::endl;

    double A = 0;
    Kmt.addAnalyzer(new TimeSeries(0.1, true));

    t1 = cv::getTickCount();
    std::cout << "Phase II: solving with TimeSeries..." << std::endl;
    Kmt.solve(5000,0.01,Kmt.getState());
    std::cout << "Phase II done. Elapsed: " << ((double)cv::getTickCount() - t1) / cv::getTickFrequency() << std::endl;



    t1 = cv::getTickCount();
    std::cout << "Phase III: TimeSeries postprocessing..." << std::endl;

    auto ts  = Kmt.getAnalyzer<TimeSeries>().back()->getAll();
    std::vector<bool> isSync(N);
    std::fill(isSync.begin(), isSync.end(), true);
    std::vector<double> maxAngle(N);
    std::fill(maxAngle.begin(), maxAngle.end(), 0.0);


    for (int i=0; i<ts.size(); i++) {
        auto row = ts[i];
        for (int j=0; j<N; j++) {
            double real_delta = fabs(row[j]-row[0]);
            if ((fmod(real_delta,2*PI) > alpha) && isSync[j]) {
                isSync[j] = false;
            }
            if (isSync[j]) {
                double angle = fmod(real_delta,2*PI);
                if (angle > PI) angle = -(2*PI-angle);
                if (fabs(angle) > fabs(maxAngle[j])) maxAngle[j] = angle;
            }
        }
    }
    std::cout << "Phase III done. Elapsed: " << ((double)cv::getTickCount() - t1) / cv::getTickFrequency() << std::endl;

    int numSync = 1;
    double Alpha;
    double maxPostive = 0, maxNegative = 0, maxFabs = 0;
    for (int i=0; i<N; i++) {
        if (i == 0) continue;
        if (isSync[i]) {
            numSync++;
            if (maxAngle[i] > 0 && maxAngle[i] > maxPostive)
                maxPostive = maxAngle[i];
            if (maxAngle[i] < 0 && maxAngle[i] < maxNegative)
                maxNegative = maxAngle[i];
            if (fabs(maxAngle[i]) > maxFabs)
                maxFabs = fabs(maxAngle[i]);
        }
    }
    Alpha = maxPostive - maxNegative;
    std::cout << "NumSynced: " << numSync << " out of " << N << std::endl;
    std::cout << "Alpha: " << Alpha << " rad; " << Alpha/PI*180.0 << " degrees" << std::endl;
    std::cout << "maxIndAlpha: " << maxFabs << " rad; " << maxFabs/PI*180.0 << " degrees" << std::endl;
    auto om_ptr = Kmt.getOmegas();
    int numSyncTheor = 1;
    for (int i=1; i<N; i++) {
        double _delta = om_ptr->at(i)-om_ptr->at(0);
        if (K > fabs(_delta)/(cos(3.0*alpha/2.0)*sin(alpha)))
            numSyncTheor++;
    }
    std::cout << "NumSynced Theor: " << numSyncTheor << std::endl;
#if 0
    std::ofstream ofs("ts.dat");
    for (auto &v : ts) {
        ofs << v[0] << " " << v[5] << " " << v[9] << std::endl;
    }
    ofs.close();
#endif
    Kmt.plotCircle();
    return 0;
}
