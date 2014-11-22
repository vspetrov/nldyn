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

#include <pthread.h>
#include <unistd.h>

#define XSTR(X) STR(X)
#define STR(X) #X

#define ALPHA PI/3

template<typename T>
std::string to_str(T value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

static void *threadFunc(void *arg);
int timeout = 1000000;



static inline double findAlpha(double K, double dmax) {
    double a = 0;
    do{
        if (K*cos(3*a/2)*sin(a) > dmax) {
            return a;
        }
    }while ((a < 1 + 1e-5) && (a += 0.0001));
    return -1;
}
int main(int argc, char **argv) {
    int N = 100;
    double alpha = ALPHA;


    std::string alpha_str(XSTR(ALPHA));
    std::replace(alpha_str.begin(), alpha_str.end(), '/', '_');

    // std::string filename=(
    // "belykh_{N:"+to_str(N)+"}_{om0:"+to_str(omega0)+
    // "}_{delta:"+to_str(delta)+"}_alpha_changing.dat"
    // );



    std::cout << "*********************************************************\n";
    std::cout << "Experiment start:" << "\n\t*** N=" << N << std::endl;


    double K = 1;


    alpha = 1.0;

    int steps = 20;
    double Kmax = 8;
    double dK = Kmax/(steps-1);
    // std::ofstream ofs("alpha_K_omega_1_delta_1_N_100.dat");
    // for (int s=0; s<steps; s++) {
    //     K = dK*s;
    //     Kuramoto Kmt(N);
    //     Kmt.initOmegasDelta(1,1);
    //     double delta_max = 0;
    //     auto om = Kmt.getOmegas();
    //     for (int i=1; i<N; i++) {
    //         if (fabs((*om)[0]-(*om)[i]) > delta_max) delta_max = fabs((*om)[0]-(*om)[i]);
    //     }

    //     Kmt.setK(K);

    //     Kmt.solve(5000,0.01,Kmt.getState());
    //     Kmt.addAnalyzer(new TimeSeries(0.1, true));
    //     Kmt.solve(1000,0.01,Kmt.getState());
    //     auto ts  = Kmt.getAnalyzer<TimeSeries>().back()->getAll();
    //     double max_angle = 0;

    //     for (int i=0; i<ts.size(); i++) {
    //         auto row = ts[i];
    //         for (int j=0; j<Kmt.getN(); j++) {
    //             double angle = fabs(row[j]-row[0]);
    //             if (fmod(angle,2*PI) > max_angle) {
    //                 max_angle = fmod(angle,2*PI);
    //             }
    //         }
    //     }
    //     std::cout << "max angle=" << max_angle
    //               << " alpha=" << alpha
    //               << " K=" << K
    //               << " AlphaTheor=" << findAlpha(K,delta_max) << std::endl;
    //     if (max_angle < 6.2) {
    //         ofs << K << " " << max_angle;
    //         double A = findAlpha(K,delta_max);
    //         if (A > 0) {
    //             ofs << " " << A;
    //         }
    //         ofs << "\n";
    //     }

    //     ofs.close();
    // }

    Kuramoto Kmt(N);
    Kmt.initOmegasDelta(1,1);
    Kmt.setK(8);
    Kmt.solve(1005,0.01,Kmt.getState());
    Kmt.plotCircle();

    return 0;
}

static void *threadFunc(void *arg) {
    std::string str;
    while (1) {
        std::cin >> str;
        std::cout << "GOT " << str << std::endl;
        if (str == "u")
            timeout /= 2;

        if (str == "d")
            timeout *= 2;
    }
}
