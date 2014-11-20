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

#define ALPHA PI/8

template<typename T>
std::string to_str(T value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

static void *threadFunc(void *arg);
int timeout = 1000000;

int main(int argc, char **argv) {
    int N = 100;
    double alpha = ALPHA;


    double mu = 0.1;

    std::string alpha_str(XSTR(ALPHA));
    std::replace(alpha_str.begin(), alpha_str.end(), '/', '_');

    // std::string filename=(
    // "belykh_{N:"+to_str(N)+"}_{om0:"+to_str(omega0)+
    // "}_{delta:"+to_str(delta)+"}_alpha_changing.dat"
    // );



    std::cout << "*********************************************************\n";
    std::cout << "Experiment start:" << "\n\t*** N=" << N
              << " Alpha=" << alpha_str
              << " mu=" << mu << std::endl;

    double K = 0;
    int Na = int(N*mu);
    int Ns = N-Na;
    bool sync = false;
    Kuramoto Kmt(Ns);
    Kmt.initOmegasDeltaRandom(1,1);
    double delta_max = 0;
    auto om = Kmt.getOmegas();
    for (int i=1; i<N; i++) {
        if (fabs((*om)[0]-(*om)[i]) > delta_max) delta_max = fabs((*om)[0]-(*om)[i]);
    }

    double KsyncTheorMu = delta_max/((1-mu)*cos(3*alpha/2)*sin(alpha)-2*mu);
    double KsyncTheor = delta_max/(cos(3*alpha/2)*sin(alpha));
    std::cout << "KsyncTheor=" << KsyncTheor << " KsyncTheorMu=" << KsyncTheorMu << std::endl;



#if 0
    Kmt.setK(K);
    Kmt.solve(5000,0.01,Kmt.getState());
    Kmt.addAnalyzer(new TimeSeries(0.1, true));
    Kmt.solve(1000,0.01,Kmt.getState());
    auto ts  = Kmt.getAnalyzer<TimeSeries>().back()->getAll();
    sync = true;
    for (int i=0; i<ts.size(); i++) {
        auto row = ts[i];
        for (int j=0; j<Kmt.getN(); j++) {
            double real_delta = fabs(row[j]-row[0]);
            if (fmod(real_delta,2*PI) > alpha) {
                sync = false;
                i = ts.size();
                break;
            }
        }
    }
    std::cout << "is sync: " << sync << std::endl;
#endif
    // K = KsyncTheorMu;
    K = 2;

    Kuramoto KmtE(N);
    KmtE.initOmegasExtend(*om, 8, 0.4); // 8, 0.4 K=4 -  cluster
                                      // 8, 2   K=4 - chimera

    {
        om = KmtE.getOmegas();
        std::ofstream ofs("omegas.dat");
        for (int i=0; i<N; i++)
            ofs << (*om)[i] << "\n";
        ofs.close();
    }

    KmtE.setK(K);
    KmtE.solve(5000,0.01,KmtE.getState());
    KmtE.plotCircle();
    sleep(5);
    KmtE.addAnalyzer(new TimeSeries(0.1, true));
    std::vector<double> s1(N),s2(N);


    std::copy(KmtE.getState().begin(),KmtE.getState().end(), s1.begin());

    const double TestTime = 1000;
    KmtE.solve(TestTime,0.01,KmtE.getState());

    std::copy(KmtE.getState().begin(),KmtE.getState().end(), s2.begin());

    std::vector<double> freqs(N);

    for (int i=0; i<N; i++) {
        freqs[i] = (s2[i]-s1[i])/TestTime;
    }

    // std::sort(freqs.begin(),freqs.end());

    std::ofstream ofs("freqs.dat");
    for (int i=0; i<N; i++)
        ofs << freqs[i] << "\n";
    ofs.close();

    {
        std::ofstream ofs("deltas.dat");
        auto ts  = KmtE.getAnalyzer<TimeSeries>().back()->getAll();
        std::vector<bool> sync_all(N);
        std::fill(sync_all.begin(), sync_all.end(), true);

        double row0_prev = fmod(ts[0][0],2*PI);
        for (int i=1; i<ts.size(); i++) {
            auto row = ts[i];

            for (int j=0; j<KmtE.getN(); j++) {
                double real_delta = fabs(row[j]-row[0]);
                if (row0_prev > PI && fmod(row[0],2*PI)<PI)
                    ofs << j << " " << fmod(real_delta,2*PI) << "\n";

                if (fmod(real_delta,2*PI) > alpha) {
                    sync_all[j] = false;
                }
            }
            row0_prev = fmod(row[0],2*PI);
        }
        ofs.close();
        int num_sync = 0;
        for (int i=0; i<N; i++){
            if (sync_all[i]) num_sync++;
        }
        std::cout << "Num sync: " << num_sync << " out of " << N << std::endl;

        pthread_t thread;
        int result=pthread_create(&thread, NULL, threadFunc,NULL);
        for (int i=0; i<ts.size(); i++) {
            auto row = ts[i];
            if (i/1*1 == i) {
                KmtE.plotCircle(row,Ns);
                usleep(timeout);
            }
        }

    }


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
