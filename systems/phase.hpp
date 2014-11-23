#ifndef _PHASE_OSC_SYSTEM_H
#define _PHASE_OSC_SYSTEM_H
#include "nld_sys.h"
const double PI = acos(-1);
class Kuramoto : public System {
private:
    std::vector<double> omega;
    double K;
    int N;
    FILE *gp;
public:
    void setK(double _K) { K = _K; }
    Kuramoto(int _size);
    ~Kuramoto() {if (gp) fclose(gp); }
    virtual void rhs(const state_t &state, state_t & out, double time);
    virtual void jac(const state_t &state,
                     matrix_t & out, double time, state_t &dfdt);
    void plotCircle();
    double getAlpha(state_t &state);
    double getAlpha();
    std::vector<double>* getOmegas() { return &omega; }
    void initOmegasDelta(double omega0, double _delta);
    void initOmegasDeltaRandom(double omega0, double delta);
    void initOmegasInterval(double omin, double omax);
    void initOmegasExtend(const std::vector<double> &om_input, double omin, double delta);
    template<typename T>
    static double getAlpha(const T &state) {
        double phase_min = fmod(state[0],2*PI);
        double phase_max = phase_min;
        for (int i=1; i<state.size(); i++) {
            double v = fmod(state[i],2*PI);
            if (v < phase_min) phase_min = v;
            if (v > phase_max) phase_max = v;
        }
        return phase_max - phase_min;
    }
    int getN() { return N;}

    template<typename T>
    void plotCircle(const T &state, int group = 1) {
        #if 0
        std::ofstream data_ofs("._phases.dat");
        const double R = 1;
        for (int i=0; i<group; i++)
            data_ofs << fmod(state[0],2*PI) << " " << R << '\n';
        data_ofs  << "\n\n";

        for (int i=group; i<N; i++)
            data_ofs << fmod(state[i],2*PI) << " " << R << '\n';
        data_ofs.close();
        if (NULL == gp) {
            gp = popen("gnuplot ", "w");
            // fprintf(gp,"plot '._circle.dat' w l, '._phases.dat' index 0 u 1:2 w p pt 7 ps 2 lc rgb 'green', '._phases.dat' index 1 u 1:2 w p pt 7 ps 2 lc rgb 'blue'\n");
            fprintf(gp,"load 'plotPhases.gpt'\n");
        } else {
            fprintf(gp,"replot\n");
        }
        fflush(gp);
        #endif
    }

};


#endif
