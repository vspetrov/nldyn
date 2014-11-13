#ifndef _PHASE_OSC_SYSTEM_H
#define _PHASE_OSC_SYSTEM_H
#include "nld_sys.h"
const double PI = acos(-1);
class Kuramoto : public System {
private:
    std::vector<double> omega;
    double K;
    int N;
public:
    void setK(double _K) { K = _K; }
    Kuramoto(int _size);
    virtual void rhs(const state_t &state, state_t & out, double time);
    virtual void jac(const state_t &state,
                     matrix_t & out, double time, state_t &dfdt);
    void plotCircle();
    double getAlpha(state_t &state);
    double getAlpha();
    std::vector<double>* getOmegas() { return &omega; }
    void initOmegasDelta(double omega0, double _delta);
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

};


#endif
