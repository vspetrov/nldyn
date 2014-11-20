#include "phase.hpp"
#include <assert.h>
#include <fstream>
#include <stdio.h>


Kuramoto::Kuramoto(int _size) : System(_size) {
    N = _size;
    omega.resize(N);
    double om0 = 1;
    double om_delta = 0.01;
    srand(11);
    gp = NULL;
    for (int i=0; i<N; i++) {
        omega[i] = om0 + om_delta*i;
        vars[i] = i*2*PI/N +fmod(rand(), PI/10);

    }

}

void Kuramoto::initOmegasDelta(double omega0, double delta_max) {

    for (int i=0; i<N; i++) {
        omega[i] = omega0 + i*delta_max/(N-1);
    }
}

void Kuramoto::initOmegasDeltaRandom(double omega0, double delta_max) {

    omega[0] = omega0;
    for (int i=1; i<N; i++) {
        omega[i] = omega0 + rand()/(double)RAND_MAX*delta_max;
    }
}

void Kuramoto::initOmegasInterval(double omin, double omax) {
    for (int i=0; i<N; i++) {
        omega[i] = omin + rand()/(double)RAND_MAX*(omax - omin);
    }
}

void Kuramoto::initOmegasExtend(const std::vector<double> &om_input, double omin, double delta){
    std::copy(om_input.begin(),om_input.end(),omega.begin());
    for (int i=om_input.size(); i<N; i++)
        omega[i] = omin + rand()/(double)RAND_MAX*delta;
}
void Kuramoto::rhs(const state_t &state, state_t & dvdt, double time) {
    double r = 0, psi = 0;
    double sin_sum = 0, cos_sum = 0;
#pragma omp parallel for default(shared) reduction(+:sin_sum,cos_sum)
    for (int i=0; i<N; i++) {
        sin_sum += sin(state[i]);
        cos_sum += cos(state[i]);
    }
    r = sqrt(sin_sum*sin_sum + cos_sum*cos_sum)/N;
    psi = atan2(sin_sum,cos_sum);

#pragma omp parallel for
    for (int i=0; i<N; i++)
        dvdt[i] = omega[i] + K*r*sin(psi - state[i]);
}


void Kuramoto::jac(const state_t &state,
             matrix_t & out, double time, state_t &dfdt) {
    assert(0);
}


double Kuramoto::getAlpha() {
    return Kuramoto::getAlpha<state_t>(vars);
}

void Kuramoto::plotCircle() {
    plotCircle(vars);
}
