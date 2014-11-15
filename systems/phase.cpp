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

    const int circle_points  = 100;
    const double angle_step = 2*PI/(circle_points - 1);
    const double R = 1;
    std::ofstream circle_ofs("._circle.dat");
    for (int i=0; i<circle_points; i++) {
        float angle = i*angle_step;
        double x = R*cos(angle);
        double y = R*sin(angle);
        circle_ofs << x << " " << y << '\n';
    }
    circle_ofs.close();
#if 0
    std::ofstream data_ofs("._phases.dat");
    data_ofs << R*cos(vars[0]) << " " << R*sin(vars[0]) << "\n\n\n";
    for (int i=1; i<N; i++)
        data_ofs << R*cos(vars[i]) << " " << R*sin(vars[i]) << '\n';
    data_ofs.close();
#endif
    std::ofstream data_ofs("._phases.dat");
    data_ofs << fmod(vars[0],2*PI) << " " << R << "\n\n\n";
    for (int i=1; i<N; i++)
        data_ofs << fmod(vars[i],2*PI) << " " << R << '\n';
    data_ofs.close();

    FILE *gp = popen("gnuplot --persist ", "w");
    // fprintf(gp,"plot '._circle.dat' w l, '._phases.dat' index 0 u 1:2 w p pt 7 ps 2 lc rgb 'green', '._phases.dat' index 1 u 1:2 w p pt 7 ps 2 lc rgb 'blue'\n");
    fprintf(gp,"load 'plotPhases.gpt'\n");
    fclose(gp);
}
