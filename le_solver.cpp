#include "le_solver.h"
#include <algorithm>

extern "C" void dgeqrf_(const int *M, const int *N, double *A, const int *LDA,
                       double *TAU, double *WORK, const int *LWORK,
                       int *INFO);

extern "C" void dorgqr_(const int *M, const int *N, const int *K,
                        double *A, const int *LDA,
                       double *TAU, double *WORK, const int *LWORK,
                       int *INFO);

static inline void Ortonormalize(std::vector<double> &vecs, const int *dim,
                                 double *tau, double *work,
                                 const int *lwork, int *info) {
    dgeqrf_(dim,dim,&*vecs.begin(),dim,tau,work,lwork,info);
    dorgqr_(dim,dim,dim,&*vecs.begin(),dim,tau,work,lwork,info);
}
std::vector<double> LyapunovExpsSolver::calcLE(double warmUpTime,
                                               double wudt,
                                               int numSteps,
                                               double stepTime,
                                               double dt,
                                               state_t &ini) {
    nld_sys->solve(warmUpTime, wudt, ini, dt);
    state_t s0 = nld_sys->getState();
    int dim = nld_sys->getDim();
    int lwork = 2*dim;
    double *work = (double*)malloc(lwork*sizeof(*work));
    double *tau = (double*)malloc(dim*sizeof(double));
    int info;

    std::vector<double> LEs(dim);
    for (int i=0; i<LEs.size(); i++) {
        LEs[i] = 0.0;
    }

    state_t s_ini(LEs);

    for (int i=0; i<numSteps; i++) {
        // std::cout << "Step " << i << " out of " << numSteps << std::endl;
        nld_sys->solve(stepTime,dt,s0,dt);
        s0 = nld_sys->getState();
        auto ts = nld_sys->getTs().getAll();
        std::vector<state_t> LE_matrix;

        for (int j=0; j<dim; j++) {
            state_t s(s_ini);
            s[j] = 1.0;
            nld_sys->solve_linearized(ts, s);
            state_t LE_row = nld_sys->getState();
            LE_matrix.push_back(LE_row);
            // std::cout << "row " << j << ": " << LE_matrix.back()  << ":" <<
                // LE_matrix.back()/L2_norm(LE_matrix.back()) << std::endl;
        }
        std::vector<double> LE_matrix_plain(dim*dim);
        for (int j=0; j<dim; j++) {
            std::copy(LE_matrix[j].begin(),LE_matrix[j].end(),LE_matrix_plain.begin()+dim*j);
        }

        Ortonormalize(LE_matrix_plain, &dim, tau, work, &lwork, &info);
        auto LE_matrix_ortonormed = array2vecvec(LE_matrix_plain, dim);

        for (int j=0; j<dim; j++) {
            // std::cout << "normed " << j << ":" << LE_matrix_ortonormed[j] << std::endl;
            LEs[j] += log(fabs(dot(LE_matrix[j],LE_matrix_ortonormed[j])));
        }
    }
    for (int i=0; i<dim; i++) {
        LEs[i] = LEs[i]/(numSteps*stepTime);
    }
    free(work);
    free(tau);
    return LEs;
}
