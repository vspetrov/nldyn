#include <iostream>
#include "nld_sys.h"
extern "C" void dgeqrf_(const int *M, const int *N, double *A, const int *LDA,
                       double *TAU, double *WORK, const int *LWORK,
                       int *INFO);

extern "C" void dorgqr_(const int *M, const int *N, const int *K,
                        double *A, const int *LDA,
                       double *TAU, double *WORK, const int *LWORK,
                       int *INFO);


int main(int argc, char **argv) {
    Lorenz lrnz(3);
    std::vector<double> ini = {1,2,3};
    lrnz.solve(10,0.001, ini, 0.01);
    auto ts = lrnz.getTs();
    // ts.plot(1,2);
#if 0
    std::vector<double> jac;
    lrnz.jac(ini, jac, 0);
    auto vecs = array2vecvec(jac,3);
    for (auto v : vecs) std::cout << v << std::endl;
    std::cout << "1x2 = " << dot(vecs[0],vecs[1]) << std::endl;
    std::cout << "1x3 = " << dot(vecs[0],vecs[2]) << std::endl;
    std::cout << "3x2 = " << dot(vecs[2],vecs[1]) << std::endl;
    const int dim = 3;
    const int lwork = 6;
    double tau[dim];
    double work[lwork];
    int info;
    dgeqrf_(&dim,&dim,&*jac.begin(),&dim,tau,work,&lwork,&info);
    assert(info == 0);
    dorgqr_(&dim,&dim,&dim,&*jac.begin(),&dim,tau,work,&lwork,&info);
    assert(info == 0);
    std::cout << std::endl;
    vecs = array2vecvec(jac,3);
    for (auto v : vecs) std::cout << v << std::endl;
    std::cout << "1x2 = " << dot(vecs[0],vecs[1]) << std::endl;
    std::cout << "1x3 = " << dot(vecs[0],vecs[2]) << std::endl;
    std::cout << "3x2 = " << dot(vecs[2],vecs[1]) << std::endl;
#endif
    return 0;
}
