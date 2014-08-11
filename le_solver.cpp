#include "le_solver.h"
#include <algorithm>


static void GramShmidt(int dim, state_t &state, std::vector<double> &norms) {
    auto st = state.begin()+dim;
    norms[0] = L2_norm(st, st+dim);
    normalize(st, st+dim, norms[0]);

    for (int i=1; i<dim; i++) {
        auto vec = st+i*dim;
        state_t m(vec,vec+dim);

        for (int j=0; j<i; j++){
            double v = std::inner_product(m.begin(),m.end(),st+j*dim,0.0);
            for (int k=0; k<dim; k++)
                *(vec+k) -= *(st+j*dim+k)*v;
        }
        norms[i] = L2_norm(vec,vec+dim);
        normalize(vec, vec+dim, norms[i]);
    }
}
std::vector<double> LyapunovExpsSolver::calcLE(double warmUpTime,
                                               double wudt,
                                               int numSteps,
                                               double stepTime,
                                               double dt,
                                               state_t &ini,
                                               double eps) {
    nld_sys->setSolveBoth(false);
    nld_sys->solve(warmUpTime, wudt, ini);
    state_t s0 = nld_sys->getState();
    int dim = nld_sys->getDim();

    bool epsCrit = (eps > 0);

    LEs.resize(dim);
    for (int i=0; i<LEs.size(); i++) {
        LEs[i] = 0.0;
    }
    auto immediateLEs(LEs);
    for (auto &v : immediateLEs) v = 1e10;

    auto diff(LEs);

    nld_sys->setSolveBoth(true);
    int dboth = nld_sys->getDimBoth();
    state_t s_both(dboth);
    std::copy(s0.begin(), s0.end(), s_both.begin());
    std::fill(s_both.begin()+dim, s_both.end(),0.0);
    for (int i=0; i<dim; i++) {
        s_both[dim+i+i*dim] = 1.0;
    }

    TimeSeries ts;

    std::vector<double> norms(dim);
    for (auto &v: norms) v = 0;

    int actualSteps = numSteps;
    for (int i=0; i<numSteps; i++) {
        nld_sys->solve(stepTime, dt, s_both);
        s_both.swap(nld_sys->getState());
        GramShmidt(dim, s_both, norms);

        for (int j=0; j<dim; j++) {
            // std::cout << "normed " << j << ":" << LE_matrix_ortonormed[j] << std::endl;
            LEs[j] += log(norms[j]);
        }

        for (int j=0; j<dim; j++) {
            double le = LEs[j]/((i+1)*stepTime);
            if (epsCrit) {
                double maxv = fabs(le) > fabs(immediateLEs[j]) ? fabs(le) : fabs(immediateLEs[j]);
                diff[j] = fabs(le-immediateLEs[j])/maxv;
            }
            immediateLEs[j] = le;

        }

        if (debugFlag) {
            ts.addPoint(immediateLEs);
        }

        if (epsCrit){
            bool can_stop = true;
            for (auto v : diff)
                if (v > eps) {
                    can_stop = false;
                    actualSteps = i+1;
                    break;
                }
            if (can_stop) {
                std::cout << "Desired accuracy achieved: iterations " <<
                    i << " out of " << numSteps << std::endl;
                break;
            }
        }

    }

    for (int i=0; i<dim; i++) {
        LEs[i] = LEs[i]/(actualSteps*stepTime);
    }
    calcKaplanYorkeDimension(LEs);
    if (debugFlag)
        ts.plotRows();
    return LEs;

}

void LyapunovExpsSolver::calcKaplanYorkeDimension(std::vector<double> &LEs)
{
    KYdim = -1;
    int k = LEs.size();
    std::sort(LEs.begin(),LEs.end(),[=](const double &v1,const double &v2){return v1 > v2; });
    if (LEs[0] < 0) return;
    double sum = std::accumulate(LEs.begin(),LEs.end(),0.0);
    while (sum < 0 && k > 0) {
        sum -= LEs[k-1];
        k--;
    }
    if (k == LEs.size())
        KYdim = -2;
    else
        KYdim = k+sum/fabs(LEs[k]);
}
