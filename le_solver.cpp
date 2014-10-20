#include "le_solver.h"
#include <algorithm>
#include "analyzers/time_series.hpp"

void LyapunovExpsSolver::GramShmidt(ublas::vector_range<state_t> &state, std::vector<double> &norms) {
    auto st = ublas::subrange(state,0,dim);
    norms[0] = ublas::norm_2(st);
    st /= norms[0];

    for (int i=1; i<dim; i++) {
        auto vec = ublas::subrange(state,i*dim,(i+1)*dim);
        for (int j=0; j<i; j++){
            projections[j] = ublas::inner_prod(vec,ublas::subrange(state,j*dim,(j+1)*dim));
            // for (int k=0; k<dim; k++)
                // *(vec+k) -= *(st+j*dim+k)*v;
        }
        for (int j=0; j<i; j++) {
            vec -= projections[j]*ublas::subrange(state,j*dim,(j+1)*dim);
        }
        norms[i] = ublas::norm_2(vec);
        vec /= norms[i];
    }
}
std::vector<double> & LyapunovExpsSolver::calcLE(double warmUpTime,
                                                 double wudt,
                                                 int numSteps,
                                                 double stepTime,
                                                 double dt,
                                                 state_t &ini,
                                                 double eps) {
    nld_sys->setSolveCombined(false);
    nld_sys->solve(warmUpTime, wudt, ini);
    state_t s0 = nld_sys->getState();
    bool epsCrit = (eps > 0);

    LEs.resize(dim);
    for (int i=0; i<LEs.size(); i++) {
        LEs[i] = 0.0;
    }
    auto immediateLEs(LEs);
    for (auto &v : immediateLEs) v = 1e10;

    auto diff(LEs);

    nld_sys->setSolveCombined(true);
    int dboth = dim*(dim+1);
    state_t s_both(dboth);
    std::copy(s0.begin(), s0.end(), s_both.begin());
    std::fill(s_both.begin()+dim, s_both.end(),0.0);
    for (int i=0; i<dim; i++) {
        s_both[dim+i+i*dim] = 1.0;
    }

    m_ImmLes.clear();

    std::vector<double> norms(dim);
    for (auto &v: norms) v = 0;

    int actualSteps = numSteps;
    for (int i=0; i<numSteps; i++) {
        nld_sys->solve(stepTime, dt, s_both);
        s_both.swap(nld_sys->getState());
        auto vecs =ublas::subrange(s_both,dim,dboth);
        GramShmidt(vecs, norms);

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

        if (keepImmediateLEs) {
            std::vector<double> ln;
            for (auto &v : norms) ln.push_back(log(v));
            ln.push_back(stepTime*(i+1));
            m_ImmLes.push_back(ln);
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
    calcKaplanYorkeDimension();
    // if (debugFlag) {
        // ts.setNoLegend(true);
        // ts.plotRows();
    // }
    return LEs;

}

void LyapunovExpsSolver::plotFiniteTimeLEs(double window)
{
    double time = m_ImmLes[0].back();
    auto imm = m_ImmLes.begin();
    TimeSeries ts;
    std::vector<double> immLes(imm->size() - 1);
    std::fill(immLes.begin(), immLes.end(), 0.0);

    while (imm != m_ImmLes.end()) {
        double t = imm->back();
        for (int i=0; i<immLes.size() - 1; i++) {
            immLes[i] += imm->at(i);
        }
        if ( t - time > window - 1e-10) {
            for (int i=0; i<immLes.size() - 1; i++) immLes[i] /= window;
            ts.addPoint(immLes);
            std::fill(immLes.begin(), immLes.end(), 0.0);
            time = t;
        }
        imm++;
    }
    ts.setNoLegend(true);
    ts.plotRows();
}


void LyapunovExpsSolver::calcKaplanYorkeDimension()
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

LyapunovExpsSolver::LyapunovExpsSolver(System *s)
{
        nld_sys = s;
        keepImmediateLEs = false;
        dim = s->getDim();
        projections.resize(dim);
}
