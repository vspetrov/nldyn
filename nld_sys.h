#ifndef _NLD_SYS
#define _NLD_SYS
#include <vector>
#include <iostream>
#include <assert.h>
#include <numeric>
#include <math.h>

typedef std::vector<double> state_t;
typedef std::vector<state_t> ts_t;
class TimeSeries {
public:
    TimeSeries();
    void addPoint(std::vector<double> &p);
    std::vector<std::vector<double> > &getAll();
    std::vector<double> row(int i);
    void clear() { ts.clear(); }
    void plot(int x, int y);
    ~TimeSeries();
private:
    ts_t ts;
};


class System {
protected:
    int dim;
    state_t vars;
    state_t _rk4[4];
    TimeSeries ts;
    std::vector<std::vector<double> > jacobian;
public:
    System(int dimension);
    virtual void rhs(std::vector<double> &state, std::vector<double> & out, double time) = 0;
    virtual void jac(std::vector<double> &state,
                     std::vector<std::vector<double> >& out, double time) = 0;
    void rk4_step(double dt, double time, state_t *state_for_jacobian = NULL);
    void rhs_dt(std::vector<double> &state, std::vector<double> & out, double time, double dt, state_t *state_for_jacobian = NULL);
    void rhs_linearized(std::vector<double> &state, state_t &lin_state, std::vector<double> & out, double time);
    void solve(double MaxTime, double dt,
          state_t ini,
          double saveTS = -1);
    void solve_linearized(ts_t &ts, state_t &ini);
    TimeSeries & getTs();
    int getDim() { return dim; }
    state_t getState() { return vars; }
};

class Lorenz : public System {
private:
    double sigma;
    double ro;
    double beta;
public:
    Lorenz(int dimension);
    virtual void rhs(std::vector<double> &state, std::vector<double> & out, double time);
    virtual void jac(std::vector<double> &state,
                     std::vector<std::vector<double> > & out, double time);
};

class Rossler : public System {
private:
    double a;
    double b;
    double c;
public:
    Rossler();
    virtual void rhs(std::vector<double> &state, std::vector<double> & out, double time);
    virtual void jac(std::vector<double> &state,
                     std::vector<std::vector<double> > & out, double time);
};

class FHN : public System {
private:
    double a;
    double e;
public:
    FHN();
    virtual void rhs(std::vector<double> &state, std::vector<double> & out, double time);
    virtual void jac(std::vector<double> &state,
                     std::vector<std::vector<double> > & out, double time);
};

inline std::vector<double>  operator/(const std::vector<double> lhs, const double rhs) {
    std::vector<double> result(lhs.size());
    for (int i=0; i<lhs.size(); i++) result[i] = lhs[i]/rhs;
    return result;
}

inline std::vector<double>  operator+(const std::vector<double> lhs, const std::vector<double> rhs) {
    std::vector<double> result(lhs.size());
    for (int i=0; i<lhs.size(); i++) result[i] = lhs[i] + rhs[i];
    return result;
}

// inline std::ostream& operator<<(std::ostream& os, std::vector<double> &v) {
    // for (int i=0; i<v.size(); i++)  os << v[i] << ' ';
    // return os;
// }
inline std::ostream& operator<<(std::ostream& os, std::vector<double> v) {
    for (int i=0; i<v.size(); i++)  os << v[i] << ' ';
    return os;
}

inline std::vector<std::vector<double> > array2vecvec(std::vector<double> &arr, int dim) {
    std::vector<std::vector<double> > rst;
    assert(arr.size() % dim == 0);
    for (int i=0; i<arr.size() / dim; i++) {
        rst.push_back(std::vector<double>(arr.begin()+dim*i,arr.begin()+dim*i+dim));
    }
    return rst;
}
inline double dot(std::vector<double> &v1, std::vector<double> &v2) {
    return std::inner_product(v1.begin(),v1.end(),v2.begin(),0.0);
}
inline double L2_norm(std::vector<double> &v1) {
    std::vector<double> v2(v1);
    return sqrt(std::inner_product(v1.begin(),v1.end(),v2.begin(),0.0));
}
#endif
