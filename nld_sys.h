#ifndef _NLD_SYS
#define _NLD_SYS
#include <vector>
#include <iostream>
#include <assert.h>
#include <numeric>
#include <math.h>

typedef std::vector<double> state_t;
typedef std::vector<state_t> ts_t;
typedef std::vector<double>::iterator it_t;
typedef std::vector<double>::const_iterator const_it_t;
class TimeSeries {
public:
    TimeSeries();
    void addPoint(std::vector<double> &p);
    std::vector<std::vector<double> > &getAll();
    std::vector<double> row(int i);
    void clear() { ts.clear(); }
    void plot(int x, int y);
    void plotRows(std::vector<int> &idx);
    void plotRows();
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
    std::vector<state_t> jacobian;
    bool solveBoth;
public:
    System(int dimension);
    virtual void rhs(const_it_t &state,
                     it_t & out, double time) = 0;

    void rhs_linearized(const_it_t &state,
                        const_it_t & lin_state,
                        it_t & out, double time);

    virtual void jac(const_it_t &state,
                     std::vector<state_t> &out, double time) = 0;
    void rk4_step(state_t &v, double dt, double time);
    void rhs_dt(const state_t &state, state_t & out, double time, double dt);
    void solve(double MaxTime, double dt,
          state_t ini,
          double saveTS);
    void solve(double MaxTime, double dt, state_t ini);
    TimeSeries & getTs();
    int getDim() { return dim; }
    state_t & getState() { return vars; }
    void setSolveBoth(bool value);
    int getDimBoth() { return dim*(dim+1);}
};

class Lorenz : public System {
private:
    double sigma;
    double ro;
    double beta;
public:
    Lorenz();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     std::vector<state_t> & out, double time);
};

class Rossler : public System {
private:
    double a;
    double b;
    double c;
public:
    Rossler();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     std::vector<state_t> & out, double time);
};

class FHN : public System {
private:
    double a;
    double e;
public:
    FHN();
    virtual void rhs(const_it_t &state, it_t & out, double time);
    virtual void jac(const_it_t &state,
                     std::vector<state_t> & out, double time);
};

inline std::vector<double>  operator/(const std::vector<double> &lhs, const double &rhs) {
    std::vector<double> result(lhs.size());
    for (int i=0; i<lhs.size(); i++) result[i] = lhs[i]/rhs;
    return result;
}

inline std::vector<double>  operator*(const std::vector<double> &lhs, const double &rhs) {
    std::vector<double> result(lhs.size());
    for (int i=0; i<lhs.size(); i++) result[i] = lhs[i]*rhs;
    return result;
}

inline std::vector<double>  operator+(const std::vector<double> &lhs, const std::vector<double> &rhs) {
    std::vector<double> result(lhs.size());
    for (int i=0; i<lhs.size(); i++) result[i] = lhs[i] + rhs[i];
    return result;
}

inline std::vector<double>  operator-(const std::vector<double> &lhs, const std::vector<double> &rhs) {
    std::vector<double> result(lhs.size());
    for (int i=0; i<lhs.size(); i++) result[i] = lhs[i] - rhs[i];
    return result;
}

inline void sum_2vec(std::vector<double> &v1, double c1,
                     std::vector<double> &v2, double c2,
                     std::vector<double> &rst) {
    for (std::vector<double>::iterator it1 = v1.begin(), it2 = v2.begin(), it_rst = rst.begin();
         it1 != v1.end(); it1++, it2++, it_rst++)
        *it_rst = (*it1)*c1+(*it2)*c2;
}

inline void sum_2vec(std::vector<double> &v1,
                     std::vector<double> &v2, double c2,
                     std::vector<double> &rst) {
    for (std::vector<double>::iterator it1 = v1.begin(), it2 = v2.begin(), it_rst = rst.begin();
         it1 != v1.end(); it1++, it2++, it_rst++)
        *it_rst = (*it1)+(*it2)*c2;
}

inline void sum_2vec(std::vector<double> &v1,
                     std::vector<double> &v2,
                     std::vector<double> &rst) {
    for (std::vector<double>::iterator it1 = v1.begin(), it2 = v2.begin(), it_rst = rst.begin();
         it1 != v1.end(); it1++, it2++, it_rst++)
        *it_rst = (*it1)+(*it2);
}

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
inline double dot(const std::vector<double> &v1, const std::vector<double> &v2) {
    return std::inner_product(v1.begin(),v1.end(),v2.begin(),0.0);
}

inline double L2_norm(const const_it_t &start, const const_it_t &end) {
    return sqrt(std::inner_product(start,end,start,0.0));
}

inline void normalize(it_t &start, const const_it_t &end, double norm) {
    for (int i=0; start+i != end; i++)
        start[i] /= norm;
}
#endif
