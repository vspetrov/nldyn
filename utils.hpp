#ifndef _NLD_UTILS_H
#define _NLD_UTILS_H

#include <vector>
#include <iostream>
#include <assert.h>
#include <numeric>
#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/odeint.hpp>
namespace ublas = boost::numeric::ublas;
namespace ode   = boost::numeric::odeint;

typedef ublas::vector<double> state_t;
typedef ublas::matrix<double> matrix_t;
typedef std::vector<double> ts_point_t;
typedef std::vector<double> ts_row_t;
typedef std::vector<ts_row_t> ts_t;
#endif
