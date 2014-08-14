#ifndef _ANALYZER_HPP
#define _ANALYZER_HPP

#include "utils.hpp"

class Analyzer {
public:
    Analyzer(){};
    virtual void addPoint(const state_t & state, const double &time) = 0;
};
#endif
