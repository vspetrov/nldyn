#ifndef _MAPPER_H
#define _MAPPER_H
#include "utils.hpp"
#include "analyzer.hpp"
#include <string>


class IBI_Collector : public Analyzer {
private:
    double threshold;
    std::vector<int> stateVarIds;
    state_t prevState;
    double prevTime;
    bool readyFlag;
    std::vector<std::vector<double> > thresholdSections;
public:
    IBI_Collector(double _threshold = 0.0);
    IBI_Collector(std::vector<int> &_stateVarIds, double _threshold = 0.0);
    virtual void addPoint(const state_t &state, const double &time);
    std::vector<double> getIBIs(int _stateVarId);
    virtual void clear() { thresholdSections.clear(); readyFlag = false;}
};
#endif
