#include "mapper.hpp"

IBI_Collector::IBI_Collector(double _threshold) {
    threshold = _threshold;
    prevTime = -1;
    readyFlag = false;
}

void IBI_Collector::addPoint(const state_t &state, const double &time) {
    for (int i=0; i<stateVarIds.size(); i++) {
        int id = stateVarIds[i];
        if (readyFlag) {
            if ((prevState[id] - threshold) < 0 &&
                (state[id] - threshold) > 0) {
                thresholdSections[i].push_back(time);
            }
            prevState[id] = state[id];
        } else {
            prevState = state;
            readyFlag = true;
        }

    }
}

IBI_Collector::IBI_Collector(std::vector<int> &_stateVarIds, double _threshold) {
    threshold = _threshold;
    stateVarIds = _stateVarIds;
    prevTime = -1;
    readyFlag  = false;
    thresholdSections.resize(_stateVarIds.size());
}


std::vector<double> IBI_Collector::getIBIs(int _stateVarId) {
    std::vector<double> ibis;
    for (int i=0; i<stateVarIds.size(); i++) {
        if (stateVarIds[i] == _stateVarId) {
            auto thresholdSectionsVarId = thresholdSections[i];
            for (int j=1; j<thresholdSectionsVarId.size(); j++) {
                ibis.push_back(thresholdSectionsVarId[j] -
                               thresholdSectionsVarId[j-1]);
            }
        }
    }
    return ibis;
}
