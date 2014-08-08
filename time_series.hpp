#ifndef _TIME_SERIES_H
#define _TIME_SERIES_H
#include "utils.hpp"

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

#endif
