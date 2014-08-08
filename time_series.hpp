#ifndef _TIME_SERIES_H
#define _TIME_SERIES_H
#include "utils.hpp"
#include <string>
class TimeSeries {
public:
    TimeSeries();
    void addPoint(std::vector<double> &p);
    std::vector<std::vector<double> > &getAll();
    std::vector<double> row(int i);
    void clear() { ts.clear(); }
    void plot(int x, int y);
    void plotRows(std::vector<int> &idx, int xaxis = -1);
    void plotRows();
    ~TimeSeries();
private:
    std::vector<std::string> row_files;
    void createRowFiles(std::vector<int> &idx, int xaxis_row_id = -1);
    void cleanupRowFiles();
    ts_t ts;
};

#endif
