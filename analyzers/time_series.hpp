#ifndef _TIME_SERIES_H
#define _TIME_SERIES_H
#include "utils.hpp"
#include "analyzer.hpp"
#include <string>

class TimeSeries : public Analyzer{
public:
    TimeSeries(double interval);
    TimeSeries();
    virtual void addPoint(const state_t &state, const double &time);
    void addPoint(const ts_row_t &point);
    ts_t &getAll();
    ts_row_t row(int i);
    virtual void clear() { ts.clear(); }
    void plot(int x, int y);
    void plotRows(std::vector<int> &idx, int xaxis = -1);
    void plotRows();
    ~TimeSeries();
    void setNoLegend(bool v) { noLegend = v; }
private:
    double m_lastPointTime;
    double m_interval;
    bool noLegend;
    std::vector<std::string> row_files;
    void createRowFiles(std::vector<int> &idx, int xaxis_row_id = -1);
    void cleanupRowFiles();
    ts_t ts;
};

#endif
