#include "time_series.hpp"
#include <string>
#include <sstream>
#include <stdio.h>

TimeSeries::TimeSeries() {
}

void TimeSeries::addPoint(std::vector<double> &p) {
    ts.push_back(p);
}

std::vector<std::vector<double> >& TimeSeries::getAll() {
    return ts;
}
std::vector<double> TimeSeries::row(int i) {
    std::vector<double> row;
    for (auto point : ts) row.push_back(point[i]);
    return row;
}

void TimeSeries::plot(int x, int y) {
    auto xr = row(x);
    auto yr = row(y);
    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp,"plot '-' w l\n");
    assert(xr.size() == yr.size());
    for (int i=0; i<xr.size(); i++) {
        fprintf(gp,"%g %g\n",xr[i],yr[i]);
    }
    fprintf(gp,"e");
    fclose(gp);
}

void TimeSeries::plotRows(std::vector<int> &idx) {
    std::vector<std::string> files;
    std::string gp_cmd="plot ";
    for (auto id : idx) {
        auto r = row(id);
        std::ostringstream oss;
        oss << id;
        std::string filename = ".nldyn.ts.row"+oss.str() + ".dat";
        files.push_back(filename);
        FILE *f = fopen(filename.c_str(),"w");
        for (auto v : r)
            fprintf(f,"%g\n",v);
        fclose(f);
    }
    for (int i=0; i<files.size() -1 ; i++)
        gp_cmd += "'"+files[i]+"' w l, ";
    gp_cmd += "'"+files.back()+"' w l\n";

    FILE *gp = popen("gnuplot -persist","w");
    fprintf(gp,"%s",gp_cmd.c_str());
    fprintf(gp,"e");
    fclose(gp);

    for (auto f : files) remove(f.c_str());
}

void TimeSeries::plotRows() {
    std::vector<int> idx;
    for (int i=0; i<ts[0].size(); i++)
        idx.push_back(i);
    plotRows(idx);
}
TimeSeries::~TimeSeries() {
}
