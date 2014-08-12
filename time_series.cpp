#include "time_series.hpp"
#include <string>
#include <sstream>
#include <stdio.h>

TimeSeries::TimeSeries() {
    noLegend = false;
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

void TimeSeries::createRowFiles(std::vector<int> &idx, int xaxis_row_id) {
    state_t xaxis;
    if (xaxis_row_id >= 0 )
        xaxis = row(xaxis_row_id);

    for (auto id : idx) {
        auto r = row(id);
        std::ostringstream oss;
        oss << id;
        std::string filename = ".nldyn.ts.row"+oss.str() + ".dat";
        row_files.push_back(filename);
        FILE *f = fopen(filename.c_str(),"w");
        if (xaxis_row_id >= 0) {
            for (int i=0; i<xaxis.size(); i++)
                fprintf(f,"%g %g\n",xaxis[i],r[i]);
        } else {
            for (auto v : r)
                fprintf(f,"%g\n",v);
        }
        fclose(f);
    }
}

void TimeSeries::cleanupRowFiles() {
    for (auto f : row_files) remove(f.c_str());
}

void TimeSeries::plotRows(std::vector<int> &idx, int xaxis) {
    createRowFiles(idx,xaxis);
    std::string extra_modifier = "";
    if (xaxis >= 0)
        extra_modifier = " u 1:2";
    std::string gp_cmd="plot ";
    for (int i=0; i<row_files.size() -1 ; i++)
        gp_cmd += "'"+row_files[i]+"'"+extra_modifier+" w l, ";
    gp_cmd += "'"+row_files.back()+"'"+extra_modifier+" w l\n";

    FILE *gp = popen("gnuplot -persist","w");
    if (noLegend) {
        fprintf(gp,"set key off\n");
    }
    fprintf(gp,"%s",gp_cmd.c_str());
    fprintf(gp,"e");
    fclose(gp);
    cleanupRowFiles();
}

void TimeSeries::plotRows() {
    std::vector<int> idx;
    for (int i=0; i<ts[0].size(); i++)
        idx.push_back(i);
    plotRows(idx);
}



TimeSeries::~TimeSeries() {
}
