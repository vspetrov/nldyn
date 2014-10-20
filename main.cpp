#include <iostream>
#include "systems/fhn.hpp"
#include "systems/lorenz.hpp"
#include "systems/ap.hpp"
#include "le_solver.h"
#include "analyzers/time_series.hpp"
#include "analyzers/mapper.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>

#include <opencv/cv.h>
int main(int argc, char **argv) {
    System *s;
#if 0
    FHN3 fhn3;
    fhn3.setDeo(0.0);
    fhn3.setDoe(0.02);
    fhn3.setDoo(0.0);
    s = &fhn3;
    s->addAnalyzer(new TimeSeries(0.1));
    s->solve(10000, 0.01, s->getState());
    std::vector<int> rows_to_plot;
    for (int i=0; i<3; i++) rows_to_plot.push_back(i*2+1);
    s->getAnalyzer<TimeSeries>().back()->plotRows(rows_to_plot,0);
    s->clearAnalyzers();
    s->addAnalyzer(new ConstSectionMapper(4,1.0));
    s->solve(50000,0.01,s->getState());
    auto ibi = s->getAnalyzer<ConstSectionMapper>().back()->getMap(4);
    for (int i=0; i<ibi->size(); i++)
        std::cout << ibi->at(i).first << std::endl;
#endif
#if 1
    std::string rst_dir_prefix = ".";
    if (argc > 1)
        rst_dir_prefix = std::string(argv[1]);

    std::string filenames[3];
    filenames[0] = rst_dir_prefix + "/ibi_xo1.dat";
    filenames[1] = rst_dir_prefix + "/ibi_xo2.dat";
    filenames[2] = rst_dir_prefix + "/ibi_xe.dat";

    std::ofstream stddev_ofs(rst_dir_prefix + "/ibi_stdev.dat");
    std::ofstream ofs_clustering(rst_dir_prefix + "/ibi_clusters.dat");
    std::ofstream ofs[3];
    for (int i=0; i<3; i++)
        ofs[i].open(filenames[i]);


    double dmin = 0;
    double dmax = 0.02;
    int steps = 100;
    double stepsize = (dmax-dmin)/(steps-1);
    for (double Doo = dmin; Doo < dmax+1E-6; Doo+=stepsize){ //0.00005
        FHN3 fhn3;
        fhn3.setDeo(0.01);
        fhn3.setDoe(0.1);
        fhn3.setDoo(Doo);
        fhn3.setEpsilon(2,0.005);
        s = &fhn3;
        std::cout << "Doo = " << Doo << std::endl;
        std::cout << "Transient skip maps" << std::endl;
        s->solve(20000,0.01,s->getState());
        typedef ConstSectionMapper MT;
        s->addAnalyzer(new MT(0,1.0));
        s->addAnalyzer(new MT(2,1.0));
        s->addAnalyzer(new MT(4,1.0));
        std::cout << "Calculating IBIs" << std::endl;
        s->solve(500000,0.01,s->getState());
        auto mappers = s->getAnalyzer<MT>();
        // std::ostringstream oss;
        // oss << Doo;
        // std::string prefix = rst_dir_prefix + "/ibi_Doo_"+oss.str();
        // mappers[0]->saveMapIntervals(4, prefix+"_xe.txt");
        // mappers[1]->saveMapIntervals(0, prefix+"_xo0.txt");
        // mappers[2]->saveMapIntervals(2, prefix+"_xo1.txt");
        std::vector<double> stddevs;
        for (int j=0; j<mappers.size(); j++){
            auto map = mappers[j]->getMap(j*2);

            double sum = 0, sq_sum = 0;
            for (int i=0; i<map->size(); i++) {
                double v = (*map)[i].first;
                ofs[j] << Doo << " " << v << std::endl;
                sum += v;
                sq_sum += v*v;
            }
            double mean = sum / map->size();
            double stdev = sqrt(sq_sum / map->size() - mean*mean);
            stddevs.push_back(stdev);
        }
        stddev_ofs << Doo;
        for (int i=0; i<3; i++) {
            stddev_ofs << " " << stddevs[i];
        }
        stddev_ofs << std::endl;


        // clustering
        auto map = mappers[2]->getMap(4);
        cv::Mat data(map->size(),1,CV_32FC1), labels;
        for (int i=0; i<map->size(); i++) {
            data.at<float>(i,0) =(float) (*map)[i].first;
        }
        cv::kmeans(data, 2, labels,
                   cv::TermCriteria(cv::TermCriteria::EPS+cv::TermCriteria::COUNT,
                                    20, 0.1), 2, cv::KMEANS_RANDOM_CENTERS);
        std::vector<float> clusters[2];
        for (int i=0; i<map->size(); i++) {
            int id = labels.at<int>(i);
            clusters[id].push_back(data.at<float>(i,0));
        }

        double cmins[2], cmaxs[2];
        double centers[2];

        for (int i=0; i<2; i++) {
            cv::minMaxLoc(clusters[i], &cmins[i], &cmaxs[i], NULL, NULL);
            centers[i] = std::accumulate(clusters[i].begin(), clusters[i].end(), 0.0)/clusters[i].size();
        }
        int large_id = clusters[0] > clusters[1] ? 0 : 1;

        ofs_clustering << Doo << " " << cmins[large_id] << " " << cmaxs[large_id] << " "
                       << centers[large_id] << " " << cmins[(large_id + 1) % 2] << " " << cmaxs[(large_id + 1) % 2] << " " <<
            centers[(large_id + 1) % 2] << std::endl;
        s->clearAnalyzers();
    }
    ofs[0].close();
    ofs[1].close();
    ofs[2].close();
    stddev_ofs.close();

    ofs_clustering.close();
    std::cout << "done" << std::endl;
#endif
    return 0;
}
