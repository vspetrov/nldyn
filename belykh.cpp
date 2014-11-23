#include <iostream>
#include "systems/phase.hpp"
#include "le_solver.h"
#include "analyzers/time_series.hpp"
#include "analyzers/mapper.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <opencv/cv.h>
#include <sstream>
#include <string>
#include <exception>

#include <pthread.h>
#include <unistd.h>
#include <mpi.h>

#include <sys/fcntl.h>
#include <sys/stat.h>

#define XSTR(X) STR(X)
#define STR(X) #X

#define ALPHA PI/8

template<typename T>
std::string to_str(T value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

#define LOG(rank, size, comm, format, ...) do {                         \
        for (int i=0; i<size; i++) {                                    \
            if (i == rank) {                                            \
                fprintf(stderr, "rank %d: "format"\n", rank, ## __VA_ARGS__); \
                usleep(1000);                                           \
            }                                                           \
        MPI_Barrier(comm);                                              \
        }                                                               \
    }while(0)

static inline
void calc_clusters(std::vector<double> &freqs, std::vector<std::pair<double,int> > &clusters) {
    if (freqs.size() == 0)
        return;

    double v = freqs[0];
    clusters.push_back(std::pair<double,int>(v,0));

    for (auto it=freqs.begin(); it != freqs.end();) {
        if (fabs(*it - v) < 0.01) {
            clusters.back().second++;
            it = freqs.erase(it);
        }else {
            it++;
        }
    }
    calc_clusters(freqs,clusters);
}

static inline
std::string arr_2_str(int *arr, int d1, int d2, int d3) {
    std::string local_rst="[ ";
    for (int i=0; i<d1; i++) {
        for (int j=0; j<d2; j++) {
            local_rst += std::to_string(arr[i*d3+j]) + " ";
        }
        if (i != d1 - 1)
            local_rst += "\n";
    }
    local_rst += "];\n";
    return local_rst;
}
int main(int argc, char **argv) {
    int N = 100;
    double mu = 0.1;



    double K = 0;
    int Na = int(N*mu);
    int Ns = N-Na;
    bool sync = false;
    Kuramoto Kmt(Ns);
    Kmt.initOmegasDeltaRandom(1,1);
    auto om = Kmt.getOmegas();

    K = 2;

    double Kmin=0;
    double Kmax=4;
    double delta_min=0;
    double delta_max=2;

    int ksteps=100;
    int delta_steps=100;

    int rank,size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert((size  % 2) == 0);
    assert((delta_steps % (size/2) == 0));
    assert(size > 0);

    int my_delta_steps = delta_steps/(size/2);
    int my_delta_i = rank/2*my_delta_steps;

    int my_k_i = (rank % 2) * ksteps/2;
    int my_k_steps = ksteps/2;
    int result;

    // LOG(rank, size, MPI_COMM_WORLD, "my_delta_i %d, my_delta_steps %d, my_k_i %d, my_k_steps %d",
    //     my_delta_i, my_delta_steps, my_k_i, my_k_steps);

    int *rst = NULL;
    if ((rank % 2) == 1) {
        rst = (int*) malloc(my_delta_steps*my_k_steps*sizeof(int));
    } else {
        rst = (int*) malloc(my_delta_steps*ksteps*sizeof(int));
    }
    int my_k_offset = ((rank % 2) == 0) ? ksteps : my_k_steps;
    for (int id=my_delta_i; id < my_delta_i + my_delta_steps; id++){
        for (int ik=my_k_i; ik < my_k_i + my_k_steps; ik++)
        {
            K = Kmin + (Kmax-Kmin)/(ksteps-1)*ik;
            double delta = delta_min + (delta_max - delta_min)/(delta_steps-1)*id;
            Kuramoto KmtE(N);
            KmtE.initOmegasExtend(*om, 8, delta); // 8, 0.4 K=4 -  cluster
            // 8, 2   K=4 - chimera

            KmtE.setK(K);
            KmtE.solve(5000,0.01,KmtE.getState());
            std::vector<double> s1(N),s2(N);
            std::copy(KmtE.getState().begin(),KmtE.getState().end(), s1.begin());
            const double TestTime = 1000;
            KmtE.solve(TestTime,0.01,KmtE.getState());
            std::copy(KmtE.getState().begin(),KmtE.getState().end(), s2.begin());

            std::vector<double> freqs(N);
            for (int i=0; i<N; i++) {
                freqs[i] = (s2[i]-s1[i])/TestTime;
            }
            std::vector<std::pair<double,int> > clusters;
            calc_clusters(freqs,clusters);
            int num_larger_2 = 0;
            for (auto &c : clusters)
                if (c.second > 1)
                    num_larger_2++;

            if (clusters.size() == 1) {
                result = 0;
            } else if (clusters.size() == 2) {
                result = 1;
            } else if (num_larger_2 == 1) {
                result = 2;
            } else {
                result = 3;
            }

            int i= id - my_delta_i;
            int j = ik - my_k_i;

            // result = (rank+1)*(id+2)+size*ik*(rank+3);
            rst[i*my_k_offset + j] = result;
        }
    }

    // std::string local_rst=arr_2_str(rst,my_delta_steps,my_k_steps,my_k_offset);
    // LOG(rank ,size, MPI_COMM_WORLD, "\n%s",local_rst.c_str());

    for (int i=0; i<my_delta_steps; i++) {
        if ((rank % 2) == 1) {
            void *buf = (void*)&rst[i*my_k_steps];
            MPI_Send(buf,my_k_steps,MPI_INT,rank-1,1,MPI_COMM_WORLD);
        } else {
            MPI_Status st;
            void *buf = (void*)&rst[my_k_steps+ksteps*i];
            MPI_Recv(buf,my_k_steps,MPI_INT,rank+1,1,MPI_COMM_WORLD,&st);
        }
    }

    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, 0, &comm);

    if ((rank % 2) == 0) {
        // std::string local_rst=arr_2_str(rst,my_delta_steps,ksteps,ksteps);
        // LOG(rank,size/2,comm,"\n%s",local_rst.c_str());
        void *sbuf = rst;
        void *rbuf = NULL;
        if (rank == 0) {
            rbuf = malloc(delta_steps*ksteps*sizeof(int));
        }
        MPI_Gather(sbuf,my_delta_steps*ksteps,MPI_INT,rbuf,my_delta_steps*ksteps, MPI_INT,0,comm);
        if (rank == 0) {
            free(rst);
            rst = (int*)rbuf;
        }
    }


    if (rank == 0) {
    // std::string local_rst = arr_2_str(rst,delta_steps,ksteps,ksteps);
        // std::cout << "RESULT:\n" << local_rst << std::endl;
        int fd = open("rst.bin", O_CREAT | O_RDWR, S_IREAD | S_IWRITE);
        write(fd,rst,delta_steps*ksteps*sizeof(int));
        fd = close(fd);
    }

    if (rst)
        free(rst);

    MPI_Comm_free(&comm);
    MPI_Finalize();

    return 0;
}
