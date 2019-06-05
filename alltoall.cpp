//
// Created by Page, Brian Andrew on 2019-06-04.
//

#include "alltoall.h"


int main (int argc, char *argv[]){
    int i, numprocs, rank, maxMessageSize = 1048576;
    double latency = 0.0, t_start = 0.0, t_stop = 0.0;
    double timer=0.0;
    double avg_time = 0.0, max_time = 0.0, min_time = 0.0;
    char * sendbuf = NULL, * recvbuf = NULL;
    int po_ret;
    size_t bufsize;
    int iterations = 10;
    int windowSize = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


    if(numprocs < 2) {
        if (rank == 0) {
            std::cout << "This test requires at least two processes" << std::endl;
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    std::vector <int> sendbuffer (maxMessageSize, rand());
    std::vector <int> receiveBuffer (maxMessageSize);

    for(int size = 1; size <= maxMessageSize; size *= 2) {
        MPI_Barrier(MPI_COMM_WORLD);
        timer = 0.0;

        for (int i=0; i < iterations; i++) {
            t_start = MPI_Wtime();
            MPI_Alltoall(sendbuffer.data(), size, MPI_INT, receiveBuffer.data(), size, MPI_INT, MPI_COMM_WORLD);
            t_stop = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);
            timer += t_stop - t_start;
        }
        latency = (double)(timer * 1e6) / iterations;

        MPI_Reduce(&latency, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&latency, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&latency, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        avg_time = avg_time/numprocs;

        double tmp = (size / 1e6) * numprocs;
        tmp *= iterations * windowSize;
        double bandwidth = tmp /(avg_time);

        if (rank == 0){
            std::cout << numprocs << "," << size << "," << min_time << "," << max_time << "," << avg_time << "," << bandwidth << std::endl;
        }
    }

    MPI_Finalize();

    return 0;
}

