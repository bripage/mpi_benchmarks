//
// Created by Page, Brian Andrew on 2019-06-04.
//

#include "alltoall.h"


int main (int argc, char *argv[]){
    int i, numprocs, rank, maxMessageSize;
    double latency = 0.0, t_start = 0.0, t_stop = 0.0;
    double timer=0.0;
    double avg_time = 0.0, max_time = 0.0, min_time = 0.0;
    char * sendbuf = NULL, * recvbuf = NULL;
    int po_ret;
    size_t bufsize;

    MPI_CHECK(MPI_Init(&argc, &argv));
    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &numprocs));


    if(numprocs < 2) {
        if (rank == 0) {
            std::cout << "This test requires at least two processes" << std::endl;
        }
        MPI_CHECK(MPI_Finalize());
        exit(EXIT_FAILURE);
    }

    std::vector <int> sendbuffer (maxMessageSize, rand());
    std::vector <int> receiveBuffer (maxMessageSize);

    for(int size = 1; size <= maxMessageSize; size *= 2) {
        MPI_Barrier(MPI_COMM_WORLD);
        timer = 0.0;

        for (int i=0; i < 10; i++) {
            t_start = MPI_Wtime();
            MPI_Alltoall(sendbuffer.data(), size, MPI_INT, receiveBuffer.data(), size, MPI_INT, MPI_COMM_WORLD);
            t_stop = MPI_Wtime();

            if (i >= options.skip) {
                timer+=t_stop-t_start;
            }
            MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
        }
        latency = (double)(timer * 1e6) / options.iterations;

        MPI_CHECK(MPI_Reduce(&latency, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0,
                             MPI_COMM_WORLD));
        MPI_CHECK(MPI_Reduce(&latency, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0,
                             MPI_COMM_WORLD));
        MPI_CHECK(MPI_Reduce(&latency, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0,
                             MPI_COMM_WORLD));
        avg_time = avg_time/numprocs;

        print_stats(rank, size, avg_time, min_time, max_time);
        MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
    }

    free_buffer(sendbuf, options.accel);
    free_buffer(recvbuf, options.accel);

    MPI_CHECK(MPI_Finalize());

    if (NONE != options.accel) {
        if (cleanup_accel()) {
            fprintf(stderr, "Error cleaning up device\n");
            exit(EXIT_FAILURE);
        }
    }

    return EXIT_SUCCESS;
}

/* vi: set sw=4 sts=4 tw=80: */
