//
//  point_to_point.c
//  
//
//  Created by David Arredondo on 10/13/20.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    int N = atoi(argv[1]);
    int numPingPongs = atoi(argv[2]);
    
    for (int npp = 0; npp < numPingPongs; npp++)
    {
        //create vector of 2^N random values
        int numDoubles = 1 << N;
        double* send_message = (double*)malloc(numDoubles*sizeof(double));
        double* recv_message = (double*)malloc(numDoubles*sizeof(double));
        for (int i = 0; i < numDoubles; i++)
        {
            send_message[i] = (double) rand();
        }
        double start, time;
        MPI_Status recv_status;
        start = MPI_Wtime();
        if (rank%2 == 0)
        {
            MPI_Send(send_message, numDoubles, MPI_DOUBLE, rank+1, 1234, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Recv(recv_message, numDoubles, MPI_DOUBLE, rank+1, 1234, MPI_COMM_WORLD, &recv_status);
            time =  (MPI_Wtime() - start)/2; // divide by two because this is the time for two messages
            //Save result
            FILE * fPtr;
            char fPath[40];
            sprintf(fPath,"Problem1/N_%d.txt",N);
            printf("%s", fPath);
            fPtr = fopen(fPath ,"a");
            if (fPtr == NULL) exit(EXIT_FAILURE);
            fprintf(fPtr,"%d, %e",N, time);
            fclose(fPtr);
        }
        else
        {
            MPI_Recv(recv_message, numDoubles, MPI_DOUBLE, rank-1, 1234, MPI_COMM_WORLD, &recv_status);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Send(send_message, numDoubles, MPI_DOUBLE, rank-1, 1234, MPI_COMM_WORLD);
        }
        free(send_message);
        free(recv_message);
    }
    MPI_Finalize();
    if (rank==0) printf("Avg ping-pong %e\n",avgTime);
    
    //Save the average time to a file
    if (rank == 0)
    {

    }

    return (0);
}
