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
    int rank, num_procs, sum, sum2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    int N = atoi(argv[1]);
    int numTests = atoi(argv[2]);
    
    //make file
    FILE * fPtr;
    char fPath[100];
    sprintf(fPath,"Problem3/ScatterAllGather/nprocs_%d/ScatterAllGather_N_%d.txt",num_procs,N);
    int numDoubles = 1 << N;
    int blockSize = (int) (numDoubles / num_procs);
    
    //send_message in rank 0 will be the data to be broadcast
    //send_message in all other ranks will update at each step in the ring to eventually hold all of the broadcast data
    double* send_message = (double*)malloc(numDoubles*sizeof(double));
    double* last_message = (double*)malloc(blockSize*sizeof(double));
    
    for (int n=0;n<numTests;n++)
    {
        double start, time;
        
        //create vector of 2^N random values
        if (rank==0)
        {
            double sum=0;
            for (int i = 0; i < numDoubles; i++)
            {
                send_message[i] = (double) rand();
                sum += send_message[i];
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        
        MPI_Scatter(send_message,blockSize,MPI_DOUBLE,last_message,blockSize,  MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        send_message[rank*blockSize] = *last_message;
        
        for (int k=1;k<num_procs;k++)
        {
            MPI_Status send_status, recv_status;
            MPI_Request send_request, recv_request;
            int send_to = (int) (rank+1)%num_procs;
            int recv_from = (int) (rank-1)%num_procs;
            MPI_Isend(last_message,blockSize,MPI_DOUBLE,send_to,k,MPI_COMM_WORLD,&send_request);
            MPI_Irecv(last_message,blockSize,MPI_DOUBLE,recv_from,k,MPI_COMM_WORLD,&recv_request);
            int idx = (int) ((rank+k)*blockSize)%num_procs;
            send_message[idx]=*last_message;
            MPI_Wait(&send_request,&send_status);
            MPI_Wait(&recv_request,&recv_status);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        time =  (MPI_Wtime() - start);
        
        printf("rank %d\n",rank);

        if (rank==1)
        {
            double sum2 = 0;
            for (int i=0;i<numDoubles;n++)
            {
                sum2+=send_message[i];
            }
            printf("sum rank %d %e\n",num_procs,summ);
        }
        
        if (rank==0)
        {
            printf("sum rank 0 %e\n",sum);
            //Save result
            fPtr = fopen(fPath ,"a");
            if (fPtr == NULL) exit(EXIT_FAILURE);
            fprintf(fPtr,"%e\n",time);
            fclose(fPtr);
        }
    }
    MPI_Finalize();

    return (0);
}
