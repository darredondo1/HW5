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
    double sum, sum2;
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
    
    //allData in rank 0 will be the data to be broadcast
    //allData in all other ranks will begin empty and update at each step in the ring to eventually hold all of the broadcast data
    double* allData = (double*)malloc(numDoubles*sizeof(double));
    double* send_message = (double*)malloc(blockSize*sizeof(double));
    double* recv_message = (double*)malloc(blockSize*sizeof(double));
    
    for (int n=0;n<numTests;n++)
    {
        double start, time;
        
        //create vector of 2^N random values
        if (rank==0)
        {
            sum=0;
            for (int i = 0; i < numDoubles; i++)
            {
                allData[i] = (double) rand();
                sum += allData[i];
            }
        }
        
        //START TIMING
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        
        //SCATTER FROM PROCESS 0
        MPI_Scatter(allData,blockSize,MPI_DOUBLE,send_message,blockSize,  MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        //UPDATE allData WITH VALUES RECEIVED FROM SCATTER
        for (int j = 0; j < blockSize; j++)
        {
            allData[rank*blockSize + j] = send_message[j];
        }
        
        //RING ALGORITHM
        for (int k=1;k<num_procs;k++)
        {
            MPI_Status send_status, recv_status;
            MPI_Request send_request, recv_request;
            int send_to = (int) (rank+1)%num_procs;
            int recv_from = (int) (rank+num_procs-1)%num_procs;
            MPI_Isend(send_message,blockSize,MPI_DOUBLE,send_to,k,MPI_COMM_WORLD,&send_request);
            MPI_Irecv(recv_message,blockSize,MPI_DOUBLE,recv_from,k,MPI_COMM_WORLD,&recv_request);
            MPI_Wait(&send_request,&send_status);
            MPI_Wait(&recv_request,&recv_status);
            
            //update send_message with the latest message
            //update allData with the latest message
            int idx = (int) ((rank+num_procs-k)%num_procs)*blockSize;
            for (int j = 0; j < blockSize; j++)
            {
                send_message[j] = recv_message[j];
                allData[idx+j] = send_message[j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        time =  (MPI_Wtime() - start);
        
        //print sums of rank 1 and 0 to check that the broadcast was successful
        if (rank==1)
        {
            sum2=0;
            for (int i=0;i<numDoubles;i++)
            {
                sum2 += allData[i];
            }
            printf("rank %d sum %e\n",rank,sum2);
        }
        
        if (rank==0)
        {
            printf("rank %d sum %e\n",rank,sum);
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
