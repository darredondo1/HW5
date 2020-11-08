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
    int numTests = atoi(argv[2]);
    
    //make file
    FILE * fPtr;
    char fPath[40];
    sprintf(fPath,"Problem3/ScatterAllGather/nprocs_%d/ScatterAllGather_N_%d.txt",num_procs,N);
    int numDoubles = 1 << N;
    double* send_message = (double*)malloc(numDoubles*sizeof(double));
    
    //create vector of 2^N random values
    if (rank==0)
    {
        for (int i = 0; i < numDoubles; i++)
        {
            send_message[i] = (double) rand();
        }
    }
    
    int nk = (int) log2(num_procs);
    
    double start, time;
    for (int n=0;n<numTests;n++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        for (int k=1;k<=nk;k++) //NUM STEPS NEEDED TO BROADCAST TO ALL NODES
        {
            int  num_sources = pow(2,k-1); //NUM PROCS AT STEP K THAT HAVE DATA
            int spacing = pow(2,nk-k+1);   //SPACING BETWEEN SOURCES
            int newspacing = pow(2,nk-k);  //HOW FAR TO SEND DATA
            MPI_Status send_status, recv_status;
            MPI_Request send_request, recv_request;
            int send = 0;
            int recv = 0;
            for (int i=0;i<num_sources;i++)
            {
                if (rank==i*spacing) //SENDERS
                {
                    send=1;
                    MPI_Isend(send_message, numDoubles, MPI_DOUBLE, rank+newspacing,k,MPI_COMM_WORLD,&send_request);
                }
                else if (rank==i*spacing+newspacing) //RECEIVERS
                {
                    recv=1;
                    MPI_Irecv(send_message, numDoubles, MPI_DOUBLE,
                              rank-newspacing,k,MPI_COMM_WORLD,&recv_request);
                }
            }
            if (send==1)
                MPI_Wait(&send_request,&send_status);
            if (recv==1)
                MPI_Wait(&recv_request,&recv_status);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        time =  (MPI_Wtime() - start);

        if (rank==0)
        {
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
