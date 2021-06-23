
#include<stdio.h>
#include <stdlib.h>
#include<string.h>
#include "mpi.h"
#include <time.h>
#include <math.h>
//using namespace std;


void delay(int number_of_seconds)
{
    // Converting time into milli_seconds
    int milli_seconds = 1000 * number_of_seconds;

    // Storing start time
    clock_t start_time = clock();

    // looping till required time is not achieved
    while (clock() < start_time + milli_seconds);
}


//this is the default code of MPI_Bcast

double MPI_Bcast_Default(int argc, char *argv[], int myrankrank, int size)
{
  int count = atoi(argv[1]); 
  double buf[count];
  double sTime, eTime, time; 

  for (int i=0; i<count; i++)
      buf[i] = (double)rand() / (double)RAND_MAX;

  // has to be called by all processes
  sTime = MPI_Wtime();
  MPI_Bcast(buf, count , MPI_DOUBLE, 0, MPI_COMM_WORLD);
  eTime = MPI_Wtime();
  time = eTime - sTime;

 
    MPI_Barrier (MPI_COMM_WORLD);
  return time;
}


//This is the default code of MPI_Reduce


double MPI_Reduce_Default(int argc, char *argv[], int myrank, int size)
{
 // printf("this is reduce default %d",n);
  int count =  atoi(argv[1]);
  //int myrank, size; 
  double sendval[count], maxval[count];
  MPI_Status status;
  double sTime, eTime, time;


  for (int i=0; i<count; i++)
  sendval[i] = (double)rand() / (double)RAND_MAX;
  sTime = MPI_Wtime();
  // reduction of 1 integer
  //sendval = myrank; // initialization
  MPI_Reduce(sendval, maxval, count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // find max of sendvals
//  printf ("%d reduce %d\n", myrank, maxval);
  eTime = MPI_Wtime();
  time = eTime - sTime;

  return time;
}



//This is the default code for MPI_Gather


double MPI_Gather_Default(int argc, char *argv[], int rank , int numtasks)
{

  double sTime,eTime, time;

 

  // Allocate message
  int arrSize = atoi(argv[1]);
  //printf("%d\n",arrSize);
  double *message = (double *) malloc(arrSize * sizeof(double));
 
  for (int i = 0; i < arrSize; i++) {
    message[i] = (double)rand() / (double)RAND_MAX;
  }


  double *recvMessage = (double *) malloc(numtasks * arrSize * sizeof(double)); //significant at the root process
  sTime = MPI_Wtime();
  MPI_Gather(message, arrSize, MPI_DOUBLE, recvMessage, arrSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  eTime = MPI_Wtime();
  time = eTime - sTime;


  return time;
}


//Default code for MPI_Alltoall

double MPI_Alltoallv_Default(int argc, char *argv[], int rank, int numTasks)
{
 //int rank, numTasks;
  double sTime,eTime,time;
  
  // Allocate message
  int totalSize = atoi(argv[1]);
  int arrSize = (totalSize/numTasks);
  double *message = (double *) malloc(totalSize * sizeof(double)); // every process sends to every other process

  // Initialize array 
  for (int i = 0; i < arrSize * numTasks; i++) {
    message[i] = (double)rand() / (double)RAND_MAX;
  }

  // every process receives arrSize elements from other processe
  double *recvMessage = (double *) malloc(totalSize * sizeof(double));
  sTime = MPI_Wtime();
  MPI_Alltoall(message, arrSize, MPI_DOUBLE, recvMessage, arrSize, MPI_DOUBLE, MPI_COMM_WORLD);
  eTime = MPI_Wtime();
  time = eTime - sTime;


  // Verify 
  MPI_Barrier (MPI_COMM_WORLD);

  return time;
}


//Optmised code for MPI_Bcast

double MPI_Bcast_Opt(int argc, char *argv[], int rank, int num_procs){
  
  //declaration of varaiables
  int src, dst, mask, relative_rank;
  int count = atoi(argv[1]);
  double eTime,sTime,time;
  double buff[count];

  int root = 0; 
 
 
  // Initialize buffer
  if(rank == 0){
  for (int i = 0; i < count; i++) {
    buff[i] = (double)rand() / (double)RAND_MAX;
  }
  }

  //getting relative ranks of process
  relative_rank = (rank >= root) ? rank - root : rank - root + num_procs;
  sTime = MPI_Wtime();
  
  mask = 0x1;
  
  //source rank for receiving process
  while (mask < num_procs) {
    if (relative_rank & mask) {
      src = rank - mask;
      if (src < 0)
        src += num_procs;
      MPI_Recv(buff, count, MPI_DOUBLE , src, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      break;
    }
    mask <<= 1;
  }
  //destination rank for sending process
  mask >>= 1;
  while (mask > 0) {
    if (relative_rank + mask < num_procs) {
      dst = rank + mask;
      if (dst >= num_procs)
        dst -= num_procs;
      MPI_Send(buff, count, MPI_DOUBLE , dst, 99, MPI_COMM_WORLD);
    }
    mask >>= 1;
  }
  eTime = MPI_Wtime();
  //printf ("\nthis is opt %d %lf %lf %lf\n", rank, eTime - sTime, buff[0], buff[1]);
  time = eTime - sTime;

  //MPI_Finalize();
  MPI_Barrier (MPI_COMM_WORLD);
  //returning time taken by process in the operation
  return time;
}



//Optimized code for MPI_Reduce

double  MPI_Reduce_Opt(int argc, char *argv[], int myrank, int size)
{
  //int myrank, size;
  int count = atoi(argv[1]);

  double sTime, eTime, time;
  int color = myrank%2;
  int newrank, newsize;
  //int newrank2, newsize2;
  int offset = size/2;

  double sendval[count], maxval[count], maxval2[count];
  MPI_Status status;

  MPI_Comm newcomm;

   //splitted the communication world based on colors

  MPI_Comm_split (MPI_COMM_WORLD, color, myrank, &newcomm);




  MPI_Comm_rank( newcomm, &newrank );
  MPI_Comm_size( newcomm, &newsize );

 
  for(int i =0; i<count; i++)
	 
  sendval[i] = (double)rand() / (double)RAND_MAX;
  sTime = MPI_Wtime();
  
  //Performing reduce in each subcommunicator
  
  if(color==0)
  MPI_Reduce(sendval, maxval,count, MPI_DOUBLE, MPI_SUM, 0, newcomm); // find sum
  else
  MPI_Reduce(sendval, maxval,count, MPI_DOUBLE, MPI_SUM, 0, newcomm); // find sum



  //then sending the reduced value from one subcommunicator to other, so that finally we can obtained a single reduced value

  if(myrank==1)
  MPI_Send (maxval,count, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  if(myrank==0){
  MPI_Recv (maxval2,count, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
  for(int i=0; i<count; i++)
  {
	  maxval[i] = maxval[i]+maxval2[i];
  }
  
  }
  eTime = MPI_Wtime();
  time = eTime - sTime;
  MPI_Barrier (MPI_COMM_WORLD);
  //printf ("%lf\n", time);

  //MPI_Finalize();
  return time;

}




//optmised code for MPI_Gather

double MPI_Gather_Opt(int argc, char *argv[], int myrank, int size)
{
  //int myrank, size;
  int arrSize = atoi(argv[1]);
  double sTime, eTime, time;
  double sendvalue[arrSize], maxvalue[arrSize];
;
  MPI_Status status;


  int color = myrank%2;
  int newrank, newsize;
  int offset = size/2;
  MPI_Comm newcomm;
  
  
  
  //splitted the communication world based on colors
  
  MPI_Comm_split (MPI_COMM_WORLD, color, myrank, &newcomm);

  MPI_Comm_rank( newcomm, &newrank );
  MPI_Comm_size( newcomm, &newsize );

  double *message = (double *) malloc(arrSize * sizeof(double));
  double *recvMessage = (double *) malloc(size * arrSize * sizeof(double));
  //double finalarray[arrSize*size];
  // initialize array
  //srand(time(NULL));
  for (int i = 0; i < arrSize; i++) {
    message[i] = (double)rand() / (double)RAND_MAX;
  }
  sTime = MPI_Wtime();
  
  //performing gather on both subcommunicators
  
  if(color == 0)
   MPI_Gather(message, arrSize, MPI_DOUBLE, recvMessage, arrSize, MPI_DOUBLE, 0, newcomm);
  else
   MPI_Gather(message, arrSize, MPI_DOUBLE, recvMessage, arrSize, MPI_DOUBLE, 0, newcomm);




  //then sending the gathered values from one sub communicator to the other to know obtain the final gathered value

  if(myrank == 1)
	  MPI_Send (recvMessage,(arrSize*offset), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

  if(myrank == 0)
  {
	   MPI_Recv(recvMessage + (offset*arrSize),offset*arrSize, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
	
  }
  eTime = MPI_Wtime();
  time = eTime - sTime;

  //printf ("%lf\n", time);
MPI_Barrier (MPI_COMM_WORLD);

   //MPI_Finalize();
   return time;

}



//Optimised version of MPI_Alltoallv

double MPI_Alltoallv_Opt(int argc, char *argv[],int rank, int numtasks)
{
 
  double sTime,eTime,time;
 

  // Allocate message
  int totalSize = atoi(argv[1]);
  int arrSize = (totalSize/numtasks);
  double *message = (double *) malloc(arrSize * sizeof(double));
  double *newMessage = (double *) malloc(arrSize * sizeof(double));

  // initialize array
  //srand(time(NULL));
  for (int i = 0; i < arrSize; i++) {
    message[i] = (double)rand() / (double)RAND_MAX;
  }


  sTime = MPI_Wtime();
  
  // scatter to all processes
  MPI_Scatter (message, 1, MPI_DOUBLE, newMessage, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  // gather at root
  double *recvMessage = (double *) malloc(totalSize * sizeof(double)); //significant at the root process


  //running MPI GATHER at all other processes.
  for(int i=0;i<numtasks;i++)
  MPI_Gather (message, arrSize, MPI_DOUBLE, recvMessage, arrSize, MPI_DOUBLE, i, MPI_COMM_WORLD);
  eTime = MPI_Wtime();
  time = eTime - sTime;
  MPI_Barrier (MPI_COMM_WORLD);


  return time;
}
void avgtime(double time, int size, int data, int mode, int OperationNo)
{


          //Declared all files to take the outputs

	FILE *fp1, *fp2, *fp3, *fp4;
	fp1 = fopen("output1.txt", "a");
	fp2 = fopen("output2.txt", "a");
	fp3 = fopen("output3.txt", "a");
	fp4 = fopen("output4.txt", "a");
	if(size == 4)
         {
	       if(OperationNo == 1)
	       fprintf(fp1, "%d 4 1 %d %lf\n",data,mode,time);
	       else if(OperationNo == 2)
	       fprintf(fp2, "%d 4 1 %d %lf\n",data,mode,time);
	       else if(OperationNo == 3)
	       fprintf(fp3, "%d 4 1 %d %lf\n",data,mode,time);
	       else 
               fprintf(fp4, "%d 4 1 %d %lf\n",data,mode,time);

	       //printf("%d 4 1 %d %lf\n", data, mode, time);
         }
         else if(size == 32)
         {
	       if(OperationNo == 1)
               fprintf(fp1, "%d 4 8 %d %lf\n",data,mode,time);
               else if(OperationNo == 2)
               fprintf(fp2, "%d 4 8 %d %lf\n",data,mode,time);
               else if(OperationNo == 3)
               fprintf(fp3, "%d 4 8 %d %lf\n",data,mode,time);
               else 
               fprintf(fp4, "%d 4 8 %d %lf\n",data,mode,time);
	       //printf("%d 4 8 %d %lf\n", data, mode, time);
         }
         else if(size == 16)
         {
	       if(OperationNo == 1)
               fprintf(fp1, "%d 16 1 %d %lf\n",data,mode,time);
               else if(OperationNo == 2)
               fprintf(fp2, "%d 16 1 %d %lf\n",data,mode,time);
               else if(OperationNo == 3)
               fprintf(fp3, "%d 16 1 %d %lf\n",data,mode,time);
               else
               fprintf(fp4, "%d 16 1 %d %lf\n",data,mode,time);

	       //printf("%d 16 1 %d %lf\n", data, mode, time);
         }
         else if(size == 128)
         {
	       if(OperationNo == 1)
               fprintf(fp1, "%d 16 8 %d %lf\n",data,mode,time);
               else if(OperationNo == 2)
               fprintf(fp2, "%d 16 8 %d %lf\n",data,mode,time);
               else if(OperationNo == 3)
               fprintf(fp3, "%d 16 8 %d %lf\n",data,mode,time);
               else
               fprintf(fp4, "%d 16 8 %d %lf\n",data,mode,time);
	       //printf("%d 16 8 %d %lf\n", data, mode, time);
         }

}

int main( int argc, char *argv[] )
{
 
   int OperationNo = atoi(argv[2]);
   int count = atoi(argv[1]);
   if(argc < 3)
   {
       printf("too less arguements");
       return 0;
   }

   if(OperationNo == 1)
   {
	 double timeDefault = 0,timeOpt = 0; 
	 int numtasks,mode = 0;
	 int D = atoi(argv[1]);
	 int myrank,size;
	 MPI_Init( &argc, &argv );
         MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
         MPI_Comm_size( MPI_COMM_WORLD, &size );

	 timeDefault = MPI_Bcast_Default(argc,argv,myrank,size);
	 if(myrank == 0)
	 avgtime(timeDefault,size,D,mode,OperationNo);
	 MPI_Barrier (MPI_COMM_WORLD); 
	 
	 mode = 1;
	 timeOpt = MPI_Bcast_Opt(argc,argv,myrank,size); 
	 if(myrank == 0)
         avgtime(timeOpt,size,D,mode,OperationNo);
         
         MPI_Finalize();
   }
  else if(OperationNo == 2)
   {
	 double timeDefault = 0,timeOpt = 0;
         int mode = 0;
         int D = atoi(argv[1]);
         int myrank,size;
         MPI_Init( &argc, &argv );
         MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
         MPI_Comm_size( MPI_COMM_WORLD, &size );
         timeDefault = MPI_Reduce_Default(argc,argv,myrank,size);
	 if(myrank == 0)
         avgtime(timeDefault,size,D,mode,OperationNo);
         MPI_Barrier (MPI_COMM_WORLD); 
	 
	 mode = 1;
         timeOpt = MPI_Reduce_Opt(argc,argv,myrank,size);
	 if(myrank == 0)
         avgtime(timeOpt,size,D,mode,OperationNo);
         
         MPI_Finalize();
   }
   else if(OperationNo == 3)
   {
	double timeDefault = 0,timeOpt = 0;
        int mode = 0;
        int D = atoi(argv[1]);
	int myrank,size;
        MPI_Init( &argc, &argv );
        MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
        MPI_Comm_size( MPI_COMM_WORLD, &size );

        timeDefault = MPI_Gather_Default(argc,argv,myrank,size);
	
	if(myrank == 0)
        avgtime(timeDefault,size,D,mode,OperationNo);
        MPI_Barrier (MPI_COMM_WORLD);

        mode = 1;
        timeOpt = MPI_Gather_Opt(argc,argv,myrank,size);
	if(myrank == 0)
        avgtime(timeOpt,size,D,mode,OperationNo);
        
        MPI_Finalize();
   }
   if(OperationNo == 4)
   {
	double timeDefault = 0,timeOpt = 0;
        int mode = 0;
        int D = atoi(argv[1]);

	int myrank,size;
	MPI_Init( &argc, &argv );
        MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
        MPI_Comm_size( MPI_COMM_WORLD, &size);

        timeDefault = MPI_Alltoallv_Default(argc, argv, myrank, size);
	
	if(myrank == 0)
        avgtime(timeDefault,size,D,mode,OperationNo);
	MPI_Barrier (MPI_COMM_WORLD);
	
	mode = 1;

        timeOpt = MPI_Alltoallv_Opt(argc, argv, myrank, size);
	//MPI_Barrier (MPI_COMM_WORLD);
	if(myrank == 0)
        avgtime(timeOpt,size,D,mode,OperationNo);
        MPI_Finalize();
   }
}

