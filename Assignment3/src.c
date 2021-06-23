#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"



int main( int argc, char *argv[])
{


FILE *fp; //fp=file poiner
	 MPI_Init(&argc, &argv);
	
	 int myrank, size; //size will take care of number of processes 
		 
	 double sTime, eTime, time,max_Time;
	
	
	//declaring the array so that all process can access it	
 	float *arr2 = NULL;  
 	float **mat= NULL;

	  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
	  MPI_Comm_size(MPI_COMM_WORLD, &size);
	
		long long int r=0, c;  
	


if(size==1)
{

if(myrank==0)
{

     //counting the number of rows and columns
        
	FILE *mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	char line[20000];
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token;
		
		c=0;
		
		token=strtok(line, ",");
		
		while(token!=NULL)
		{
			
			c++;
			token= strtok(NULL, ",");
		}
		
		r++;
		
	}
	
	fclose(mf);
	
	r=r-1, c=c-2;
	//printf("The number of row and column is is %d %d \n", r, c);
	
	

       //declaring the matrix of size m[r][c]
	
	mat = (float **)malloc(r * sizeof(float *));
    	
        for (int i=0; i<r; i++)
        mat[i] = (float *)malloc(c * sizeof(float));
	
	
	
	
	// Now storing the temperatures only in matrix
	
	mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	
	int i=0,j=0;
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token1;
		
		j=0;
		
		token1=strtok(line, ",");
		
		while(token1!=NULL)
		{
			if(i>0 && j>1)
		        {
				float d;
	   			sscanf(token1, "%f", &d);
				//printf(" %0.2f",d);
				mat[i-1][j-2] = d;
			}
			
			token1= strtok(NULL, ",");
			
			j++;
		}
		
		i++;
		//printf("\n");
	}
	
	fclose(mf);
	
	
        //starting time
	sTime = MPI_Wtime(); 
         
         //To find the minimum of every year, we must at first declare an array having the size as total number of columns
	
	 arr2 = (float*)malloc(c * sizeof(float));
	 
	 float min;
	 
	 for (long long int j = 0; j <  c; j++)
	 {
	  	min=1000;
      		for (long long int i = 0; i < r; i++)
      		{
         		if(mat[i][j]<min)
         		min=mat[i][j];
         	}
         	
         	arr2[j]=min;
          }

	
	//Now we need to find the global minimum
	
	min=1000;
	
	for(long long int i=0;i<c;i++)
	{
		if(i==0)
		{
			fp=fopen("output.txt", "w");
	 		fprintf(fp, "%0.2f,",arr2[i]);
	 		fclose(fp);
	 	}
	 	else if(i==c-1)
	 	{
	 		fp=fopen("output.txt", "a");
	 		fprintf(fp, "%0.2f",arr2[i]);
	 		fclose(fp);
	 	}
	 	else
	 	{
	 		fp=fopen("output.txt", "a");
	 		fprintf(fp, "%0.2f,",arr2[i]);
	 		fclose(fp);
	 	}
	 	
	 	
	 	
		if(arr2[i]<min)
		min=arr2[i];
	}
	
	fp=fopen("output.txt", "a");
	fprintf(fp, "\n");
	fclose(fp);
	
	fp=fopen("output.txt", "a");
	fprintf(fp, "%0.2f\n",min);
	fclose(fp);
	
  eTime = MPI_Wtime();
  time = eTime - sTime;

  // obtain max time
  MPI_Reduce (&time, &max_Time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  if (myrank==0) 
  {
  	fp=fopen("output.txt", "a");
  	fprintf(fp, "%lf\n", max_Time);
  	fclose(fp);
  }	

   } //end of myrank==0
	


}




//When number of process is 2

if(size==2)
{
	if(myrank==0)
{
        //counting the number of rows and columns
        
	FILE *mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	char line[20000];
	
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token;
		
		c=0;
		
		token=strtok(line, ",");
		
		while(token!=NULL)
		{
			
			c++;
			token= strtok(NULL, ",");
		}
		
		r++;
		
	}
	
	fclose(mf);
	
	r=r-1, c=c-2;
	//printf("The number of row and column is is %d %d \n", r, c);
	
	

       //declaring the matrix of size m[r][c]
	
	mat = (float **)malloc(r * sizeof(float *));
    	
        for (int i=0; i<r; i++)
        mat[i] = (float *)malloc(c * sizeof(float));
	
	
	
	
	// Now storing the temperatures only in matrix
	
	mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	
	int i=0,j=0;
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token1;
		
		j=0;
		
		token1=strtok(line, ",");
		
		while(token1!=NULL)
		{
			if(i>0 && j>1)
		        {
				float d;
	   			sscanf(token1, "%f", &d);
				//printf(" %0.2f",d);
				mat[i-1][j-2] = d;
			}
			
			token1= strtok(NULL, ",");
			
			j++;
		}
		
		i++;
		//printf("\n");
	}
	
	fclose(mf);
}

         MPI_Barrier(MPI_COMM_WORLD);
	//starting time
	sTime = MPI_Wtime();

if(myrank==0)
{

	
	float **snd_buf = (float **)malloc( (c-c/2) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/2); i++)
        snd_buf[i] = (float *)malloc(r * sizeof(float));
        


	 arr2= (float*)malloc(c * sizeof(float));
	
	//total size of the data that is being sent
	long long int outs=r*4;
	
	 int position=0;
	
	
	//since we will be sending maximum c/2 number of columns so we have taken c
	MPI_Request request[c];

	
			//we will send the ha0.2f of the matrix to process 2
			for(long long int j=c/2;j<c;j++)
			{
				position=0; //reassigning position after each and every send
				
				for(long long int i=0;i<r;i++)
				{
					MPI_Pack(&mat[i][j], 1 , MPI_FLOAT,snd_buf[j-c/2],outs,&position,MPI_COMM_WORLD);
				}
			}
			
			
			//sending all the buffers
			for(long long int j=c/2;j<c;j++)
			{
				MPI_Isend (snd_buf[j-c/2], outs , MPI_PACKED, 1 /*dest*/ , j /*tag*/ , MPI_COMM_WORLD,&request[j]);
			}
			
	
	

	
	
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=0;j<c/2;j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(mat[i][j]<min1)
				 	min1=mat[i][j];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
			
	} //end of if (myrank==0) 	
	
	
	
	
	
	
	
	//Broadcasting the value of r and c to other process from the root process after reading the data
	 MPI_Bcast(&r, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	 MPI_Bcast(&c, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	 
	
	//total size of the data that is beign sent
	long long int outs=r*8;
	
	 int position=0;
	
	
	
	MPI_Request request1[c];
  	MPI_Status status[c];
  	
  	
		if(myrank!=0)
		{ 
		 	arr2= (float*)malloc(c * sizeof(float));
		 }
	
	
	float **recv_buf = (float **)malloc( (c-c/2) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/2); i++)
        recv_buf[i] = (float *)malloc(r * sizeof(float));
	
	
	
	float **buf = (float **)malloc( (c-c/2) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/2); i++)
        buf[i] = (float *)malloc(r * sizeof(float));

	
	
	//receiving the data
	
	if(myrank==1)
		{
		
			 for(long long int j=c/2;j<c;j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/2], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=c/2; j<c;j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/2],outs,&position,&buf[j-c/2][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=c/2;j<c;j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/2][i]<min1)
				 	min1=buf[j-c/2][i];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		}
		
		
		
		



		//now we will send this array recieve which storred the minimum from process2 to process1 using MPI_contiguous 
		
		MPI_Datatype type;
  		MPI_Type_contiguous( (c-c/2) , MPI_FLOAT, &type );
  		MPI_Type_commit(&type);
		
		if(myrank==1)
		{
			
	
				MPI_Send (&arr2[c/2], 1 , type, 0 /*dest*/ , 100 /*tag*/ , MPI_COMM_WORLD);
		
		}
		
		//receiving the data
	
		if(myrank==0)
		{
				
				 MPI_Recv(&arr2[c/2], 1, type, 1 /*src*/ , 100 /*tag*/, MPI_COMM_WORLD,&status[0]);
				
			
		} 	

		
		
		//printing the entire array
		
		//now we will print the minimum array using rank0 i.e the root
		
		if(myrank==0)
		{
		float min2=10000;
		int count=0; //to keeo track of number of elements in the array
		
			for(long long int i=0;i<c;i++)                        //printing the ha0.2f that is filled by myrank1
			{
				if(i==0)
			{
				fp=fopen("output.txt", "w");
		 		fprintf(fp, "%0.2f,",arr2[i]);
		 		fclose(fp);
		 	}
		 	else if(i==c-1)
		 	{
		 		fp=fopen("output.txt", "a");
		 		fprintf(fp, "%0.2f",arr2[i]);
		 		fclose(fp);
		 	}
		 	else
		 	{
		 		fp=fopen("output.txt", "a");
		 		fprintf(fp, "%0.2f,",arr2[i]);
		 		fclose(fp);
		 	}
				
				
				if(arr2[i]<min2)
				min2=arr2[i];
			}
			
			
		fp=fopen("output.txt", "a");
		fprintf(fp,"\n");
		fclose(fp);
		
		
		fp=fopen("output.txt", "a");
		fprintf(fp,"%0.2f\n",min2);
		fclose(fp);
		
			
		}
		
		
  eTime = MPI_Wtime();
  time = eTime - sTime;

  // obtain max time
  MPI_Reduce (&time, &max_Time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank==0) 
  {	
  	fp=fopen("output.txt", "a");
	fprintf(fp,"%lf\n", max_Time);
	fclose(fp);
  }
		



} //end of if(size==2)



if(size==4)
{
	if(myrank==0)
{

     //counting the number of rows and columns
        
	FILE *mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	char line[20000];
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token;
		
		c=0;
		
		token=strtok(line, ",");
		
		while(token!=NULL)
		{
			
			c++;
			token= strtok(NULL, ",");
		}
		
		r++;
		
	}
	
	fclose(mf);
	
	r=r-1, c=c-2;
	//printf("The number of row and column is is %d %d \n", r, c);
	
	

       //declaring the matrix of size m[r][c]
	
	mat = (float **)malloc(r * sizeof(float *));
    	
        for (int i=0; i<r; i++)
        mat[i] = (float *)malloc(c * sizeof(float));
	
	
	
	
	//// Now storing the temperatures only in matrix
	
	mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	
	int i=0,j=0;
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token1;
		
		j=0;
		
		token1=strtok(line, ",");
		
		while(token1!=NULL)
		{
			if(i>0 && j>1)
		        {
				float d;
	   			sscanf(token1, "%f", &d);
				//printf(" %0.2f",d);
				mat[i-1][j-2] = d;
			}
			
			token1= strtok(NULL, ",");
			
			j++;
		}
		
		i++;
		//printf("\n");
	}
	
	fclose(mf);
	
}	
 
     MPI_Barrier(MPI_COMM_WORLD); 
     //starting time
	sTime = MPI_Wtime();
         
 if(myrank==0)
 {        
         //To find the minimum of every year, we must at first declare an array having the size as total number of columns
	
	 float* arr;
	 arr = (float*)malloc(c * sizeof(float));

	
	float **snd_buf = (float **)malloc( (c-c/4) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/4); i++)
        snd_buf[i] = (float *)malloc(r * sizeof(float));
        
        


	 arr2= (float*)malloc(c * sizeof(float));
	
	//total size of the data that is beign sent
	long long int outs=r*4;
	
	 int position=0;
	
	
	//since we will be sending maximum (c-c/4) number of columns so we have taken c as maximum
	MPI_Request request[c];
  	
	
	
	
			//we are packing the remaining portions after c/4
			for(long long int j=c/4;j<c;j++)
			{
				position=0; //reassigning position after each and every send
				
				for(long long int i=0;i<r;i++)
				{
					MPI_Pack(&mat[i][j], 1 , MPI_FLOAT,snd_buf[j-c/4],outs,&position,MPI_COMM_WORLD);
				}
			}
			
			
			
			//sending to P1 and P2 from P0
			for(int i=1;i<3;i++)
			{
				
				for(long long int j=(i*c/4); j<( (i*c/4) + c/4) ; j++)
				{
					MPI_Isend (snd_buf[j-c/4], outs , MPI_PACKED, i /*dest*/ , j /*tag*/ , MPI_COMM_WORLD,&request[j]);
				}
			}
			
			
			
			//sending to P3
			for(long long int j=(3*(c/4)); j<c ; j++)
			{
				MPI_Isend (snd_buf[j-c/4], outs , MPI_PACKED, 3 /*dest*/ , j /*tag*/ , MPI_COMM_WORLD,&request[j]);
			}
			
			
			//min1 will store the minimum for rank0
			float min1=1000;
			
			for(long long int j=0;j<c/4;j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(mat[i][j]<min1)
				 	min1=mat[i][j];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		
			
			
			
	

} //end of myrank==0







	//Broadcasting the value of r and c to other process from the root process after reading the data
	 MPI_Bcast(&r, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	 MPI_Bcast(&c, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	 
	 

        float **recv_buf = (float **)malloc( (c-c/4) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/4); i++)
        recv_buf[i] = (float *)malloc(r * sizeof(float));
	
	
	
	float **buf = (float **)malloc( (c-c/4) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/4); i++)
        buf[i] = (float *)malloc(r * sizeof(float));
        
        
     
		
		
		if(myrank!=0)
		{ 
		 	arr2= (float*)malloc(c * sizeof(float));
		 }
	
        
        MPI_Request request1[c];
  	MPI_Request request2[c];
  	MPI_Request request3[c];
  	MPI_Status status[c];
  	MPI_Status status2[c];

	
	//total size of the data that is beign sent
	long long int outs=r*4;
	
	 int position=0;
		
		
	
	
	//receiving the data at P1
	
	if(myrank==1)
		{
		
			 for(long long int j=c/4;j<(c/4 + c/4);j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/4], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=c/4;j<(c/4 + c/4);j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/4],outs,&position,&buf[j-c/4][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
			//min1 will store the minimum for rank1
			float min1=1000;
			
		
			
			for(long long int j=c/4;j<(c/4 + c/4); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/4][i]<min1)
				 	min1=buf[j-c/4][i];
				 }
				 
				 arr2[j]=min1;
				// printf(" %0.2f", arr2[j]);
				 
				 
			}
			
		}
		
		
		
		//receiving the data at P2
	
	if(myrank==2)
		{
		
			 for(long long int j=2*c/4; j<( (2*c/4) + c/4); j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/4], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=2*c/4;j<( (2*c/4) + c/4); j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/4],outs,&position,&buf[j-c/4][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
			
		
		
		
			
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=2*c/4; j<( (2*c/4) + c/4); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/4][i]<min1)
				 	min1=buf[j-c/4][i];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		}
		
		
		
		
		
		
		
		//receiving the data at P3
	
	if(myrank==3)
		{
		
			 for(long long int j=3*c/4; j<c; j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/4], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=3*c/4; j<c; j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/4],outs,&position,&buf[j-c/4][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
			
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=3*c/4; j<c; j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/4][i]<min1)
				 	min1=buf[j-c/4][i];
				 }
				 
				 arr2[j]=min1;
				// printf("%0.2f ", arr2[j]);
				 
			}
			
		}
		
		
		
		




		//now we will send this array recieve which storred the minimum from process2 to process1 using MPI_contiguous 
		
		
		//for process 1 and 2, size will be c/4
		MPI_Datatype type;
  		MPI_Type_contiguous( c/4 , MPI_FLOAT, &type );
  		MPI_Type_commit(&type);
  		
  		
  		
  		//for process 1 and 2, size will be c/4, for process 3 it will be c- 3*(c/4) since the last process also has the remainder rows
  		MPI_Datatype type1;
  		MPI_Type_contiguous( c-(3*c/4) , MPI_FLOAT, &type1 );
  		MPI_Type_commit(&type1);
  		  
  		
  		
		
		
		//now sending the data from P1
		if(myrank==1)
		{
			
				
				MPI_Isend (&arr2[c/4], 1 , type, 0 /*dest*/ , 100 /*tag*/ , MPI_COMM_WORLD ,&request2[1]);
		
		}
		
		
		//now sending the data from P2
		if(myrank==2)
		{
			
	
				MPI_Isend (&arr2[2*c/4], 1 , type, 0 /*dest*/ , 200 /*tag*/ , MPI_COMM_WORLD, &request2[2]);
		
		}
		
		
		//now sending the data from P3
		if(myrank==3)
		{
			
	
				MPI_Isend (&arr2[3*c/4], 1 , type1, 0 /*dest*/ , 300 /*tag*/ , MPI_COMM_WORLD, &request2[3]);
		
		}
		
		
		
		//receiving the data at rank 0
	
		if(myrank==0)
		{
				 
				
				 MPI_Irecv(&arr2[1*(c/4)], 1, type, 1 /*src*/ , 100 /*tag*/, MPI_COMM_WORLD,&request3[1]);
				 MPI_Wait(&request3[1], &status2[1]);
				 
				 MPI_Irecv(&arr2[2*(c/4)], 1, type, 2 /*src*/ , 200 /*tag*/, MPI_COMM_WORLD,&request3[2]);
				 MPI_Wait(&request3[2], &status2[2]);
				 
				 MPI_Irecv(&arr2[3*(c/4)], 1, type1, 3 /*src*/ , 300 /*tag*/, MPI_COMM_WORLD,&request3[3]);
				 MPI_Wait(&request3[3], &status2[3]);
				
			
		} 	

		
		
		//printing the entire array
		
		//now we will print the minimum array using rank0 i.e the root
		
		if(myrank==0)
		{
		float min2=10000;
		int count=0; //to keep track of number of elements in the array
		
			for(long long int i=0;i<c;i++)                        //printing the ha0.2f that is filled by myrank1
			{
				if(i==0)
				{
					fp=fopen("output.txt", "w");
			 		fprintf(fp, "%0.2f,",arr2[i]);
			 		fclose(fp);
			 	}
			 	else if(i==c-1)
			 	{
			 		fp=fopen("output.txt", "a");
			 		fprintf(fp, "%0.2f",arr2[i]);
			 		fclose(fp);
			 	}
			 	else
			 	{
			 		fp=fopen("output.txt", "a");
			 		fprintf(fp, "%0.2f,",arr2[i]);
			 		fclose(fp);
			 	}
				
				
				if(arr2[i]<min2)
				min2=arr2[i];
			}
			
			
		fp=fopen("output.txt", "a");
		fprintf(fp,"\n");
		fclose(fp);
		
		
		fp=fopen("output.txt", "a");
		fprintf(fp,"%0.2f\n",min2);
		fclose(fp);
		
			
		}
		
		
  eTime = MPI_Wtime();
  time = eTime - sTime;

  // obtain max time
  MPI_Reduce (&time, &max_Time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (myrank==0) 
  {	
  	fp=fopen("output.txt", "a");
	fprintf(fp,"%lf\n", max_Time);
	fclose(fp);
  }
		



} //end of if(size==4)
	
	
	
	

//when number of process is 8	
	
if(size==8)
{	
	

	  
if(myrank==0)
{	  

     //counting the number of rows and columns
        
	FILE *mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	char line[20000];
	
	
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token;
		
		c=0;
		
		token=strtok(line, ",");
		
		while(token!=NULL)
		{
			
			c++;
			token= strtok(NULL, ",");
		}
		
		r++;
		
	}
	
	fclose(mf);
	
	r=r-1, c=c-2;
	//printf("The number of row and column is is %d %d \n", r, c);
	
	
	

       //declaring the matrix of size m[r][c]
	
	mat = (float **)malloc(r * sizeof(float *));
    	
        for (int i=0; i<r; i++)
        mat[i] = (float *)malloc(c * sizeof(float));
	
	
	
	
	// // Now storing the temperatures only in matrix
	
	mf=fopen(argv[1], "r");
	if(mf==NULL)
	{
		perror("Unable to open the file");
		exit(1);
	}
	
	
	int i=0,j=0;
	
	while(fgets(line,sizeof(line), mf))
	{
		char *token1;
		
		j=0;
		
		token1=strtok(line, ",");
		
		while(token1!=NULL)
		{
			if(i>0 && j>1)
		        {
				float d;
	   			sscanf(token1, "%f", &d);
				//printf(" %0.2f",d);
				mat[i-1][j-2] = d;
			}
			
			token1= strtok(NULL, ",");
			
			j++;
		}
		
		i++;
		//printf("\n");
	}
	
	fclose(mf);
	
}


 MPI_Barrier(MPI_COMM_WORLD);
//starting time
sTime = MPI_Wtime();


if(myrank==0)
{
	
		
	float **snd_buf = (float **)malloc( (c-c/8) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/8); i++)
        snd_buf[i] = (float *)malloc(r * sizeof(float));
        


	
	
	 	arr2= (float*)malloc(c * sizeof(float));
	
	//total size of the data that is beign sent
	long long int outs=r*4;
	
	 int position=0;
	
	
	//since we will be sending maximum c/2 number of columns so we have taken c
	MPI_Request request[c];
  
	
	
	//packing and sending the data
	
			
			for(long long int j=c/8;j<c;j++)
			{
				position=0; //reassigning position after each and every send
				
				for(long long int i=0;i<r;i++)
				{
					MPI_Pack(&mat[i][j], 1 , MPI_FLOAT,snd_buf[j-c/8],outs,&position,MPI_COMM_WORLD);
				}
			}
			
			
			
			//sending to P1,P2,P3,P4,P5,P6 from P0
			for(int i=1;i<7;i++)
			{
				
				for(long long int j=(i*c/8); j<( (i*c/8) + c/8) ; j++)
				{
					MPI_Isend (snd_buf[j-c/8], outs , MPI_PACKED, i /*dest*/ , j /*tag*/ , MPI_COMM_WORLD,&request[j]);
				}
			}
			
			
			
			//sending to P7
			for(long long int j=(7*(c/8)); j<c ; j++)
			{
				MPI_Isend (snd_buf[j-c/8], outs , MPI_PACKED, 7 /*dest*/ , j /*tag*/ , MPI_COMM_WORLD,&request[j]);
			}
			
			
			//min1 will store the minimum for rank0
			float min1=1000;
			
			for(long long int j=0;j<c/8;j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(mat[i][j]<min1)
				 	min1=mat[i][j];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
			
			
	
			
} 	//end of myrank==0	
	
	
	
	
	
	

	//Broadcasting the value of r and c to other process from the root process after reading the data
	 MPI_Bcast(&r, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	 MPI_Bcast(&c, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	
	
	float **recv_buf = (float **)malloc( (c-c/8) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/8); i++)
        recv_buf[i] = (float *)malloc(r * sizeof(float));
	
	
	
	float **buf = (float **)malloc( (c-c/8) * sizeof(float *));
    	
        for (long long int i=0; i<(c-c/8); i++)
        buf[i] = (float *)malloc(r * sizeof(float));


	
	
		
		
		if(myrank!=0)
		{ 
		 	arr2= (float*)malloc(c * sizeof(float));
		 }
	
		
	//total size of the data that is beign sent
	long long int outs=r*4;
	
	 int position=0;
		
	MPI_Request request1[c];
  	MPI_Request request2[c];
  	MPI_Request request3[c];
  	MPI_Status status[c];
  	MPI_Status status2[c];
	
	//receiving the data at P1
	
	if(myrank==1)
		{
		
			 for(long long int j=c/8;j<(c/8 + c/8);j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/8], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=c/8;j<(c/8 + c/8);j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/8],outs,&position,&buf[j-c/8][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
			//min1 will store the minimum for rank1
			float min1=1000;
			
		
			
			for(long long int j=c/8;j<(c/8 + c/8); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/8][i]<min1)
				 	min1=buf[j-c/8][i];
				 }
				 
				 arr2[j]=min1;
				// printf(" %0.2f", arr2[j]);
				 
				 
			}
			
		}
		
		
		
		//receiving the data at P2
	
	if(myrank==2)
		{
		
			 for(long long int j=2*c/8; j<( (2*c/8) + c/8); j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/8], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=2*c/8;j<( (2*c/8) + c/8); j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/8],outs,&position,&buf[j-c/8][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
		
		
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=2*c/8; j<( (2*c/8) + c/8); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/8][i]<min1)
				 	min1=buf[j-c/8][i];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		}
		
		
		
		
		
		//receiving the data at P3
	
	if(myrank==3)
		{
		
			 for(long long int j=3*c/8; j<( (3*c/8) + c/8); j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/8], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=3*c/8;j<( (3*c/8) + c/8); j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/8],outs,&position,&buf[j-c/8][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
		
		
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=3*c/8; j<( (3*c/8) + c/8); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/8][i]<min1)
				 	min1=buf[j-c/8][i];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		}
		
		
		
		
		
		
		//receiving the data at P4
	
	if(myrank==4)
		{
		
			 for(long long int j=4*c/8; j<( (4*c/8) + c/8); j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/8], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=4*c/8;j<( (4*c/8) + c/8); j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/8],outs,&position,&buf[j-c/8][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
		
		
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=4*c/8; j<( (4*c/8) + c/8); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/8][i]<min1)
				 	min1=buf[j-c/8][i];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		}
		
		
		
		
		//receiving the data at P5
	
	if(myrank==5)
		{
		
			 for(long long int j=5*c/8; j<( (5*c/8) + c/8); j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/8], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=5*c/8;j<( (5*c/8) + c/8); j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/8],outs,&position,&buf[j-c/8][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
		
		
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=5*c/8; j<( (5*c/8) + c/8); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/8][i]<min1)
				 	min1=buf[j-c/8][i];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		}
		
		
		
		//receiving the data at P6
	
	if(myrank==6)
		{
		
			 for(long long int j=6*c/8; j<( (6*c/8) + c/8); j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/8], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=6*c/8;j<( (6*c/8) + c/8); j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/8],outs,&position,&buf[j-c/8][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
		
		
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=6*c/8; j<( (6*c/8) + c/8); j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/8][i]<min1)
				 	min1=buf[j-c/8][i];
				 }
				 
				 arr2[j]=min1;
				 
				 
			}
			
		}
		
		
		
		
		
		
		
		
		
		
		//receiving the data at P7
	
	if(myrank==7)
		{
		
			 for(long long int j=7*c/8; j<c; j++)
			 {		
				 MPI_Irecv(recv_buf[j-c/8], outs, MPI_PACKED, 0 /*src*/ , j /*tag*/, MPI_COMM_WORLD, &request1[j]);
				 MPI_Wait(&request1[j], &status[j]);
			 }
			
			for(long long int j=7*c/8; j<c; j++)
			{
				position=0;
				 for(long long int i=0;i<r;i++)
				 {
				 	MPI_Unpack(recv_buf[j-c/8],outs,&position,&buf[j-c/8][i], 1/*outcount*/, MPI_FLOAT, MPI_COMM_WORLD);
				 }
			}
			
			
			
		
			//min1 will store the minimum for rank1
			float min1=1000;
			
			for(long long int j=7*c/8; j<c; j++)
			{
				min1=1000;
				
				for(long long int i=0;i<r;i++)
				 {
				 	if(buf[j-c/8][i]<min1)
				 	min1=buf[j-c/8][i];
				 }
				 
				 arr2[j]=min1;
				// printf("%0.2f ", arr2[j]);
				 
			}
			
		}
		
		
		
		




		//now we will send this array recieve which storred the minimum from process2 to process1 using MPI_contiguous 
		
		
		//for process 1 to 6, size will be c/8
		MPI_Datatype type;
  		MPI_Type_contiguous( c/8 , MPI_FLOAT, &type );
  		MPI_Type_commit(&type);
  		
  		
  		
  		//for process 1 to 6, size will be c/8, for process 7 it will be c- 7*(c/4) since the last process also has the remainder rows
  		MPI_Datatype type1;
  		MPI_Type_contiguous( c-(7*c/8) , MPI_FLOAT, &type1 );
  		MPI_Type_commit(&type1);
  		  
  		
  		
		
		
		//now sending the data from P1
		if(myrank==1)
		{
			
				
				MPI_Isend (&arr2[c/8], 1 , type, 0 /*dest*/ , 100 /*tag*/ , MPI_COMM_WORLD ,&request2[1]);
		
		}
		
		
		//now sending the data from P1
		if(myrank==2)
		{
			
	
				MPI_Isend (&arr2[2*c/8], 1 , type, 0 /*dest*/ , 200 /*tag*/ , MPI_COMM_WORLD, &request2[2]);
		
		}
		
		
		
		//now sending the data from P1
		if(myrank==3)
		{
			
	
				MPI_Isend (&arr2[3*c/8], 1 , type, 0 /*dest*/ , 300 /*tag*/ , MPI_COMM_WORLD, &request2[3]);
		
		}
		
		
		//now sending the data from P1
		if(myrank==4)
		{
			
	
				MPI_Isend (&arr2[4*c/8], 1 , type, 0 /*dest*/ , 400 /*tag*/ , MPI_COMM_WORLD, &request2[4]);
		
		}
		
		
		
		//now sending the data from P1
		if(myrank==5)
		{
			
	
				MPI_Isend (&arr2[5*c/8], 1 , type, 0 /*dest*/ , 500 /*tag*/ , MPI_COMM_WORLD, &request2[5]);
		
		}
		
		
		
		//now sending the data from P1
		if(myrank==6)
		{
			
	
				MPI_Isend (&arr2[6*c/8], 1 , type, 0 /*dest*/ , 600 /*tag*/ , MPI_COMM_WORLD, &request2[6]);
		
		}
		
		
		//now sending the data from P1
		if(myrank==7)
		{
			
	
				MPI_Isend (&arr2[7*c/8], 1 , type1, 0 /*dest*/ , 700 /*tag*/ , MPI_COMM_WORLD, &request2[7]);
		
		}
		
		
		
		
		
		//receiving the data
	
		if(myrank==0)
		{
				 
				
				 MPI_Irecv(&arr2[1*(c/8)], 1, type, 1 /*src*/ , 100 /*tag*/, MPI_COMM_WORLD,&request3[1]);
				 MPI_Wait(&request3[1], &status2[1]);
				 
				 MPI_Irecv(&arr2[2*(c/8)], 1, type, 2 /*src*/ , 200 /*tag*/, MPI_COMM_WORLD,&request3[2]);
				 MPI_Wait(&request3[2], &status2[2]);
				 
				 
				 MPI_Irecv(&arr2[3*(c/8)], 1, type, 3 /*src*/ , 300 /*tag*/, MPI_COMM_WORLD,&request3[3]);
				 MPI_Wait(&request3[3], &status2[3]);
				 
				 
				 
				 MPI_Irecv(&arr2[4*(c/8)], 1, type, 4 /*src*/ , 400 /*tag*/, MPI_COMM_WORLD,&request3[4]);
				 MPI_Wait(&request3[4], &status2[4]);
				 
				 
				 MPI_Irecv(&arr2[5*(c/8)], 1, type, 5 /*src*/ , 500 /*tag*/, MPI_COMM_WORLD,&request3[5]);
				 MPI_Wait(&request3[5], &status2[5]);
				 
				 MPI_Irecv(&arr2[6*(c/8)], 1, type, 6 /*src*/ , 600 /*tag*/, MPI_COMM_WORLD,&request3[6]);
				 MPI_Wait(&request3[6], &status2[6]);
				 
				 
				 MPI_Irecv(&arr2[7*(c/8)], 1, type1, 7 /*src*/ , 700 /*tag*/, MPI_COMM_WORLD,&request3[7]);
				 MPI_Wait(&request3[7], &status2[7]);
				
			
		} 	

		
		
		//printing the entire array
		
		//now we will print the minimum array using rank0 i.e the root
		
			if(myrank==0)
		{
		float min2=10000;
		int count=0; //to keeo track of number of elements in the array
		
			for(long long int i=0;i<c;i++)                        //printing the ha0.2f that is filled by myrank1
			{
				if(i==0)
				{
					fp=fopen("output.txt", "w");
			 		fprintf(fp, "%0.2f,",arr2[i]);
			 		fclose(fp);
			 	}
			 	else if(i==c-1)
			 	{
			 		fp=fopen("output.txt", "a");
			 		fprintf(fp, "%0.2f",arr2[i]);
			 		fclose(fp);
			 	}
			 	else
			 	{
			 		fp=fopen("output.txt", "a");
			 		fprintf(fp, "%0.2f,",arr2[i]);
			 		fclose(fp);
			 	}
				
				if(arr2[i]<min2)
				min2=arr2[i];
			}
			
			
		fp=fopen("output.txt", "a");
		fprintf(fp,"\n");
		fclose(fp);
		
		
		fp=fopen("output.txt", "a");
		fprintf(fp,"%0.2f\n",min2);
		fclose(fp);
		
			
		}
		
		
  eTime = MPI_Wtime();
  time = eTime - sTime;
  // obtain max time
  MPI_Reduce (&time, &max_Time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
 // printf("%d %lf\n",myrank,max_Time);
  if (myrank==0) 
  {	
  	fp=fopen("output.txt", "a");
	fprintf(fp,"%lf\n", max_Time);
	fclose(fp);
  }
		



} //end of if(size==8)

		
	    MPI_Finalize();
  	    return 0;
	
	
}
	

