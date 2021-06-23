// Timing codes

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include "mpi.h"

int main( int argc, char *argv[]) //the first argument takes the value of number of data points, and the second one num_time_steps 
{
FILE *fp; //fp=file poiner

for(int m=1;m<=3;m++)
{

if(m==1)  //using MPI_Isend()
{

		MPI_Init(&argc, &argv);

		int t=5;


		while(t)
		{  
		 
		 
		  int myrank, size; //size will take care of number of processes 
		 
		  double sTime, eTime, time,max_time;

		  

		  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
		  MPI_Comm_size(MPI_COMM_WORLD, &size);

		  int count = atoi (argv[1]);
		      
		  //count=count/size;  //Since number of data points per process is provided, so no need to do this

		  count= sqrt(count); //count is actually the dimension of the matrix

		   //original matrix
		   //double m[count][count];

		  
		  //the matrix in which the updated values will be stored
		  //double up_m[count][count];
		  
		  
		  //declaring the original matrix and the updated matrix dynamically
		  
		    double **m = (double **)malloc(count * sizeof(double *)); 
		    for (int i=0; i<count; i++)
		    { 
			 m[i] = (double *)malloc(count * sizeof(double)); 
		    }
	
	
	  	    double **up_m = (double **)malloc(count * sizeof(double *)); 
		    for (int i=0; i<count; i++)
		    { 
			 up_m[i] = (double *)malloc(count * sizeof(double)); 
		    }
		  
		  
		  
		  //declaring the 4 buffers for a process for recieving from the neighbouring process
		  double buf1[count];
		  double buf2[count];
		  double buf3[count];
		  double buf4[count];
		  
		  
		//setting the value of all the buffers to 0
		  
		for(int i=0;i<count;i++)
		{
			buf1[i]=0;
			buf2[i]=0;
			buf3[i]=0;
			buf4[i]=0;
		}
		  
		  
		  MPI_Request request[size];
		  MPI_Request request1[size];
		  MPI_Status status[size];

		  
		//filling the matrix by randomly
		  
		for(int i=0;i<count;i++)
		{
			for(int j=0;j<count;j++)
			{
				m[i][j]=rand()%100;  				
			}

		}


		

		//the distance 

		int d=sqrt(size);
		int num_time_steps=atoi(argv[2]);

		int s=sqrt(size);   //s is the dimension of the matrix of processes
		
		/* Say we have 9 processes then s=3
		
		 p0 p1 p2
		 p3 p4 p5
		 p6 p7 p8
		 
		 */
		  


		//for loop should start from here, for 50 iterations

		 sTime = MPI_Wtime();
		 
		 for(int iter=0;iter<num_time_steps;iter++)
		 {
		 

			//Now we will have to deal with send

			//Now what we will do is that, we will check whether there exists a process whose rank is within the limits and if so then we will forward our 	              elements to them

			if(myrank%s!=0)
			{
				//we will send the left most column of this matrix to the buffer3 of that process
				for(int i=0;i<count;i++)
				{
					for(int j=0;j<1;j++)
					{
						MPI_Isend (&m[i][j], 1 , MPI_DOUBLE, myrank-1, i, MPI_COMM_WORLD,&request[myrank-1]);
					}
				}
				
			} 	
			  

			//we will send the right most column to the buffer1 of that process
			
			if( (myrank+1)%s !=0)
			{

				for(int i=0;i<count;i++)
				{
					for(int j=count-1;j<count;j++)
					{
						MPI_Isend (&m[i][j], 1 , MPI_DOUBLE,myrank+1, i, MPI_COMM_WORLD,&request[myrank+1]);
					}
				}
				

			}


			if(myrank+d<size)
			{
				//we will send the last row to the buffer2 of that process
				for(int i=count-1;i<count;i++)
				{
					for(int j=0;j<count;j++)
					{
						MPI_Isend (&m[i][j], 1 , MPI_DOUBLE,myrank+d, j, MPI_COMM_WORLD,&request[myrank+d]);
					}
				}

			}


			if(myrank-d>=0)
			{
				//we will send the last row to the buffer2 of that process
				for(int i=0;i<1;i++)
				{
					for(int j=0;j<count;j++)
					{
						MPI_Isend (&m[i][j], 1 , MPI_DOUBLE,myrank-d, j, MPI_COMM_WORLD,&request[myrank-d]);
					}
				}

			}



			 


			///Now we will have to deal with recieve

			if(myrank%s!=0)
		{
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Irecv(&buf1[i], 1, MPI_DOUBLE, myrank-1, i, MPI_COMM_WORLD, &request1[myrank-1]);
			 	MPI_Wait(&request1[myrank-1], &status[myrank-1]);
			 }
		}

		if(myrank-d>=0)
		{
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Irecv(&buf2[i], 1, MPI_DOUBLE, myrank-d, i, MPI_COMM_WORLD, &request1[myrank-d]);
			 	MPI_Wait(&request1[myrank-d], &status[myrank-d]);
			 }
		}

		if((myrank+1)%s!=0)
		{
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Irecv(&buf3[i], 1, MPI_DOUBLE, myrank+1, i, MPI_COMM_WORLD, &request1[myrank+1]);
			 	MPI_Wait(&request1[myrank+1], &status[myrank+1]);
			 }
		}


		if(myrank+d<size)
		{
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Irecv(&buf4[i], 1, MPI_DOUBLE, myrank+d, i, MPI_COMM_WORLD, &request1[myrank+d]);
			 	 MPI_Wait(&request1[myrank+d], &status[myrank+d]);
			 }
		}


			//Now our recieve is complete from all the neighbours, we can insert MPI_Barrier here if we want

			

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////
			////
			////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			

			//Now we must calculate the new matrix
			 
			  //lets first calculate for the center elements who are not in the halo region, i.e. they are getting data from the same process itself
			  
			  for(int i=1;i<count-1;i++)
			  {
				for(int j=1;j<count-1;j++)
				{
					up_m[i][j]=( m[i][j-1]+m[i][j+1]+m[i-1][j]+m[i+1][j] )/4;
				}
			  }
			 
			 
			 //Now lets calculate for each corner
			 
			 
			 //top left
			 
			 if(buf1[0]==0 && buf2[0]==0)
			 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/2;
			 else if(buf1[0]==0 || buf2[0]==0)
			 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/3;
			 else
			 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/4;
			 
			 
			 
			 //top right
			 
			 if(buf2[count-1]==0 && buf3[0]==0)
			 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/2;
			 else if(buf2[count-1]==0 || buf3[0]==0)
			 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/3;
			 else
			 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/4;
			 
			 
			 
			 
			 //bottom left
			 
			 if(buf1[count-1]==0 && buf4[0]==0)
			 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/2;
			 else if(buf1[count-1]==0 || buf4[0]==0)
			 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/3;
			 else
			 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/4;
			 
			 
			 //bottom right
			 if(buf3[count-1]==0 && buf4[count-1]==0)
			 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/2;
			 else if(buf3[count-1]==0 || buf4[count-1]==0)
			 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/3;
			 else 
			 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/4;




			//Now lets calculate for the remaining HALO regions except the corner cases


			//top row

			for(int j=1;j<count-1;j++)
			{
				if(buf2[j]==0)
				up_m[0][j]=( m[0][j-1] + buf2[j] + m[0][j+1] +m[1][j] )/3;
				else
				up_m[0][j]=( m[0][j-1] + buf2[j] + m[0][j+1] +m[1][j] )/4;
			}

			//bottom row

			for(int j=1;j<count-1;j++)
			{
				if(buf4[j]==0)
				up_m[count-1][j]=( m[count-1][j-1]+m[count-2][j] + m[count-1][j+1]+buf4[j]  )/3;
				else
				up_m[count-1][j]=( m[count-1][j-1]+m[count-2][j] + m[count-1][j+1]+buf4[j]  )/4;
			}

			//leftmost column

			for(int i=1;i<count-1;i++)
			{
				if(buf1[i]==0)
				up_m[i][0]=( buf1[i]+m[i-1][0]+ m[i][1]+ m[i+1][0] )/3;
				else
				up_m[i][0]=( buf1[i]+m[i-1][0]+ m[i][1]+ m[i+1][0] )/4;	
			}


			//rightmost column

			for(int i=1;i<count-1;i++)
			{	
				if(buf3[i]==0)
				up_m[i][count-1]=( m[i][count-2]+m[i-1][count-1]+buf3[i]+m[i+1][count-1] )/3;
				else
				up_m[i][count-1]=( m[i][count-2]+m[i-1][count-1]+buf3[i]+m[i+1][count-1] )/4;
			}


			//So till this stage we are done with all the calculations, now we need to copy this to our data to original matrix





			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////
			////
			////
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			for(int i=0;i<count;i++)
			{
				for(int j=0;j<count;j++)
				{
					m[i][j]=up_m[i][j];
				}

			}
			 
			 
			 
		} //end of iter for loop		 
		 
		 	
		  eTime = MPI_Wtime();     
		  
		  time = eTime - sTime;

		 
		  MPI_Reduce (&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, size-1, MPI_COMM_WORLD);
		  
		  
		  //now we will create the logic of the file
		  
		  
		  
		  
		  //for P=16
		  if(size==16 && myrank == size-1)
		  {
		  
			  if(atoi(argv[1])==256 && t==5)
			  {
			        fp=fopen("Data_P16.txt", "w");
			  	fprintf(fp, "%lf\n", max_time); 
			  	fclose(fp);
			  }
			  
			  else
			  {		
			  	fp=fopen("Data_P16.txt", "a");
			  	if (myrank == size-1) fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  }
			 
			 // MPI_Barrier(MPI_COMM_WORLD); 
		  }	  
		  
		
		  
		   
		   
		   //for P=36
		  if(size==36 && myrank == size-1)
		  {
		  
			  if(atoi(argv[1])==256 && t==5)
			  {
			        fp=fopen("Data_P36.txt", "w");
			  	fprintf(fp, "%lf\n", max_time); 
			  	fclose(fp);
			  }
			  
			  else
			  {		
			  	fp=fopen("Data_P36.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  }
			  //MPI_Barrier(MPI_COMM_WORLD); 
		  }	  
		  
		  
		 
		//for P=49
		
		 
		  if(size==49 && myrank == size-1)
		  {
		  
			  if(atoi(argv[1])==256 && t==5)
			  {
			        fp=fopen("Data_P49.txt", "w");
			  	fprintf(fp, "%lf\n", max_time); 
			  	fclose(fp);
			  }
			  
			  else
			  {		
			  	fp=fopen("Data_P49.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  }
			  //MPI_Barrier(MPI_COMM_WORLD); 
		  }	  
		
		
		 
		  if(size==64 && myrank == size-1)
		  {
		  
			  if(atoi(argv[1])==256 && t==5)
			  {
			        fp=fopen("Data_P64.txt", "w");
			  	fprintf(fp, "%lf\n", max_time); 
			  	fclose(fp);
			  }
			  
			  else
			  {		
			  	fp=fopen("Data_P64.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  }
			  //MPI_Barrier(MPI_COMM_WORLD); 
		  }	  
		
		
		t--;

		}//end of while loop


//		MPI_Finalize();


}//end of if 1


if(m==2)  //using MPI_Pack()
{
	
//MPI_Init(&argc, &argv);

  int t=5;
  
  
  while(t--)
 {  
	 
	 
	  int myrank, size; //size will take care of number of processes 
	 
	  double sTime, eTime, time,max_time;

	  

	  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
	  MPI_Comm_size(MPI_COMM_WORLD, &size);

	  int count = atoi (argv[1]);
  	      
  	  //count=count/size;  //Since number of data points per process is provided, so no need to do this
    
	  count= sqrt(count); //count is actually the dimension of the matrix

	   //original matrix
	  //double m[count][count];

	  
	  //the matrix in which the updated values will be stored
	  //double up_m[count][count];
	   
	   
	   //declaring the original matrix and the updated matrix dynamically
		  
		    double **m = (double **)malloc(count * sizeof(double *)); 
		    for (int i=0; i<count; i++)
		    { 
			 m[i] = (double *)malloc(count * sizeof(double)); 
		    }
	
	
	  	    double **up_m = (double **)malloc(count * sizeof(double *)); 
		    for (int i=0; i<count; i++)
		    { 
			 up_m[i] = (double *)malloc(count * sizeof(double)); 
		    }
	   
	  
	  MPI_Request request[size];
  	  MPI_Request request1[size];
  	  MPI_Status status[size];

	  
	//filling the matrix by randomly
	  
	for(int i=0;i<count;i++)
	{
		for(int j=0;j<count;j++)
		{
			m[i][j]=rand()%100;   				
		}

	}


	 

	//the distance 

	int d=sqrt(size);
        int num_time_steps=atoi(argv[2]);
        
        int s=sqrt(size);

	
	//for loop should start from here, for 50 iterations
	
	 sTime = MPI_Wtime();
         
         for(int iter=0;iter<num_time_steps;iter++)
         {
         
        
          //declaring the 4 final buffers for a process for recieving from the neighbouring process in which the recv_buffer unpacks
	  
	  double buf1[count];
	  double buf2[count];
	  double buf3[count];
	  double buf4[count];
         
         
         //declaring the sending buffers
	  
	  double snd_buf1[count];
	  double snd_buf2[count];
	  double snd_buf3[count];
	  double snd_buf4[count];
	  
	  
	  //declaring the recieving buffers
	  
	  double recv_buf1[count];
	  double recv_buf2[count];
	  double recv_buf3[count];
	  double recv_buf4[count];
	  
	 
	  
	//setting the value of all the buffers to 0
	  
	for(int i=0;i<count;i++)
	{
		buf1[i]=0; buf2[i]=0; buf3[i]=0; buf4[i]=0;
		snd_buf1[i]=0; snd_buf2[i]=0; snd_buf3[i]=0; snd_buf4[i]=0; 
		recv_buf1[i]=0; recv_buf2[i]=0; recv_buf3[i]=0; recv_buf4[i]=0;
	}
	 
	 
	 
	  
	   int position1=0,position2=0,position3=0,position4=0;
         
         
         //this indicates the size of the buffer
         int outs=count*8;

		//Now we will have to deal with send

		//Now what we will do is that, we will check whether there exists a process whose rank is within the limits and if so then we will forward our 	              elements to them

		if(myrank%s!=0)
		{
			//we will send the left most column of this matrix using buf1 of this process to the buffer3 of that process
			for(int i=0;i<count;i++)
			{
				for(int j=0;j<1;j++)
				{
					MPI_Pack(&m[i][j], 1 , MPI_DOUBLE,snd_buf1,outs,&position1,MPI_COMM_WORLD);
				}
			}
			
			
			MPI_Isend (snd_buf1, count , MPI_PACKED, myrank-1, 1, MPI_COMM_WORLD,&request[myrank-1]);
		} 	
		  

		

		if(myrank-d>=0)
		{
			//we will send the last row to the buffer2 of that process
			for(int i=0;i<1;i++)
			{
				for(int j=0;j<count;j++)
				{
					MPI_Pack(&m[i][j], 1 , MPI_DOUBLE,snd_buf2,outs,&position2,MPI_COMM_WORLD);
				}
			}

			MPI_Isend (snd_buf2, count , MPI_PACKED, myrank-d, 2, MPI_COMM_WORLD,&request[myrank-d]);

		}

				
		


                //we will send the right most column to the buffer1 of that process
                
		if( (myrank+1)%s !=0)
		{
	
			for(int i=0;i<count;i++)
			{
				for(int j=count-1;j<count;j++)
				{
					MPI_Pack(&m[i][j], 1 , MPI_DOUBLE,snd_buf3,outs,&position3,MPI_COMM_WORLD);
					
				}
			}
			
			MPI_Isend (snd_buf3, count , MPI_PACKED, myrank+1, 3, MPI_COMM_WORLD,&request[myrank+1]);
		}


		if(myrank+d<size)
		{
			//we will send the last row to the buffer2 of that process
			for(int i=count-1;i<count;i++)
			{
				for(int j=0;j<count;j++)
				{
					MPI_Pack(&m[i][j], 1 , MPI_DOUBLE,snd_buf4,outs,&position4,MPI_COMM_WORLD);
				}
			}
		
			MPI_Isend (snd_buf4, count , MPI_PACKED, myrank+d, 4, MPI_COMM_WORLD,&request[myrank+d]);
		}


		


		

		
		//re-initialising positions for using in recieve
		position1=0; position2=0; position3=0; position4=0;

		//Now we will have to deal with recieve

		if(myrank%s!=0)
		{
		
			 MPI_Irecv(recv_buf1, count, MPI_PACKED, myrank-1, 3, MPI_COMM_WORLD, &request1[myrank-1]);
			 MPI_Wait(&request1[myrank-1], &status[myrank-1]);
			
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Unpack(recv_buf1,outs,&position1,&buf1[i],1, MPI_DOUBLE, MPI_COMM_WORLD);
			 }
		}

		if(myrank-d>=0)
		{
			 MPI_Irecv(recv_buf2, count, MPI_PACKED, myrank-d, 4, MPI_COMM_WORLD, &request1[myrank-d]);
			 MPI_Wait(&request1[myrank-d], &status[myrank-d]);
			
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Unpack(recv_buf2,outs,&position2,&buf2[i],1, MPI_DOUBLE, MPI_COMM_WORLD);
			 }
		}

		if((myrank+1)%s!=0)
		{
			 MPI_Irecv(recv_buf3, count, MPI_PACKED, myrank+1, 1, MPI_COMM_WORLD, &request1[myrank+1]);
			 MPI_Wait(&request1[myrank+1], &status[myrank+1]);
			 
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Unpack(recv_buf3,outs,&position3,&buf3[i],1, MPI_DOUBLE, MPI_COMM_WORLD);
			 }
		}


		if(myrank+d<size)
		{
		
			 MPI_Irecv(recv_buf4, count, MPI_PACKED, myrank+d, 2, MPI_COMM_WORLD, &request1[myrank+d]);
			 MPI_Wait(&request1[myrank+d], &status[myrank+d]);
			 
			 for(int i=0;i<count;i++)
			 {
			 	MPI_Unpack(recv_buf4,outs,&position4,&buf4[i],1, MPI_DOUBLE, MPI_COMM_WORLD);
			 }
		}



		//Now our recieve is complete from all the neighbours, we can insert MPI_Barrier here if we want

		

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////
		////
		////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		

		//Now we must calculate the new matrix
		 
		  //lets first calculate for the center elements who are not in the halo region, i.e. they are getting data from the same process itself
		  
		  for(int i=1;i<count-1;i++)
		  {
			for(int j=1;j<count-1;j++)
			{
				up_m[i][j]=( m[i][j-1]+m[i][j+1]+m[i-1][j]+m[i+1][j] )/4;
			}
		  }
		 
		 
		 //Now lets calculate for each corner
		 
		 
		 //top left
		 
		 if(buf1[0]==0 && buf2[0]==0)
		 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/2;
		 else if(buf1[0]==0 || buf2[0]==0)
		 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/3;
		 else
		 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/4;
		 
		 
		 
		 //top right
		 
		 if(buf2[count-1]==0 && buf3[0]==0)
		 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/2;
		 else if(buf2[count-1]==0 || buf3[0]==0)
		 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/3;
		 else
		 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/4;
		 
		 
		 
		 
		 //bottom left
		 
		 if(buf1[count-1]==0 && buf4[0]==0)
		 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/2;
		 else if(buf1[count-1]==0 || buf4[0]==0)
		 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/3;
		 else
		 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/4;
		 
		 
		 //bottom right
		 if(buf3[count-1]==0 && buf4[count-1]==0)
		 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/2;
		 else if(buf3[count-1]==0 || buf4[count-1]==0)
		 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/3;
		 else 
		 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/4;




		//Now lets calculate for the remaining HALO regions except the corner cases


		//top row

		for(int j=1;j<count-1;j++)
		{
			if(buf2[j]==0)
			up_m[0][j]=( m[0][j-1] + buf2[j] + m[0][j+1] +m[1][j] )/3;
			else
			up_m[0][j]=( m[0][j-1] + buf2[j] + m[0][j+1] +m[1][j] )/4;
		}

		//bottom row

		for(int j=1;j<count-1;j++)
		{
			if(buf4[j]==0)
			up_m[count-1][j]=( m[count-1][j-1]+m[count-2][j] + m[count-1][j+1]+buf4[j]  )/3;
			else
			up_m[count-1][j]=( m[count-1][j-1]+m[count-2][j] + m[count-1][j+1]+buf4[j]  )/4;
		}

		//leftmost column

		for(int i=1;i<count-1;i++)
		{
			if(buf1[i]==0)
			up_m[i][0]=( buf1[i]+m[i-1][0]+ m[i][1]+ m[i+1][0] )/3;
			else
			up_m[i][0]=( buf1[i]+m[i-1][0]+ m[i][1]+ m[i+1][0] )/4;	
		}


		//rightmost column

		for(int i=1;i<count-1;i++)
		{	
			if(buf3[i]==0)
			up_m[i][count-1]=( m[i][count-2]+m[i-1][count-1]+buf3[i]+m[i+1][count-1] )/3;
			else
			up_m[i][count-1]=( m[i][count-2]+m[i-1][count-1]+buf3[i]+m[i+1][count-1] )/4;
		}


		//So till this stage we are done with all the calculations, now we need to copy this to our data to original matrix





		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////
		////
		////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		for(int i=0;i<count;i++)
		{
			for(int j=0;j<count;j++)
			{
				m[i][j]=up_m[i][j];
			}

		}
		 
		 	 
	} //end of iter for loop		 
	 
	 	
	       
	  eTime = MPI_Wtime();
	  time = eTime - sTime;

	 
	  MPI_Reduce (&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, size-1, MPI_COMM_WORLD);
  	
  	 
  	  
  	   if(size==16 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P16.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  //MPI_Barrier(MPI_COMM_WORLD); 
		  }
  	
	  
  	   if(size==36 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P36.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  //MPI_Barrier(MPI_COMM_WORLD); 
		  } 
		  
		  
	   
  	   if(size==49 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P49.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
		//MPI_Barrier(MPI_COMM_WORLD); 	  
		  }
	

	  
  	   if(size==64 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P64.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  	//MPI_Barrier(MPI_COMM_WORLD); 
			  
		  }
		

 }//end of while loop
 
 
 // MPI_Finalize();
  
  
  
}//end if2


if(m==3) //using MPI_type_contiguous()
{
	

	 int t=5;
  
  
  while(t--)
 {  
	 
	 
	  int myrank, size; //size will take care of number of processes 
	 
	  double sTime, eTime, time,max_time;

	

	  MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;
	  MPI_Comm_size(MPI_COMM_WORLD, &size);

	  int count = atoi (argv[1]);
  	    
  	  //count=count/size;  //Since number of data points per process is provided, so no need to do this
    
	  count= sqrt(count); //count is actually the dimension of the matrix

	   //original matrix
	  // double m[count][count];

	  
	  //the matrix in which the updated values will be stored
	  // double up_m[count][count];
	    
	  
	  //declaring the original matrix and the updated matrix dynamically
		  
		    double **m = (double **)malloc(count * sizeof(double *)); 
		    for (int i=0; i<count; i++)
		    { 
			 m[i] = (double *)malloc(count * sizeof(double)); 
		    }
	
	
	  	    double **up_m = (double **)malloc(count * sizeof(double *)); 
		    for (int i=0; i<count; i++)
		    { 
			 up_m[i] = (double *)malloc(count * sizeof(double)); 
		    }
	  
	  
	  MPI_Request request[size];
  	  MPI_Request request1[size];
  	  MPI_Status status[size];

	  
	//filling the matrix by randomly
	  
	for(int i=0;i<count;i++)
	{
		for(int j=0;j<count;j++)
		{
			m[i][j]=rand()%100;  				
		}

	}


	

	//the distance 

	int d=sqrt(size);
        int num_time_steps=atoi(argv[2]);
        
        int s=sqrt(size);

	
	//for loop should start from here, for 50 iterations
	
	 sTime = MPI_Wtime();
         
         for(int iter=0;iter<num_time_steps;iter++)
         {
         
        
          //declaring the 4 final buffers for a process for recieving from the neighbouring process in which the recv_buffer unpacks
	  
	  double buf1[count];
	  double buf2[count];
	  double buf3[count];
	  double buf4[count];
         
         	  MPI_Datatype type;
         	  MPI_Datatype newtype;
         
         //declaring the sending buffers
	  
	  double snd_buf1[count];
	  double snd_buf2[count];
	  double snd_buf3[count];
	  double snd_buf4[count];
	  
	  
	 
	  
	//setting the value of all the buffers to 0
	  
	for(int i=0;i<count;i++)
	{
		buf1[i]=0; buf2[i]=0; buf3[i]=0; buf4[i]=0;
		snd_buf1[i]=0; snd_buf2[i]=0; snd_buf3[i]=0; snd_buf4[i]=0; 
		
	}
	 
	 		 

	 
	  
	   int position1=0,position2=0,position3=0,position4=0;
           MPI_Type_contiguous( count, MPI_DOUBLE, &type );
           MPI_Type_commit(&type);
           MPI_Type_vector(count,1,count, MPI_DOUBLE,&newtype);
           MPI_Type_commit(&newtype);
       

		//Now we will have to deal with send

		//Now what we will do is that, we will check whether there exists a process whose rank is within the limits and if so then we will forward our 	              elements to them

		if(myrank%s!=0)
		{
			//we will send the left most column of this matrix using buf1 of this process to the buffer3 of that process
		
			
			//MPI_Isend (&m[0][0], 1 , newtype, myrank-1, 1, MPI_COMM_WORLD,&request[myrank-1]);
			
			//we will send the left most column of this matrix using buf1 of this process to the buffer3 of that process
			for(int i=0;i<count;i++)
			{
				for(int j=0;j<1;j++)
				{
					//MPI_Pack(&m[i][j], 1 , MPI_DOUBLE,snd_buf1,outs,&position1,MPI_COMM_WORLD);
					
					snd_buf1[i]=m[i][j];
				}
			}
			
			
			MPI_Isend (snd_buf1, 1 , type, myrank-1, 1, MPI_COMM_WORLD,&request[myrank-1]);
		} 	
		  
		

		

		if(myrank-d>=0)
		{
			//we will send the first row to the buffer2 of that process
			
			
/*			for(int i=0;i<1;i++)
			{
				for(int j=0;j<count;j++)
				{
					snd_buf2[j]=m[i][j];
				}
			}
*/			
			MPI_Isend (&m[0][0], 1 , type, myrank-d, 2, MPI_COMM_WORLD,&request[myrank-d]);

		}

				
		


                //we will send the right most column to the buffer1 of that process
                
		if( (myrank+1)%s !=0)
		{
	
			//MPI_Isend (&m[0][count], 1 , newtype, myrank-1, 1, MPI_COMM_WORLD,&request[myrank-1]);
			
			for(int i=0;i<count;i++)
			{
				for(int j=count-1;j<count;j++)
				{
					snd_buf3[i]=m[i][j];
					
				}
			}
			
			MPI_Isend (snd_buf3, 1 , type, myrank+1, 3, MPI_COMM_WORLD,&request[myrank+1]);
			
			
		}

		
		

		if(myrank+d<size)
		{
		
		
			//we will send the last row to the buffer2 of that process
/*			for(int i=count-1;i<count;i++)
			{
				for(int j=0;j<count;j++)
				{
					snd_buf4[j]=m[i][j];
				}
			}
			
*/			//we will send the last row
			MPI_Isend (&m[count-1][0], 1 , type, myrank+d, 4, MPI_COMM_WORLD,&request[myrank+d]);
			
		}


		


		 

		
		

		//Now we will have to deal with recieve

		if(myrank%s!=0)
		{
		
			 MPI_Irecv(buf1, 1, type, myrank-1, 3, MPI_COMM_WORLD, &request1[myrank-1]);
			 MPI_Wait(&request1[myrank-1], &status[myrank-1]);
			 
		}

		
		if(myrank-d>=0)
		{
			 MPI_Irecv(buf2, 1, type, myrank-d, 4, MPI_COMM_WORLD, &request1[myrank-d]);
			MPI_Wait(&request1[myrank-d], &status[myrank-d]);
			 
		}
		 
	 

		if((myrank+1)%s!=0)
		{
			 MPI_Irecv(buf3, 1, type, myrank+1, 1, MPI_COMM_WORLD, &request1[myrank+1]);
			MPI_Wait(&request1[myrank+1], &status[myrank+1]);
			 
			
		}

		
		if(myrank+d<size)
		{
		
			 MPI_Irecv(buf4, 1, type, myrank+d, 2, MPI_COMM_WORLD, &request1[myrank+d]);
			 MPI_Wait(&request1[myrank+d], &status[myrank+d]);
			 
		}

		

		

		//Now our recieve is complete from all the neighbours, we can insert MPI_Barrier here if we want



		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////
		////
		////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		

		//Now we must calculate the new matrix
		 
		  //lets first calculate for the center elements who are not in the halo region, i.e. they are getting data from the same process itself
		  
		  for(int i=1;i<count-1;i++)
		  {
			for(int j=1;j<count-1;j++)
			{
				up_m[i][j]=( m[i][j-1]+m[i][j+1]+m[i-1][j]+m[i+1][j] )/4;
			}
		  }
		 
		 
		 //Now lets calculate for each corner
		 
		 
		 //top left
		 
		 if(buf1[0]==0 && buf2[0]==0)
		 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/2;
		 else if(buf1[0]==0 || buf2[0]==0)
		 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/3;
		 else
		 up_m[0][0]=(buf1[0]+buf2[0]+m[0][1]+m[1][0])/4;
		 
		 
		 
		 //top right
		 
		 if(buf2[count-1]==0 && buf3[0]==0)
		 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/2;
		 else if(buf2[count-1]==0 || buf3[0]==0)
		 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/3;
		 else
		 up_m[0][count-1]=( buf2[count-1]+buf3[0]+m[0][count-2]+m[1][count-1] )/4;
		 
		 
		 
		 
		 //bottom left
		 
		 if(buf1[count-1]==0 && buf4[0]==0)
		 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/2;
		 else if(buf1[count-1]==0 || buf4[0]==0)
		 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/3;
		 else
		 up_m[count-1][0]=( buf1[count-1]+buf4[0]+m[count-2][0]+m[count-1][1] )/4;
		 
		 
		 //bottom right
		 if(buf3[count-1]==0 && buf4[count-1]==0)
		 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/2;
		 else if(buf3[count-1]==0 || buf4[count-1]==0)
		 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/3;
		 else 
		 up_m[count-1][count-1]=( buf3[count-1]+buf4[count-1]+m[count-2][count-1]+m[count-1][count-2] )/4;




		//Now lets calculate for the remaining HALO regions except the corner cases


		//top row

		for(int j=1;j<count-1;j++)
		{
			if(buf2[j]==0)
			up_m[0][j]=( m[0][j-1] + buf2[j] + m[0][j+1] +m[1][j] )/3;
			else
			up_m[0][j]=( m[0][j-1] + buf2[j] + m[0][j+1] +m[1][j] )/4;
		}

		//bottom row

		for(int j=1;j<count-1;j++)
		{
			if(buf4[j]==0)
			up_m[count-1][j]=( m[count-1][j-1]+m[count-2][j] + m[count-1][j+1]+buf4[j]  )/3;
			else
			up_m[count-1][j]=( m[count-1][j-1]+m[count-2][j] + m[count-1][j+1]+buf4[j]  )/4;
		}

		//leftmost column

		for(int i=1;i<count-1;i++)
		{
			if(buf1[i]==0)
			up_m[i][0]=( buf1[i]+m[i-1][0]+ m[i][1]+ m[i+1][0] )/3;
			else
			up_m[i][0]=( buf1[i]+m[i-1][0]+ m[i][1]+ m[i+1][0] )/4;	
		}


		//rightmost column

		for(int i=1;i<count-1;i++)
		{	
			if(buf3[i]==0)
			up_m[i][count-1]=( m[i][count-2]+m[i-1][count-1]+buf3[i]+m[i+1][count-1] )/3;
			else
			up_m[i][count-1]=( m[i][count-2]+m[i-1][count-1]+buf3[i]+m[i+1][count-1] )/4;
		}


		//So till this stage we are done with all the calculations, now we need to copy this to our data to original matrix





		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////
		////
		////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		for(int i=0;i<count;i++)
		{
			for(int j=0;j<count;j++)
			{
				m[i][j]=up_m[i][j];
			}

		}
		 
		
		 
	} //end of iter for loop		 
	 
	 	
	       
	  eTime = MPI_Wtime();
	  time = eTime - sTime;

	 
	  MPI_Reduce (&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, size-1, MPI_COMM_WORLD);
  	  
  	  
  	 
  	  
  	     if(size==16 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P16.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  //MPI_Barrier(MPI_COMM_WORLD); 
		  }
	  
	  
	  
  	  
  	     if(size==36 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P36.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			  //MPI_Barrier(MPI_COMM_WORLD); 
		  }
		  
		  
	  
	 
  	  
  	     if(size==49 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P49.txt", "a");
			  	fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			 // 	MPI_Barrier(MPI_COMM_WORLD); 
			  
		  }
	  
	  
	 
  	  
  	     if(size==64 && myrank == size-1)
		  {
		  
			  	fp=fopen("Data_P64.txt", "a");
			        fprintf(fp, "%lf\n", max_time);
			  	fclose(fp);
			 // MPI_Barrier(MPI_COMM_WORLD); 
		  }
	  

 }//end of while loop
 
 
 
 



} //end if3




}//end of for loop

    MPI_Finalize();
  return 0;

}


