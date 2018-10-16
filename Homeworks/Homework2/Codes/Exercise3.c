#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "math.h"
#include <mpi.h>
#include <string.h>

MPI_Status status;

#define N 1000
#define iterations 800000
#define h 1.0/(N+1)


//r(x) = (-(10(x-0.5)^2))
double r_function(double x){
    double result = -(10*(pow((x-0.5),2)));
    return result;
}

//Just to compare with the exact result (generate error)
//u(x) = (40(x-0.5)^2-10)
double u_function(double x){
    double result = (40*(pow((x-0.5),2))- 10);
    return result;
}


//f(x) = 80 - 10(40(x - 0.5)^2-10)(x - 0.5)^2
double f_function(double x){
    double result = 80 - 10*(40*(pow((x - 0.5),2))-10) *(pow((x - 0.5),2));
    return result;
}


// Main Function 
int main (int argc, char **argv) {

	//To store the errors per iteration and then be able to plot it.
    double error_final[iterations];

    //Variables needed for the MPI implementation
    int tag = 10;
    int size = 0;
    int rank = 0;
    int rc = 0;

    //MPI Initialization
    rc = MPI_Init(NULL,NULL);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //The size variable define the number of processes. 
    int P = size;

    //Linear data distribution with overlap (Lecture slides)
    int a = 2;
    int M = N + (P-1)*a;
    int L = M/P;
    int R = M % P;
    int Ip;

    //Compute Ip ( How many numbers/elements will have each processor.)
    //Lecture slides
    if ( rank < R){
        Ip = L+1;
    }else{
        Ip = L;
    }

    double x[Ip];
    //Init value for each process. Depending of the p value.
    int int_value = (L-a)*rank;
    int aux = 0;
    if (R >= rank){
        aux = rank;
    }else{
        aux = R;
    }
    int_value = int_value + aux;

    //Compute the x_n values.
    for (int i = 0; i<Ip; i++ ){
        x[i] = h*(int_value+i+1);
    }

    //Compute the f(xn) i r(xn) per each value.
    double fx[Ip];
    double rx[Ip];
    double u_exact[Ip];

    //Compute the three functions needed per all values of x (per each process)
    for (int i = 0; i < Ip; i++){
        fx[i] = f_function(x[i]);
        rx[i] = r_function(x[i]);
        u_exact[i] = u_function(x[i]);
    }

    double u[Ip];
    double u_error[Ip];
    double new_u[Ip];

    //Init
    for (int i = 0; i < Ip; i++){
        u[i] = 0;
        u_error[i] = 0;
        new_u[i] = 0;
    }

    //Going through all the iterations.
    //The convergence is slowly so you need to do a lot of iterations.
    for (int iter = 0; iter <iterations;iter ++){
                
        //Red-black communication to fill the ghost cells. 
        if((rank % 2 ) == 0){
        	//To fill the ghost cells you will need to identify the first process and the last one 
        	//(only send/fill the specific ghost cells per each case)
            //Even
            if (rank != (size -1)){
                MPI_Send(&u[Ip-2],1, MPI_DOUBLE, rank+1,tag, MPI_COMM_WORLD);
                MPI_Recv(&u[Ip-1],1, MPI_DOUBLE, rank+1,tag, MPI_COMM_WORLD,&status);
            }
            if (rank != 0){
                MPI_Send(&u[1],1, MPI_DOUBLE, rank-1,tag, MPI_COMM_WORLD);
                MPI_Recv(&u[0],1, MPI_DOUBLE, rank-1,tag, MPI_COMM_WORLD,&status);
            }
        }else{
            //Odd
            if (rank != 0){
                MPI_Recv(&u[0],1, MPI_DOUBLE, rank-1,tag, MPI_COMM_WORLD,&status);
                MPI_Send(&u[1],1, MPI_DOUBLE, rank-1,tag, MPI_COMM_WORLD);
            }
            if (rank != (size-1)){
                MPI_Recv(&u[Ip-1],1, MPI_DOUBLE, rank+1,tag, MPI_COMM_WORLD,&status);
                MPI_Send(&u[Ip-2],1, MPI_DOUBLE, rank+1,tag, MPI_COMM_WORLD);
            }
        }
   
   		//Compute u(x) but not filling the ghost cells.
        int init = 1;
        for (int r = init; r <= Ip -2 ; r++){
            //The sequence u[k] is given by...
            new_u[r] = (u[r-1]+u[r+1]-(pow(h,2))*(fx[r]))/(2.0-((pow(h,2))*rx[r]));
        }

        //Update u
        for (int i = 0; i < Ip; i++){
            u[i] = new_u[i];
        }
      
        //Compute error comparing with the exact function.
        for (int i = 0; i < Ip; i++){
            u_error[i] = pow(u[i] - u_exact[i],2);
        }

        //Compute the sum of the u_error:
        double sum_error = 0;
        double globalsum = 0;
        for ( int i = 0; i < Ip; i++){
            sum_error = sum_error + u_error[i];
        }

        //Now we need to sum all the errors per all the processes.
        MPI_Allreduce(&sum_error,&globalsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        //Now, global sum has all the errors regarding/iteration.
        //This is important to be able to plot the error.
        error_final[iter] = globalsum/N;
    }

    //Now we have all the errors. We are going to store the errors in a file and then plot them.
    FILE *error_file;
    error_file = fopen("error.txt","wb");
    for(int i = 0; i < iterations; i++){
        fprintf(error_file,"%.5f\n",error_final[i]);
    }
    fclose(error_file);

    //After all the iterations we are going to create a file per each process.
    //The file will store the u values, and then we are going to print them through an external python code.
    FILE *fptr;
    char name[50];
    //The name of the file will have the name of the process.
    sprintf(name,"%d",rank);
    char direction[4] = ".txt";
    strcat(name,direction);
    printf("Name: %s\n",name);
    fptr = fopen(name,"wb");
   	
   	//Store the u(x) per different x cases.
    if (rank == 0){
        for(int i = 0; i<(Ip-1); i++){
            fprintf(fptr,"%.5f",u[i]);
            fprintf(fptr,",%.5f\n",x[i]);
        }
    }else{
        //Take into account the ghost cells
        if(rank == size-1){
            for(int i = 1; i<Ip; i++){
                fprintf(fptr,"%.5f",u[i]);
                fprintf(fptr,",%.5f\n",x[i]);
            }
        }else{
            for (int i = 1; i< Ip-1; i++){
                fprintf(fptr,"%.5f",u[i]);
                fprintf(fptr,",%.5f\n",x[i]);
            }
        }
    }
    //Close the file
    fclose(fptr);
    //Close the MPI library/communication.
    MPI_Finalize();
}
