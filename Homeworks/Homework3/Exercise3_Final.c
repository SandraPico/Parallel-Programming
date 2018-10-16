#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

//For the mpi library
MPI_Status status;

//Compare function for the quicksort algorithm.
int compareFunction (const void * a, const void * b){
    if( (*(double*)a > *(double*)b )){
        return 1;
    }else if(*(double*)a < *(double*)b ){
        return -1;
    }else{
        return 0;
    }
}

//Print the result through the command line.
void PrintResult(double* data, int N){
    int i = 0;
    for (i=0;i<N;i++){
        printf("%f\n",data[i]); 
    }
}

//Merge the two arrays
void merge_Function(double * subArray, int slice, double*neighbour_subArray, int neighbour_slice, int type ){

    double *mergeList = (double*)malloc((neighbour_slice+slice)*sizeof(double));
    int total_slice = neighbour_slice + slice;
    int index_1 = 0;
    int index_2 = 0;

    for (int i = 0; i < total_slice; i++ ){
        if (index_2 < neighbour_slice && index_1 < slice){
            if (neighbour_subArray[index_2] > subArray[index_1]){
                mergeList[i] = subArray[index_1];
                index_1  = index_1 +1;
            }else{
                mergeList[i] = neighbour_subArray[index_2];
                index_2 = index_2 +1;
            }
        }else if (!(index_1 < slice)){
            mergeList[i] = neighbour_subArray[index_2];
            index_2 = index_2+1;
        }else{
            mergeList[i] = subArray[index_1];
            index_1 = index_1 +1;
        }
    }
    //Low-high merge
    if (type == 0){
        int j = 0;
        for (j = 0; j < slice; j++){
            subArray[j] = mergeList[j+neighbour_slice];
        }
    }else{
        int x = 0;
        for(x= 0; x< slice; x++){
            subArray[x] = mergeList[x];
        }
    }     
}

//Function to create the file with the results in order to be able to plot them. 
void CreateFile(double *data, int N){
    FILE * f;
    f = fopen("Result.txt","a");
    for (int i = 0; i < N; i ++){
        fprintf(f, "%f/", data[i]);
    }
    fclose(f);
}

//Exchange the proper information taking into account the even-odd phase 
void Exchange_Information(int type,int P, int N, int rank,double* subArray, int slice){
    if (type == 1){
        int neighbour_slice = (N+P-rank+1-1)/P;
        double *neighbour_subArray = (double *) malloc(neighbour_slice*sizeof(double));
        MPI_Send(subArray, slice, MPI_DOUBLE, rank+1, 10, MPI_COMM_WORLD);
        MPI_Recv(&neighbour_subArray[0], neighbour_slice, MPI_DOUBLE, rank+1, 10, MPI_COMM_WORLD, &status);
        merge_Function(subArray,slice,neighbour_subArray, neighbour_slice, 1);
    }else{
        int neighbour_slice = (N+P-rank+1-1)/P;
        double *neighbour_subArray = (double *) malloc(neighbour_slice*sizeof(double));
        MPI_Recv(&neighbour_subArray[0], neighbour_slice, MPI_DOUBLE, rank-1, 10, MPI_COMM_WORLD, &status);
        MPI_Send(subArray, slice, MPI_DOUBLE, rank-1, 10, MPI_COMM_WORLD);
        merge_Function(subArray,slice,neighbour_subArray, neighbour_slice, 0);
    }
}

//Even Odd Transposation algorithm.
void EvenOdd(int N, double *data,int size, int rank, double*subArray, int P, int slice, int Print){
    int even_process = (rank%2 == 0);
    int even_phase = 1;
    int i = 0;
    //Time needed to be sure that everything is completely sorted.
    for( i = 0 ; i < P ; i++ ) {
        if ((even_phase && (rank%2 == 0)) || (!even_phase && !(rank%2 ==0))){
            if (rank < (P-1)){
                Exchange_Information(1,P,N,rank,subArray,slice);
            }
        }
        else if (rank > 0){
            Exchange_Information(0,P,N,rank,subArray,slice);
        }
        even_phase = !even_phase;
    }
    //Concatenate all the arrays into one (data)
    MPI_Gather(subArray,slice,MPI_DOUBLE,data,slice,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if ( rank == 0){
        CreateFile(data,N);
        if (Print == 1){
            PrintResult(data,N);
        }
    }
}

int main(int argc, char **argv) {

    //Initialize to print the cpu time.
    //3.2 ask for ploting the efficiency
    clock_t start, end;
    double cpu_time_used;

    int Print = 0;

    //Start the clock
    start = clock();
    
    //Read the N number as an argument
    int N = atoi(argv[1]);
    
    //Variables for the MPI interface
    int tag = 10;
    int size = 0;
    int rank = 0;
    int rc = 0;
    //MPI Initialization
    rc = MPI_Init(NULL,NULL);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //Total data
    double *data = (double*)malloc(N*sizeof(double));
    //In order to print the cpu time needed.
    double *cpu_array = (double*)malloc(size*sizeof(double));
    // The size variable define the number of processes. 
    int P = size;  

    //Calculate the slide regarding each particular case (for each process)
    int slice = (N-rank-1+P)/P;

    double *subArray = (double*)malloc(slice*sizeof(double));
    //Fill each subArray(per each processor) with random values. (between 0-1)
    srand(rank);
    int x = 0;
    for (x = 0; x < slice; x++){
        subArray[x] = (double)rand()/RAND_MAX;
    }

    //Local sort
    qsort(subArray,slice,sizeof(double),compareFunction);

    //Here, all the processes have their own data already sorted.
    //If we only have 1 processor, then everything is done though the quicksort algorithm. 
    if (P == 1){
        CreateFile(subArray,N);
        MPI_Finalize();
    }else{
        //Start OddEven algorithm to sort
        EvenOdd(N,data,size,rank,subArray,P,slice,Print);
        MPI_Finalize();
    }
    end = clock();
    //Calculate the maximum cpu time.
    cpu_array[rank] = ((double) (end - start)) / CLOCKS_PER_SEC;
    double max_cpu = 0;
    if ( rank == 0 ){
        for (int i = 0; i < P ; i ++){
            if (cpu_array[i] > max_cpu){
                max_cpu = cpu_array[i];
            }
        }
        //Print the cpu time needed to finish the program.
        printf("The cpu time used is: %f\n",max_cpu);
    }
    return 0;
}
