#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>

MPI_Status status;

// ----------------------------------
// MANDELBROT ALGORITHM PARAMETERS 
//
// width  = Width of area to calculate (pixels)
// height = Height of area to calculate (pixels)
// N      = Number of different colors
// b      = Bound value for Mandelbrot set. Modyfying this value acts as a zoom.
//
#define width  6000
#define height 6000
double N = 256;
double b = 2;
// ----------------------------------

int color[height][width] = {0};
int final_color[height][width] = {0};
int num_workers = 0;
int slice_width = 0;
FILE *fp;

// Computations cal_pixel 
int cal_pixel(double complex d, double b, double N){
	int count = 1;
	double complex z = 0;
	while ((cabs(z)<b) && (count < N)){
		z = (z*z)+ d;
		count= count + 1;
	}
	return count;
}

// Create color.txt for plotting of result 
void CreateFile(){
    fp = fopen("color.txt", "w");
    for (int i = 0; i<height; i++){
        for (int j=0; j<width; j++){
            fprintf(fp,"%d", final_color[i][j]);
            fprintf(fp," ");
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

// Add color array to final array 
void ArrayAssignment(int color[height][width]){
	for (int i = 0; i < height;i++){
     	for (int j = 0; j < width; j++){
     		final_color[i][j] += color[i][j];
     	}
     }
}

// Mandelbrot calculations 
void CalculateColumn(int offset){
    
    double dx = 2*b/(width-1);
    double dy = 2*b/(height-1);
    double dreal = 0;
    double dimag = 0;
    double complex d = 0;
    int yoff = 0;

    for (int x = 0; x < slice_width; x++){
        dreal = (x+offset)*dx - b;
        for (int y = 0; y < height; y++){
            dimag = (y+yoff)*dy-b;
            d = dreal + 1I*dimag;
            color[(x+offset)][y] = cal_pixel(d,b,N);
        }
    }
}

// Main function 
int main (int argc, char **argv) {

    int offset = 0;
    int size = 0;
    int rank = 0;
    int rc = 0;
    int id, tasks, tag, i;
    int source;
    
    // Initiate MPI environment, obtain number of workers and rank.
    rc = MPI_Init(NULL,NULL);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Number of workers is the total minus one (the master). Each worker will have a column to calculate.
    num_workers = size-1;
    slice_width = width/num_workers; 

    // Master
    if (rank == 0){
        // Send offset in y to each worker.
        for(int dest = 1; dest <= num_workers; dest++){
            MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
            offset = offset + slice_width;
        }

        // Recieve all results from workers and join arrays for final result.
        for(int i = 1; i <= num_workers; i++){
            source = i;
            MPI_Recv(&color,height*width,MPI_INT,source,1,MPI_COMM_WORLD,&status);
           	ArrayAssignment(color);
        }

        // Create final file with results.
     	CreateFile();
    
    // Workers
    }else{
        source = 0;
        // Recieve offset in y from which to start calculating, and send result array back.
        MPI_Recv(&offset, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
        CalculateColumn(offset);
        MPI_Send(&color, height*width, MPI_INT, 0, 1, MPI_COMM_WORLD);

    }
    MPI_Finalize();
}
