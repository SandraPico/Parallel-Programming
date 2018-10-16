/***********************************************************************************************+
* @author: Miquel Sven Larsson and Sandra Pico Oristrell
* @data: 7 May 2018
* @description: Parallel region growing
* @Use makefile to compile and: mpirun --oversubscribe -np 5 parallel_region_growing to execute.
***********************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "stack.h"
#include <math.h>
#include <mpi.h>

/*************************************+
* Global Variables
***************************************/
//Threshold
float threshold = 0.1;
//Dimensions image 
int n = 20;
int m = 20;
int p = 20;
int size_array = 20*20*20;

//MPI global variables
int size = 0;
int tag = 10;
int rank = 0;
MPI_Status status;

//MPI grid global variables
int dims[2];
int coords[2];
int periods[2] = {0,0};

//Directions in the MPI grid
int north;
int south;
int east;
int west;

// Cartesian communicator
MPI_Comm cart_comm;             


/***********************************************************+
* @Name: coord2index
* @Description: convert 3D coordinates to 1D index value
* @Parameters:(coordinates(i,j,k) and total local dimensions(m,n,p))
* @Output: 1D index value
************************************************************/
int coord2index(int i, int j, int k, int m, int n, int p) {
    return n*p*i + p*j + k;
}


/***********************************************************+*********
* @Name: index2coord
* @Description: convert 1D index value to 3D cordinates
* @Parameters:(1D index, total dimensions(m,n,p) and coordinates array)
* @Output: 1D index value
*********************************************************************/
void index2coord(int index,int* coord, int m, int n, int p){      
    coord[0] = floor(index/(n*p));
    coord[1] = floor((index-(coord[0]*n*p))/p);
    coord[2] = index - (coord[0]*n*p) - (coord[1]*p);
}


/***********************************************************+************************
* @Name: AddNeighbors
* @Description: add into the stack structure all the neighbours of (i,j,k) point
* @Parameters: point coordinates(i,j,k), total dimensions(m,n,p) and regions array.
* @Output: add the neighbours points into the needs_check stack structure
************************************************************************************/
void AddNeighbors(int i, int j, int k, stack * needs_check, int * regions, int m, int n, int p) {
        //Index in the regions array.
    int index;
    //Element variable for the stack structure.
    element * el;

    /******************************************************
    * The function needs to check the borders of the image 
    * 0 <= i <= m
    * 0 <= j <= n
    * 0 <= k <= p
    ******************************************************/

    /*****************************************************
    * i  borders
    ******************************************************/
    if (i >= 1) {
        index = coord2index(i-1,j,k,m,n,p);
        if (regions[index] == 0) {
            el = (element *) malloc(sizeof(element));
            element_init(el, i-1, j, k);
            stack_push(needs_check, el);
        }
    }

    if (i < m-1) {
        index = coord2index(i+1,j,k,m,n,p);
        if (regions[index] == 0) {
            el = (element *) malloc(sizeof(element));
            element_init(el, i+1, j, k);
            stack_push(needs_check, el);
        }
    }

    /*****************************************************
    * j  borders
    ******************************************************/
    if (j >= 1) {
        index = coord2index(i,j-1,k,m,n,p);
        if (regions[index] == 0) {
            el = (element *) malloc(sizeof(element));
            element_init(el, i, j-1, k);
            stack_push(needs_check, el);
        }
    }

    if (j < n-1) {
        index = coord2index(i,j+1,k,m,n,p);
        if (regions[index] == 0) {
            el = (element *) malloc(sizeof(element));
            element_init(el, i, j+1, k);
            stack_push(needs_check, el);
        }
    }

    /*************************************************
    * k borders
    *************************************************/
    if (k >= 1) {
        index = coord2index(i,j,k-1,m,n,p);
        if (regions[index] == 0) {
            el = (element *) malloc(sizeof(element));
            element_init(el, i, j, k-1);
            stack_push(needs_check, el);
        }
    }

    if (k < p-1) {
        index = coord2index(i,j,k+1,m,n,p);
        if (regions[index] == 0) {
            el = (element *) malloc(sizeof(element));
            element_init(el, i, j, k+1);
            stack_push(needs_check, el);
        }
    }
}



/***********************************************************+****************************************
* @Name: InclusionCriteria
* @Description: check if two different points fulfill the inclusion criteria(based on the threshold)
* @Parameters: point coordinates(i,j,k), seed coordinates(si,sj,sk), image and total dimensions(m,n,p).
* @Output: Boolean variable indicating if the inclusion criteria is fulfilled.
****************************************************************************************************/
bool InclusionCriteria(int i,int j,int k,int si,int sj,int sk,double * img, int m, int n , int p){
    
    //Seed index
    int index_seed = coord2index(si,sj,sk,m,n,p); 
    //Point index
    int index_pixel = coord2index(i,j,k,m,n,p); 

    //Inclusion criteria based on threshold comparison.
    if (fabs(img[index_pixel] - img[index_seed])< threshold){
        return true;
    }else{
        return false;
    }
}

/***********************************************************+****************************************
* @Name: RegionGrowing 
* @Description: Execute the region growing algorithm based on the initial seed point and the image
* @Parameters: total dimensions(m,n,p), image, seed point(si,sj,sk) and regions array.
* @Output: regions array updated with new region values.
****************************************************************************************************/
void RegionGrowing(double * img, int si, int sj, int sk, int * regions, int m, int n , int p, int size_array) {
    
    //Check the region number in the region array
    int max = -1;
    for (int i = 0; i < size_array; i++){
        if (regions[i] > max){
            max = regions[i];
        }
    }
    //region number
    int r = max+1;

    int index, i, j, k;
    stack needs_check; 
    element * el;

    // Initialize the stack (needs_check structure)
    stack_init(&needs_check);

    //Calculate the index of the seed and update the regions array of the seed position.
    index = coord2index(si,sj,sk,m,n,p);
    regions[index] = ((rank*100) + r);

    // Add the neighbors of the seed point to the stack.
    AddNeighbors(si, sj, sk, &needs_check,regions,m,n,p);


    //Until there are elements in needs_check structure...
    while (needs_check.n_elements > 0) {
        // Pop a coordinate from the stack.
        el = stack_pop(&needs_check);
        i = el->i; j = el->j; k = el->k;
        free(el);

        //Calculate the 1D index of the point to check.
        index = coord2index(i, j, k,m,n,p);

        //Is already checked?
        if (regions[index] != 0) continue;

        //Check if the inclusion criteria is satisfied
        if (InclusionCriteria(i,j,k,si,sj,sk,img,m,n,p)){
            regions[index] = ((rank*100) + r);
            AddNeighbors(i, j, k, &needs_check, regions,m,n,p);
        }
    }
}


/***********************************************************+****************************************
* @Name: FindIndex
* @Description: Find the proper index in the regions array based on the zero index  
* @Parameters: index zero and regions array
* @Output: Real index in the regions array
****************************************************************************************************/
int FindIndex(int index_zero,int* regions){
    int count_zeros = 0;
    int i = 0;
    index_zero ++;
    while (count_zeros < index_zero){
        if (regions[i] == 0){
            count_zeros++;
        }
        i++;
    }
    return i-1;
}

/***********************************************************+****************************************
* @Name: Count_0
* @Description: Count how many 0 (pixels without a region assigned) there's in regions array
* @Parameters: regions array
* @Output: int counter number.
****************************************************************************************************/
int Count_0(int* regions, int local_size){
    int count = 0;
    for (int i = 0; i < local_size; i++){
        if (regions[i] == 0){
            count = count + 1;
        }
    }
    return count;
}



/***********************************************************+****************************************
* @Name: CreateFile
* @Description: Generate the results file based on regions information
* @Parameters:  Regions array
* @Output: Generate Results.txt file
****************************************************************************************************/
void CreateFile(int* regions){
    FILE * f;
    f = fopen("Results_Parallel.txt","a");
    for (int i = 0; i < size_array; i ++){
        fprintf(f, "%d\n", regions[i]);
    }
    fclose(f);
}


/***********************************************************+****************************************
* @Name: MPI_init
* @Description: Initialize the MPI Grid and the MPI Communicator
* @Output: MPI set up
****************************************************************************************************/
void mpi_init(){

    //Initilialize MPI library
    int rc = 0;
    rc = MPI_Init(NULL,NULL);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //Create the dimensions for the grid , dims[1] = rows, dims[0] = columns.
    MPI_Dims_create(size, 2, dims);

    //Create the MPI coord grid (based on the dimensions)
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&cart_comm);
    //Assign coord (in the MPI grid) per each rank.
    MPI_Cart_coords(cart_comm,rank,2,coords);


    //To know if the process have a north,south,west and/or east communications.
    //1,1 : Direction 1: column, 1 jump.
    //0,1 : Direction 0: row, 1 jump.
    MPI_Cart_shift(cart_comm,1,1,&north,&south);
    MPI_Cart_shift(cart_comm,0,1,&west,&east);
}



/***********************************************************+***************************************************
* @Name: CombineRegions
* @Description: Exchange information between the north-south and east-west processes (compare the border image)
* @Parameters:  regions array, local dimensions (m,n,p) and image per each process
***************************************************************************************************************/
void CombineRegions(int * regions,int m,int n, int p, double * sub_img){
    
    //Receive from west...
    if (west>= 0){  
       
        //Receive neighbour image border
        double neigh_west_image[n*p];
        MPI_Recv(neigh_west_image,n*p,MPI_DOUBLE,west,tag,cart_comm,&status);
        
        //Receive neighbour regions border
        int neigh_west_region[n*p];
        MPI_Recv(neigh_west_region,n*p,MPI_INTEGER,west,tag,cart_comm,&status);

        //Update local image and region borders
        double local_west_image[n*p];
        int local_west_region[n*p];
        for(int i = 0; i < n*p; i++){
            local_west_image[i]  = sub_img[i];
            local_west_region[i] = regions[i];
        }

        //Compare the two borders (neighbour and local) using the InclusionCriteria.
        int aux = 0;
        for(int i = 0; i < n*p; i++){
            if (neigh_west_region[i] != local_west_region[i]){
                if (fabs(neigh_west_image[i]-local_west_image[i]) < threshold){
                    //Update border array
                    aux = local_west_region[i];
                    for (int j=0; j< n*p; j++){
                        if (local_west_region[j]==aux){
                            local_west_region[j]=neigh_west_region[i];
                        }
                    }
                    //Update whole local region
                    for (int j=0; j< m*n*p; j++){
                        if (regions[j]==aux){
                            regions[j]=neigh_west_region[i];
                        }
                    }
                }
            }
        }
    }


    //Receive from north...
    if (north >= 0){

        //Receive neighbour image border
        double neigh_north_image[m*p];
        MPI_Recv(neigh_north_image,m*p,MPI_DOUBLE,north,tag,cart_comm,&status);  

        //Receive neighbour regions border
        int neigh_north_region[m*p];
        MPI_Recv(neigh_north_region,m*p,MPI_INTEGER,north,tag,cart_comm,&status);

        //Update local image and region borders
        double local_north_image[m*p];
        int local_north_region[m*p];
        int coord[3];
        int k = 0;
        for(int i = 0; i < m*n*p; i++){
            index2coord(i,coord,m,n,p);
            if (coord[1] == 0){
                local_north_image[k] = sub_img[i];
                local_north_region[k] = regions[i];
                k++;
            }
        }

        //Compare the two borders (neighbour and local) using the InclusionCriteria.
        int aux = 0;
        for(int i = 0; i < m*p; i++){
            if (neigh_north_region[i] != local_north_region[i]){
                if (fabs(neigh_north_image[i]-local_north_image[i])< threshold){
                    //Update border array
                    aux = local_north_region[i];
                    for (int j=0; j< m*p; j++){
                        if (local_north_region[j]==aux){
                            local_north_region[j]=neigh_north_region[i];
                        }
                    }
                    //Update whole local region
                    for (int j=0; j< m*n*p; j++){
                        if (regions[j]==aux){
                            regions[j]=neigh_north_region[i];
                        }
                    }
                }
            }
        }
    } 

    //Send to east process..
    if (east >= 0){

        //Update local image and region border array
        double local_east_image[n*p];
        int local_east_region[n*p];
        int k = 0;
        for(int i = (m-1)*n*p; i < m*n*p; i++ ){
            local_east_image[k] = sub_img[i];
            local_east_region[k] = regions[i];
            k++;
        }

        //Send it to east process..
        MPI_Send(local_east_image,n*p,MPI_DOUBLE,east,tag,cart_comm);
        MPI_Send(local_east_region,n*p,MPI_INTEGER,east,tag,cart_comm);
    }


    if (south >= 0){
        //Update local image and region border array
        double local_south_image[m*p];
        int local_south_region[m*p];
        int coord[3];
        int k = 0;
        for(int i = 0; i < m*n*p; i++){        
            index2coord(i,coord,m,n,p);
            if (coord[1] == n-1){
                local_south_image[k] = sub_img[i];
                local_south_region[k] = regions[i];
                k++;
            }
        }

        //Send it to south process..
        MPI_Send(local_south_image,m*p,MPI_DOUBLE,south,tag,cart_comm);
        MPI_Send(local_south_region,m*p,MPI_INTEGER,south,tag,cart_comm);
    }
}


/***********************************************************+***************************************************
* @Name: Distribute_Image
* @Description: Divide the image and send it to all the processes.
* @Parameters:  image and sub_image
* @Output: Every process will have their own sub-image.
***************************************************************************************************************/
void Distribute_image(double * img, double * sub_image){
    
    int rows = (size);
    int columns = (int)(size_array/(size));
    
    //If it's the master, divide the image and distribute it to each node...
    if (rank == 0){        
        
        double sub_blocks[rows][columns];
        //m - x, n - y , p - z
        int m_split = (int)(m/dims[0]);
        int n_split = (int)(n/dims[1]);
        int nodes = dims[0]*dims[1];
        int node = 0;
        int i = 0;
        for(int x_min = 0; x_min < m; x_min = x_min + m_split){
            for(int y_min = 0; y_min < n; y_min = y_min + n_split){
                for (int x = 0; x < m; x++){
                    for (int y = 0; y < n ; y++){
                        for (int z =0 ; z < p; z++){
                            if((x>= x_min) && (x<(x_min+m_split)) && (y>=y_min) && (y < (y_min + n_split))){
                                int index = coord2index(x,y,z,m,n,p);
                                sub_blocks[node][i] = img[index];
                                i++;
                            }
                        }
                    }
                }
                node++;
                i = 0;
            }
        }
        
        //Master send parts of the image to the other processors.
        for (int i = 1; i < nodes; i++){
            //data, count (data size), type, destination, tag, communicator.
            MPI_Send(sub_blocks[i],columns,MPI_DOUBLE,i,tag,cart_comm);
        }
        
        //Update master information.
        for (int i = 0; i < columns; i++){
            sub_image[i] = sub_blocks[0][i];
        }
    }else{
        //Receive sub-image from the master process...
        MPI_Recv(sub_image,columns,MPI_DOUBLE,0,tag,cart_comm,&status);
    }
}

/***********************************************************+***************************************************
* @Name: Report_results
* @Description: Report the final results (how many elements per each region)
* @Parameters:  final regions array
***************************************************************************************************************/
void Report_results(int* final_regions){
    
    //Max number of regions in an image : 100
    int distinct_regions[100];
    int count_elements[100];

    for (int i = 0; i<100; i++){
        distinct_regions[i]=-1;
        count_elements[i]=0;
    }

    for (int i = 0; i<size_array; i++){
        bool in_list_already= false;
        int j=0;
        while (distinct_regions[j]>0){
            if (distinct_regions[j]==final_regions[i]){
                in_list_already = true;
                break;
            }
            j++;
        }
        if (in_list_already==false){
            distinct_regions[j]=final_regions[i];
        }
        count_elements[j]++;
    }

    int i=0;

    //Print final results
    while (distinct_regions[i]>0){
        printf("Region: %d - Num Elements: %d \n", distinct_regions[i], count_elements[i]);
        i++;
    }
    
}

/***********************************************************+***************************************************+********
* @Name: Gather_Regions
* @Description: Gather all the information from the processes into a final array. Call the generate and report functions
* @Parameters:  local regions (per each process) and local size.
* @Output: Every process will have their own sub-image.
***************************************************************************************************************+*********/
void Gather_Regions(int* regions, int size_local){    

    //If its the master process..
    if(rank == 0){

        //Receive from all the processes and generate the final regions.
        int final_regions[n*m*p];
        int final_sub_blocks[size][size_local];
        for(int i = 0; i < size; i++){
            int region_local[size_local];
            if (i>0){
                MPI_Recv(region_local,size_local,MPI_INTEGER,i,tag,cart_comm,&status);
                for (int j = 0; j< size_local; j++){
                    final_sub_blocks[i][j] = region_local[j];
                }
            }else{
                for (int j = 0; j< size_local; j++){
                    final_sub_blocks[i][j] = regions[j];
                }
            }
        }
        int m_split = (int)(m/dims[0]);
        int n_split = (int)(n/dims[1]);
        int nodes = dims[0]*dims[1];
        int node = 0;
        int i = 0;

        for(int x_min = 0; x_min < m; x_min = x_min + m_split){
            for(int y_min = 0; y_min < n; y_min = y_min + n_split){
                for (int x = 0; x < m; x++){
                    for (int y = 0; y < n ; y++){
                        for (int z =0 ; z < p; z++){
                            if((x>= x_min) && (x<(x_min+m_split)) && (y>=y_min) && (y < (y_min + n_split))){
                                int index = coord2index(x,y,z,m,n,p);
                                final_regions[index] = final_sub_blocks[node][i];
                                i++;
                            }
                        }
                    }
                }
                node++;
                i = 0;
            }
        }

        //Create the file
        CreateFile(final_regions);
        //Report results
        Report_results(final_regions);
    }else{
        //Send the sub-regions array to the master process...(rank = 0)
        MPI_Send(regions,size_local, MPI_INTEGER,0,tag,cart_comm);
    }
}


/***********************************************************+****************************************
* @Name: Main
* @Description: Execute the region growing algorithm based on the file data.
****************************************************************************************************/
int main(int argc, char** argv) {

    //Read data
    int array[n*m*p];
    FILE *fptr;
    if ((fptr = fopen("cube_ball.txt", "r")) == NULL){
        printf("Error! opening file");
        exit(1);         
    }
    for (int i = 0; i < n*m*p; i++){
        fscanf(fptr, "%d", &array[i]);
    }
    fclose(fptr);
    double img[n*m*p];
    for(int i = 0; i < n*m*p; i++){
        img[i] = (double)array[i];
    }

    //MPI Init
    mpi_init();
    
    int n_local = (int)(n/dims[1]);
    int m_local = (int)(m/dims[0]);
    int p_local = p;

    //Initial data size
    if (rank==0){
        printf("Image Dimensions (m, n, p): %d, %d, %d\n", m,n,p);
        printf("MPI Dimensions: (%d,%d)\n\n",dims[1],dims[0]);
    }

    //Calculate sub-image and sub-image dimensions
    int columns = (int)(size_array/(size));
    double sub_image[size_array/(size)];
    Distribute_image(img,sub_image);

    //Regions/process
    int regions[n_local*m_local*p_local];
    int local_size = n_local*m_local*p_local;
    for (int i = 0 ; i < local_size; i ++){
        regions[i] = 0;
    }

    //Region growing call per each process 
    int coord[3];
    int count_zero = Count_0(regions,local_size);
    while (count_zero != 0){
        int index_zero = rand()%count_zero;
        int index_real = FindIndex(index_zero,regions);     
        index2coord(index_real,coord, m_local, n_local, p_local);
        RegionGrowing(sub_image,coord[0],coord[1],coord[2],regions, m_local, n_local, p_local, local_size); 
        count_zero = Count_0(regions,local_size);
    }
    //Wait to all processes to finish.
    MPI_Barrier(cart_comm);
    
    //All the processes have their own sub-regions. We just need to combine them (final regions)
    //Compare borders of each sub-image using the inclusion criteria.
    CombineRegions(regions,m_local,n_local, p_local,sub_image);
    
    //Wait
    MPI_Barrier(cart_comm);

    //Gather information to the master and report the final results.
    Gather_Regions(regions,local_size);    
    MPI_Finalize();
    return 0;
}







