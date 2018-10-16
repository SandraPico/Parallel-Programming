/***********************************************************+
* @author: Miquel Sven Larsson and Sandra Pico Oristrell
* @data: 7 May 2018
* @description: Sequential region growing
* @reference: http://notmatthancock.github.io/2017/10/09/region-growing-wrapping-c.html
* @Comand to execute and compile it: gcc region_growing_sequential.c stack.h && ./a.out
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "stack.h"
#include <math.h>

/*************************************+
* Global Variables
***************************************/

//Thresdhold for the inclusion criteria
float threshold = 0.1;
//Size of the total image array
int size_array = 20*20*20;

/***********************************************************+
* @Name: coord2index
* @Description: convert 3D coordinates to 1D index value
* @Parameters:(coordinates(i,j,k) and total dimensions(m,n,p))
* @Output: 1D index value
************************************************************/
int coord2index(int i, int j, int k, int m, int n, int p) {
    return  n*p*i + p*j + k;
}

/***********************************************************+*********
* @Name: index2coord
* @Description: convert 1D index value to 3D cordinates
* @Parameters:(1D index, total dimensions(m,n,p) and coordinates array)
* @Output: 1D index value
*********************************************************************/
void index2coord(int index, int m , int n, int p, int* coord){        
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
void AddNeighbors(int i, int j, int k, int m, int n, int p, stack * needs_check, int * regions) {
    
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
bool InclusionCriteria(int i,int j,int k,int si,int sj,int sk,double * img, int m, int n, int p){
    
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
void RegionGrowing(int m, int n, int p, double * img, int si, int sj, int sk, int * regions) {
    
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
    regions[index] = r;

    // Add the neighbors of the seed point to the stack.
    AddNeighbors(si, sj, sk, m, n, p, &needs_check, regions);

    //Until there are elements in needs_check structure...
    while (needs_check.n_elements > 0) {       
        // Pop a coordinate from the stack.
        el = stack_pop(&needs_check);
        i = el->i; j = el->j; k = el->k;
        free(el);

        //Calculate the 1D index of the point to check.
        index = coord2index(i, j, k, m, n, p);

        //Is already checked?
        if (regions[index] != 0) continue;

        //Check if the inclusion criteria is satisfied
        if (InclusionCriteria(i,j,k,si,sj,sk,img,m,n,p)){
            regions[index] = r;
            AddNeighbors(i, j, k, m, n, p, &needs_check, regions);
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
int Count_0(int* regions){
    int count = 0;
    for (int i = 0; i < size_array; i++){
        if (regions[i] == 0){
            count = count + 1;
        }
    }
    return count;
}


/***********************************************************+****************************************
* @Name: PrintRegions
* @Description: Function that prints the total number of elements in each region
* @Parameters:  regions array, image.
* @Output: Print the region information
****************************************************************************************************/
void PrintRegions(int* regions, double * img){
    
    int max = -1;
    for (int i = 0; i < size_array; i++){
        if (regions[i] > max){
            max = regions[i];
        }
    }

    for (int k=1; k<max+1; k++){
        int c = 0;
        for (int j = 0; j < size_array; j++){
            if (regions[j] == k){
                c++;
            }
        }
        printf("Region: %d - Elements: %d \n", k, c);
    }

}

/***********************************************************+****************************************
* @Name: CreateFile
* @Description: Generate the results file based on regions information
* @Parameters:  Regions array
* @Output: Generate Results.txt file
****************************************************************************************************/
void CreateFile(int* regions){
    FILE * f;
    f = fopen("Results_Sequential.txt","a");
    for (int i = 0; i < size_array; i ++){
        fprintf(f, "%d\n", regions[i]);
    }
    fclose(f);
}

/***********************************************************+****************************************
* @Name: Main
* @Description: Execute the region growing algorithm based on the file data.
****************************************************************************************************/
int main() {
    
    //Dimensions of the 3D image array.
    int n = 20;
    int m = 20;
    int p = 20;

    //Read from the file
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

    //Create and initialize the regions array
    int regions[n*m*p];
    for (int i = 0 ; i < n*m*p; i ++){
        regions[i] = 0;
    }

    //Seed coordinates array
    int seed_coord[3];
    
    //Count how many pixels without a region there are in the regions array.
    int count_zero = Count_0(regions);
    while (count_zero != 0){

        //Seed point calculation

        //Find index zero
        int index_zero = rand()%count_zero;
        //Find real index into the regions array
        int index_real = FindIndex(index_zero,regions);
        //Convert index into 3D coordinates of the point
        index2coord(index_real,m,n,p,seed_coord);
        
        //Region growing algorithm
        RegionGrowing(m,n,p,img,seed_coord[0],seed_coord[1],seed_coord[2],regions); 
        
        //Count how many zeros do we have again
        count_zero = Count_0(regions);
    }

    //Print regions information
    PrintRegions(regions,img);

    //Generate the output file
    CreateFile(regions);
    return 0;
}







