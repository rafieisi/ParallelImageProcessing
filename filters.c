/* ------------
 * This code is provided solely for the personal and private use of 
 * students taking the CSC367H5 course at the University of Toronto.
 * Copying for purposes other than this use is expressly prohibited. 
 * All forms of distribution of this code, whether as given or with 
 * any changes, are expressly prohibited. 
 * 
 * Authors: Bogdan Simion, Felipe de Azevedo Piovezan
 * 
 * All of the files in this directory and all subdirectories are:
 * Copyright (c) 2019 Bogdan Simion
 * -------------
*/

#include "filters.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct common_work_t
{
    const filter *f;
    const int32_t *original_image;
    int32_t *output_image;
    int32_t *inter;
    int32_t width;
    int32_t height;
    int32_t min;
    int32_t max; 
    parallel_method method;
    int32_t max_threads;
    int32_t padding;
    int32_t size;
    int32_t work_chunk;
    pthread_barrier_t barrier;
    pthread_mutex_t lock;
    int32_t chunknumber;
    pthread_mutex_t lastchunk; 
    int32_t lastchunkstarted;
} common_work;

typedef struct work_t
{
    common_work *common;
    int32_t id;
    int min;
    int max;
} work;
/************** FILTER CONSTANTS*****************/
/* laplacian */
int8_t lp3_m[] =
    {
        0, 1, 0,
        1, -4, 1,
        0, 1, 0,
    };
filter lp3_f = {3, lp3_m};

int8_t lp5_m[] =
    {
        -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1,
        -1, -1, 24, -1, -1,
        -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1,
    };
filter lp5_f = {5, lp5_m};

/* Laplacian of gaussian */
int8_t log_m[] =
    {
        0, 1, 1, 2, 2, 2, 1, 1, 0,
        1, 2, 4, 5, 5, 5, 4, 2, 1,
        1, 4, 5, 3, 0, 3, 5, 4, 1,
        2, 5, 3, -12, -24, -12, 3, 5, 2,
        2, 5, 0, -24, -40, -24, 0, 5, 2,
        2, 5, 3, -12, -24, -12, 3, 5, 2,
        1, 4, 5, 3, 0, 3, 5, 4, 1,
        1, 2, 4, 5, 5, 5, 4, 2, 1,
        0, 1, 1, 2, 2, 2, 1, 1, 0,
    };
filter log_f = {9, log_m};

/* Identity */
int8_t identity_m[] = {1};
filter identity_f = {1, identity_m};

filter *builtin_filters[NUM_FILTERS] = {&lp3_f, &lp5_f, &log_f, &identity_f};

/* Normalizes a pixel given the smallest     common->barrier=NULL;d largest integer values
 * in the image */
void normalize_pixel(int32_t *target, int32_t pixel_idx, int32_t smallest, 
        int32_t largest)
{
    if (smallest == largest)
    {
        return;
    }
    
    target[pixel_idx] = ((target[pixel_idx] - smallest) * 255) / (largest - smallest);
}
/*************** COMMON WORK ***********************/
/* Process a single pixel and returns the value of processed pixel
 * TODO: you don't have to implement/use this function, but this is a hint
 * on how to reuse your code.
 * */

/*
 * translator will map the index of the target matrix that does not have any zero paddings to 
 * the intermeddiate matrix that contains the intermeddiate matrix
 * 
 * returns the index of the inter matrix
 */
int translator(int index, int width, int height, int padding){
    int starting_pad = width+3*padding; //the starting_pad is the amount of zeroes in the beginning of the matrix
    int jump = (int)((float)index/width); // to see how much we should jumb to get to the next row without reading the zeroes
    return index + (jump*2*padding) + starting_pad; // each jump contains 2 padding(left and right one)
}

/*
 *applies filter to one pixel(index) of the matrix and returns the value of the Laplacian filter
 */
int32_t applyfilter(const filter *f, int *inter,
        int32_t width, int32_t height,
        int rowindex, int columnindex, int padding)
{
    int32_t sum = 0;
    //since the inter matrix already contains the zero matrix, we just need to look at the neighborings pixels
    for(int col=(-1)*padding;col<=padding;col++){
        for(int row=(-1)*padding; row<=padding;row++){
            //adding row will take a look at the left and the right neighbors
            //adding ((width+padding*2)*col will look at the top and bottom neighbors
            sum += (inter[translator(rowindex*width+columnindex,width,height,padding)+row+((width+padding*2)*col)]*f->matrix[padding+row+(f->dimension*(col+padding))]);
        }
        
    }
    return sum;
}
/*********SEQUENTIAL IMPLEMENTATIONS ***************/
/* TODO: your sequential implementation goes here.
 * IMPORTANT: you must test this thoroughly with lots of corner cases and 
 * check against your own manual calculations on paper, to make sure that your code
 * produces the correct image.
 * Correctness is CRUCIAL here, especially if you re-use this code for filtering 
 * pieces of the image in your parallel implementations! 
 */
void apply_filter2d(const filter *f, 
        const int32_t *original, int32_t *target,
        int32_t width, int32_t height)
{
    int padding=(int)(f->dimension/2);
    //we are going to use a matrix that also includes the zero paddings
    int32_t *inter=malloc(sizeof(int32_t)*(width+(padding*2))*(height+(padding*2)));
    target = malloc(sizeof(int32_t)*(width)*(height));
    //making all of the elements in the inter equal to 0
    for (int i=0;i<(width+padding+padding)*(height+padding+padding);i++){
        inter[i]=0;
    }
    
    //we have to check for the index that are going to be 0 and ignore them
    int index=0;//for keeping track of the index of the original matrix during copying
    for (int i=0;i<(width+padding+padding)*(height+padding+padding);i+=1){
        if((i<(width+padding+padding)*padding)||i>=((width+padding+padding)*(height+padding))){
            continue;//takes care of top and bottom rows
        }
        if(i%(width+padding+padding)<padding){
            continue;//takes care of the left zero paddings
        }
        if(i%(width+padding+padding)>=padding+width){
            continue;//takes care of the right zero paddings
        }
        
        inter[i]=original[index];
        index++;
    }
    
    int32_t max=0;
    int32_t min=0;

    unsigned short first_iteration = 1; //for initializing the max and min to the first filter value for later checks
    //looping through each index to find its filter value
    for (int i=0; i<height; i++){
        for(int j=0;j<width; j++){
            int sum = applyfilter(f,inter,width,height,i,j,padding);
            
            target[i*width+j] = sum;
            if(first_iteration){
                min = sum;
                max = sum;
                first_iteration = 0;
                continue;
            }
            
            if(max < sum){
                max = sum;
            }

            if(min > sum){
                min = sum;
            }
        }
    }
    //after finishing the filter process and finding the max and min, we want to normalize all of the pixels
    for (int i=0; i<height; i++){
        for(int j=0;j<width; j++){
            normalize_pixel(target, i*width+j, min, max);
        }
    }
    free(inter);
}

/***************** WORK QUEUE *******************/
/* TODO: you don't have to implement this. It is just a suggestion for the
 * organization of the code.
 */
void* queue_work(void *works)
{
    work * w=works;    
    const filter *f=w->common->f;
    int *inter=w->common->inter;
    int padding=w->common->padding;
    int height=w->common->height;
    int width=w->common->width;
    int *target=w->common->output_image;
    //the total number of chunks in one row
    //note that we have to take the ceiling for the case that it does not fit one chunk perfectly
    int chunksperrow= (height+w->common->work_chunk-1)/w->common->work_chunk;
    //the total number of chunks in one col
    //note that we have to take the ceiling for the case that it does not fit one chunk perfectly
    int chunkspercol= (width+w->common->work_chunk-1)/w->common->work_chunk;
    int numofchunks=chunksperrow*chunkspercol;
    //each thread will be assigned to a specific chunk number
    int chunknumber=w->id;
    int chunk=w->common->work_chunk;
    pthread_mutex_lock(&w->common->lastchunk);
    w->common->lastchunkstarted=chunknumber;
    pthread_mutex_unlock(&w->common->lastchunk);

    int min=0;
    int max=0;
    unsigned short first_iteration = 1; //for initializing the max and min to the first filter value for later checks
    //this loops continue as longs as there are chunks left for the threads to calculate
    while (chunknumber<numofchunks){
        for( int row= chunknumber*chunksperrow; row < (chunknumber+1)*chunksperrow * chunk;row++){
            //this is specially for the case that the chunk cannot fit perfectly
            if(row > height){
                break;
            }
            for (int col= chunknumber*chunkspercol;col < (chunknumber+1)*chunkspercol * chunk;col++){
                if ((row < height) && (col < width)){
                    target[row*width+col] = applyfilter(f,inter,width,height,row,col,padding);
                    if(first_iteration){
                        min = target[row*width+col];
                        max = target[row*width+col];
                        first_iteration = 0;
                        continue;
                    }

                    if(max < target[row*width+col]){
                        max = target[row*width+col];
                    }

                    if(min > target[row*width+col]){
                        min = target[row*width+col];
                    }  
                }
            }
        }

        pthread_mutex_lock(&w->common->lastchunk);
        w->common->lastchunkstarted++;//will process the next chunk left
        chunknumber=w->common->lastchunkstarted;
        pthread_mutex_unlock(&w->common->lastchunk);  
    }
    pthread_mutex_lock(&w->common->lock);
            
    if(w->common->min > min){
        w->common->min = min;
    }

    if(w->common->max < max){
        w->common->max = max;
    }
    pthread_mutex_unlock(&w->common->lock); 

    pthread_barrier_wait(&w->common->barrier);
    
    //same process as above in order to normalize each pixel
    chunknumber=w->id;
    pthread_mutex_lock(&w->common->lastchunk);
    w->common->lastchunkstarted=chunknumber;
    pthread_mutex_unlock(&w->common->lastchunk);
    
    while (chunknumber<numofchunks){

        for( int row= chunknumber*chunksperrow; row < (chunknumber+1)*chunksperrow * chunk;row++){
            if(row > height){break;}
            for (int col= chunknumber*chunkspercol;col < (chunknumber+1)*chunkspercol * chunk;col++){
                if ((row < height) && (col < width)){
                    normalize_pixel(target, row*width+col, min, max);
                }
            }
        }
    
    
        pthread_mutex_lock(&w->common->lastchunk);
        w->common->lastchunkstarted++;
        chunknumber=w->common->lastchunkstarted;
        pthread_mutex_unlock(&w->common->lastchunk); 
    }
     
    return NULL;
}

/****************** ROW/COLUMN SHARDING ************/
/* TODO: you don't have to implement this. It is just a suggestion for the
 * organization of the code.
 */

/* Recall that, once the filter is applied, all threads need to wait for
 * each other to finish before computing the smallest/largets elements
 * in the resulting matrix. To accomplish that, we declare a barrier variable:
 *      pthread_barrier_t barrier;
 * And then initialize it specifying the number of threads that need to call
 * wait() on it:
 *      pthread_barrier_init(&barrier, NULL, num_threads);
 * Once a thread has finished applying the filter, it waits for the other
 * threads by calling:
 *      pthread_barrier_wait(&barrier);
 * This function only returns after *num_threads* threads have called it.
 */
void* sharding_work(void *works)
{
    /* Your algorithm is essentially:
     *  1- Apply the filter on the image
     *  2- Wait for all threads to do the same
     *  3- Calculate global smallest/largest elements on the resulting image
     *  4- Scale back the pixels of the image. For the non work queue
     *      implementations, each thread should scale the same pixels
     *      that it worked on step 1.
     */
    work * w=works;
    int thread_num=w->id;
    int width=w->common->width;
    int height=w->common->height;
    int padding=w->common->padding;
    int *inter=w->common->inter;
    parallel_method method = w->common->method;
    int size = w->common->size;
    const filter *f=w->common->f;
    int *target=w->common->output_image;
    //start and end will keep track of the starting index and the ending index of their work either row or col
    int start = 0;
    int end = 0;
    int increment = 0;
    
    int min=0;
    int max=0;
    unsigned short first_iteration = 1;//for initializing the max and min to the first filter value for later checks
    if(method == SHARDED_ROWS){
        start = thread_num * size;
        end = (thread_num+1) * size; //stops at the next thread work
        increment = 1;
        //since this is sharded rows, we want to divide the rows and go through each col(0 to width) 
        for (int i = start ; i < end ; i += increment){
            for(int j=0;j<width;j++){
                if (i < height){//if the last row had more rows than height
                    target[i*width+j] = applyfilter(f,inter,width,height,i,j,padding);
                    if(first_iteration){
                        min = target[i*width+j];
                        max = target[i*width+j];
                        first_iteration = 0;
                        continue;
                    }

                    if(max < target[i*width+j]){
                        max = target[i*width+j];
                    }

                    if(min > target[i*width+j]){
                        min = target[i*width+j];
                    }  
                }  
            }
        }
    }
    if(method == SHARDED_COLUMNS_COLUMN_MAJOR){
        start = thread_num * size;
        end = (thread_num + 1) * size; //stops at the next thread work
        increment = 1;
        //since this is SHARDED_COLUMNS_COLUMN_MAJOR, we want to divide the columns and go through each row(0 to height)
        for(int j=start;j<end;j++){
            for (int i = 0 ; i < height ; i++){
                if ((i*width+j) < (width*height)){//for the threads that were assigned more than height and width
                    target[i*width+j] = applyfilter(f,inter,width,height,i,j,padding);
                    

                    if(first_iteration){
                        min = target[i*width+j];
                        max = target[i*width+j];
                        first_iteration = 0;
                        continue;
                    }

                    if(max < target[i*width+j]){
                        max = target[i*width+j];
                    }

                    if(min > target[i*width+j]){
                        min = target[i*width+j];
                    }
                }
            }
        }
    }


    if(method == SHARDED_COLUMNS_ROW_MAJOR){
        start = thread_num * size;
        end = (thread_num+1) * size;
        increment = 1;
        //since this is SHARDED_COLUMNS_ROW_MAJOR, we want to divide the columns and iterate through each row of the columns first
        //then we can move on to the next row
        for (int i = 0 ; i < height ; i++){
            for (int j = start ; j < end ; j += increment){
                if (j < width){ 
                    target[i*width+j] = applyfilter(f,inter,width,height,i,j,padding);
                    if(first_iteration){
                        min = target[i*width+j];
                        max = target[i*width+j];
                        first_iteration = 0;
                        continue;
                    }

                    if(max < target[i*width+j]){
                        max = target[i*width+j];
                    }

                    if(min > target[i*width+j]){
                        min = target[i*width+j];
                    }  
                }  
            }
        }  
    }
    if(method == WORK_QUEUE){
        queue_work(w);
        return NULL;
    }

    //assigning the final min and max to the struct
    pthread_mutex_lock(&w->common->lock);
    
    if(w->common->min > min){
        w->common->min = min;
    }

    if(w->common->max < max){
        w->common->max = max;
    }

    

    pthread_mutex_unlock(&w->common->lock);
    
    //After waiting for all of the threads, we are going to normalize the pixels
    //for normalization we are going to iterate the same depeding on the method and call normalize_pixel()
    pthread_barrier_wait(&w->common->barrier);
    min = w->common->min;
    max = w->common->max;
    if(method == SHARDED_ROWS){
        start = thread_num * size;
        end = (thread_num+1) * size;
        increment = 1;
        
        for (int i = start ; i < end ; i += increment){
            for(int j=0;j<width;j++){
                //i*width+j => turns row index and height index into one row index(for one dimensional matrix)
                if (i < height){
                    normalize_pixel(target, i*width+j, min, max);
                }
            }
        }
    }
    
    if(method == SHARDED_COLUMNS_COLUMN_MAJOR){
        start = thread_num * size;
        end = (thread_num + 1) * size;
        increment = 1;
        for(int j=start;j<end;j++){
            for (int i = 0 ; i < height ; i++){
                if ((i*width+j) < (width*height)){
                     //i*width+j => turns row index and height index into one row index(for one dimensional matrix)
                    normalize_pixel(target, i*width+j, min, max);
                }
            }
        }
    }

    if(method == SHARDED_COLUMNS_ROW_MAJOR){
        start = thread_num * size;
        end = (thread_num+1) * size;
        increment = 1;
        
        for (int i = 0 ; i < height ; i++){
            for (int j = start ; j < end ; j += increment){
                if (j < width){
                     //i*width+j => turns row index and height index into one row index(for one dimensional matrix)
                    normalize_pixel(target, i*width+j, min, max); 
                }  
            }
        }
    }

    return NULL;
}

/***************** MULTITHREADED ENTRY POINT ******/
/* TODO: this is where you should implement the multithreaded version
 * of the code. Use this function to identify which method is being used
 * and then call some other function that implements it.
 */
void apply_filter2d_threaded(const filter *f,
        const int32_t *original, int32_t *target,
        int32_t width, int32_t height,
        int32_t num_threads, parallel_method method, int32_t work_chunk)
{
    
    int padding=(int)(f->dimension/2);
    int32_t *inter=malloc(sizeof(int32_t)*(width+(padding*2))*(height+(padding*2)));
    common_work *common=malloc(sizeof(common_work));
    target = malloc(sizeof(int32_t)*(width)*(height));
    
    common->f=f;
    common->original_image=original;
    common->output_image=target;
    common->inter=inter;
    common->width=width;
    common->height=height;
    common->method=method;
    common->max_threads=num_threads;
    common->padding=padding;
    //depending on the methodology, our size for each thread will differ
    if(method == SHARDED_ROWS){
        if (num_threads>height){
            num_threads=height;
	    }
        common->size=(int)((height+num_threads-1)/num_threads);
    }
    else if(method == SHARDED_COLUMNS_COLUMN_MAJOR || method == SHARDED_COLUMNS_ROW_MAJOR){
        if (num_threads>width){
            num_threads=width;
        }
        common->size=(int)((width+num_threads-1)/num_threads);
    }else if (method == WORK_QUEUE){
        int numofchunks=(((width+work_chunk-1)/work_chunk)*((height+work_chunk-1)/work_chunk));
        if ((num_threads)>numofchunks){
            num_threads=numofchunks;
        }
    }
    
    common->max_threads=num_threads;
    common->work_chunk=work_chunk;
    pthread_barrier_init(&common->barrier,NULL,common->max_threads);
    pthread_mutex_init(&common->lastchunk, NULL);
    pthread_mutex_init(&common->lock, NULL);
    //making all of the elements in the target equal to 0
    for (int i=0;i<(width+padding+padding)*(height+padding+padding);i++){
        inter[i]=0;
    }
    //we have to check for the index that are going to be 0
    int index=0;//for keeping track of the index of the original during copying
    for (int i=0;i<(width+padding+padding)*(height+padding+padding);i+=1){
        if((i<(width+padding+padding)*padding)||i>=((width+padding+padding)*(height+padding))){
            continue;//takes care of top and bottom rows
        }
        if(i%(width+padding+padding)<padding){
            continue;//takes care of the left zero paddings
        }
        if(i%(width+padding+padding)>=padding+width){
            continue;//takes care of the right zero paddings
        }
        
        inter[i]=original[index];
        index++;
    }

    //initializing min and max
    common->min=applyfilter(f,inter,width,height,0,0,padding);
    common->max=applyfilter(f,inter,width,height,0,0,padding);

    pthread_t thread[num_threads];
	work thread_args[num_threads];
    
    
	for (int i=0; i<common->max_threads; i++){
		thread_args[i].id=i;
        thread_args[i].common=common;
		pthread_create(&thread[i],NULL, sharding_work,(void *) &thread_args[i]);	
	}

	for(int i=0;i<common->max_threads;i++){
		pthread_join(thread[i],NULL);
	}
    
    pthread_mutex_destroy(&common->lock);
    free(common);
    free(inter);
      
    pthread_mutex_destroy(&common->lastchunk);
}
