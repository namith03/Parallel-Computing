//
// Sorts a list using multiple threads
//

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#define MAX_THREADS     65536
#define MAX_LIST_SIZE   270000000

#define DEBUG 0

// Thread variables
//
// VS: ... declare thread variables, mutexes, condition varables, etc.,
// VS: ... as needed for this assignment 
//

typedef struct node_barrier {
    pthread_mutex_t cnt_lck;
    pthread_cond_t sync;
    pthread_cond_t release;
    int count;
} node_barrier;

node_barrier bnode[MAX_THREADS];
pthread_t p_threads[MAX_THREADS];

// routine to initiliaze barrier nodes
void mylog_init_barrier(node_barrier *b) {
    for(int i = 0; i < MAX_THREADS; ++i) {
        b[i].count = 0;
        pthread_mutex_init(&(b[i].cnt_lck), NULL);
        pthread_cond_init(&(b[i].sync), NULL);
        pthread_cond_init(&(b[i].release), NULL);
    }
}

int *ptr, *thread_ids;
int np, q;

// Global variables
int num_threads;        // Number of threads to create - user input 
int list_size;          // List size
int *list;          // List of values
int *work;          // Work array
int *list_orig;         // Original list of values, used for error checking

// Print list - for debugging
void print_list(int *list, int list_size) {
    int i;
    for (i = 0; i < list_size; i++) {
        printf("[%d] \t %16d\n", i, list[i]); 
    }
    printf("--------------------------------------------------------------------\n"); 
}

// Comparison routine for qsort (stdlib.h) which is used to 
// a thread's sub-list at the start of the algorithm
int compare_int(const void *a0, const void *b0) {
    int a = *(int *)a0;
    int b = *(int *)b0;
    if (a < b) {
        return -1;
    } else if (a > b) {
        return 1;
    } else {
        return 0;
    }
}

// Return index of first element larger than or equal to v in sorted list
// ... return last if all elements are smaller than v
// ... elements in list[first], list[first+1], ... list[last-1]
//
//   int idx = first; while ((v > list[idx]) && (idx < last)) idx++;
//
int binary_search_lt(int v, int *list, int first, int last) {
   
    // Linear search code
    // int idx = first; while ((v > list[idx]) && (idx < last)) idx++; return idx;

    int left = first; 
    int right = last-1; 

    if (list[left] >= v) {
        if(DEBUG)   printf("val = %d, first = %d, last = %d, ret = %d\n", v, first, last, left);
        return left;
    }
    if (list[right] < v) {
        if(DEBUG)   printf("val = %d, first = %d, last = %d, ret = %d\n", v, first, last, right+1);        
        return right+1;
    }
    int mid = (left+right)/2; 
    while (mid > left) {
        if (list[mid] < v) {
        left = mid; 
    } else {
        right = mid;
    }
    mid = (left+right)/2;
    }
    if(DEBUG)   printf("val = %d, first = %d, last = %d, ret = %d\n", v, first, last, right);        
    return right;
}
// Return index of first element larger than v in sorted list
// ... return last if all elements are smaller than or equal to v
// ... elements in list[first], list[first+1], ... list[last-1]
//
//   int idx = first; while ((v >= list[idx]) && (idx < last)) idx++;
//
int binary_search_le(int v, int *list, int first, int last) {

    // Linear search code
    // int idx = first; while ((v >= list[idx]) && (idx < last)) idx++; return idx;
 
    int left = first; 
    int right = last-1; 

    if (list[left] > v) return left; 
    if (list[right] <= v) return right+1;
    int mid = (left+right)/2; 
    while (mid > left) {
        if (list[mid] <= v) {
        left = mid; 
    } else {
        right = mid;
    }
    mid = (left+right)/2;
    }
    return right;
}

void *parallel_sort(void *temp) {
    int my_id = *(int*)temp;
    if(DEBUG) {
        printf("my_id = %d\n", my_id);
        printf("q = %d\n", q);
    }

    int i, level; 
    int my_list_size;

    int my_own_blk, my_own_idx;
    int my_blk_size, my_search_blk, my_search_idx, my_search_idx_max;
    int my_write_blk, my_write_idx;
    int my_search_count; 
    int idx, i_write; 

    int thread_grp_sz, base, div;
    int thread_id, index;

    my_list_size = ptr[my_id+1]-ptr[my_id];
    if(DEBUG)   printf("start = %d, end = %d, my_list_size = %d\n", ptr[my_id], ptr[my_id+1], my_list_size);
    // printf("Before initial sort, in thread %d :", my_id);
    // for(int i = ptr[my_id]; i < ptr[my_id+1]; ++i) {
    //     printf("before thread_id = %d, id = %d, val = %d \n", my_id, i, list[i]);
    // }
    qsort(&list[ptr[my_id]], my_list_size, sizeof(int), compare_int);
    // printf("After initial sort, in thread %d :", my_id);
    // for(int i = ptr[my_id]; i < ptr[my_id+1]; ++i) {
    //     printf("after thread_id = %d, id = %d, val = %d \n", my_id, i, list[i]);
    // }

    
    //TODO: add barrier logic, before we start merging sorted sub-lists
    level = 0;
    // thread_grp_sz = (1 << (level+1)); 
    thread_grp_sz = num_threads;
    base = (my_id / thread_grp_sz) * thread_grp_sz; div = 2;
    // thread_id = (my_id - base);
    thread_id = my_id;
    if(num_threads != 1) {
        do {
            index = base + thread_id / div;
            if(DEBUG)   printf("initial barrier thread_id = %d, thread_pos = %d, index = %d, thread_grp_sz = %d\n", my_id, thread_id, index, thread_grp_sz);
            if(thread_id%div == 0) {
                pthread_mutex_lock(&(bnode[index].cnt_lck));
                bnode[index].count++;
                while(bnode[index].count < 2) {
                    pthread_cond_wait(&(bnode[index].sync), &(bnode[index].cnt_lck));
                }
                pthread_mutex_unlock(&(bnode[index].cnt_lck));
            } else {
                pthread_mutex_lock(&(bnode[index].cnt_lck));
                bnode[index].count++;
                if(bnode[index].count == 2) {
                    pthread_cond_signal(&(bnode[index].sync));
                }
                while(pthread_cond_wait(&(bnode[index].release), &(bnode[index].cnt_lck)) != 0);
                pthread_mutex_unlock(&(bnode[index].cnt_lck));
                break;
            }
            base += thread_grp_sz / div, div *= 2;
        } while(div <= thread_grp_sz);

        div /= 2;
        while(div > 1) {
            base -= thread_grp_sz / div;
            // thread_id = (my_id - base);
            index = base + thread_id / div;
            pthread_mutex_lock(&(bnode[index].cnt_lck));
            bnode[index].count = 0;
            if(DEBUG)   printf("level = %d, thread_id = %d, releasing barrier node = %d\n", level, my_id, index);
            pthread_cond_signal(&(bnode[index].release));
            pthread_mutex_unlock(&(bnode[index].cnt_lck));
            div /= 2;
        }
        if(DEBUG)   printf("thread_id = %d exited initial barrier = %d\n", my_id, level);
    }

    // int s = my_id == 0 ? 8 : 0;
    // int e = my_id == 0 ? 16 : 8;
    // for(int i = s; i < e; ++i) {
    //     printf("thread_id = %d, idx = %d, search val = %d\n", my_id, i, list[i]);
    // }

    for(int level = 0; level < q; ++level) {
        my_blk_size = np * (1 << level); 

        my_own_blk = ((my_id >> level) << level);
        my_own_idx = ptr[my_own_blk];

        my_search_blk = ((my_id >> level) << level) ^ (1 << level);
        my_search_idx = ptr[my_search_blk];
        my_search_idx_max = my_search_idx+my_blk_size;

        my_write_blk = ((my_id >> (level+1)) << (level+1));
        my_write_idx = ptr[my_write_blk];

        idx = my_search_idx;
        
        my_search_count = 0;

        if(DEBUG) {
            printf("level = %d, thread_id = %d, my_search_idx = %d, my_search_idx_max = %d, my_search_blk = %d, my_write_idx = %d, my_write_blk = %d, my_own_idx = %d, my_own_blk = %d\n", level, my_id, my_search_idx, my_search_idx_max, my_search_blk, my_write_idx, my_write_blk, my_own_idx, my_own_blk);
        }

        // Binary search for 1st element
        if (my_search_blk > my_own_blk) {
               idx = binary_search_lt(list[ptr[my_id]], list, my_search_idx, my_search_idx_max); 
        } else {
               idx = binary_search_le(list[ptr[my_id]], list, my_search_idx, my_search_idx_max); 
        }
        my_search_count = idx - my_search_idx;
        i_write = my_write_idx + my_search_count + (ptr[my_id]-my_own_idx); 
        if(DEBUG)   printf("level = %d, thread_id = %d, first idx = %d, first id = %d, val = %d, first write id = %d, my_search_count = %d\n", level, my_id, idx, ptr[my_id], list[ptr[my_id]], i_write, my_search_count);
        work[i_write] = list[ptr[my_id]];

        // Linear search for 2nd element onwards
        for (i = ptr[my_id]+1; i < ptr[my_id+1]; i++) {
            if (my_search_blk > my_own_blk) {
                while ((list[i] > list[idx]) && (idx < my_search_idx_max)) {
                    idx++; my_search_count++;
                }
            } else {
                while ((list[i] >= list[idx]) && (idx < my_search_idx_max)) {
                    idx++; my_search_count++;
                }
            }
            i_write = my_write_idx + my_search_count + (i-my_own_idx); 
            if(DEBUG)   printf("level = %d, thread_id = %d, id = %d, val = %d, write id = %d, my_search_count = %d\n", level, my_id, i, list[i], i_write, my_search_count);
            work[i_write] = list[i];
        }        

        //TODO: add barrier logic, before copying the work to list
        // thread_grp_sz = (1 << (level+1 + (level+1 == q ? 0 : 1)));
        thread_grp_sz = num_threads;
        base = (my_id / thread_grp_sz) * thread_grp_sz;
        div = 2;
        // thread_id = (my_id - base);
        thread_id = my_id;

        if(num_threads != 1) {
            do {
                index = base + thread_id / div;
                if(DEBUG)   printf("level = %d, thread_id = %d, thread_pos = %d, index = %d, base = %d, thread_grp_sz = %d\n", level, my_id, thread_id, index, base, thread_grp_sz);
                if(thread_id%div == 0) {
                    pthread_mutex_lock(&(bnode[index].cnt_lck));
                    bnode[index].count++;
                    while(bnode[index].count < 2) {
                        pthread_cond_wait(&(bnode[index].sync), &(bnode[index].cnt_lck));
                    }
                    pthread_mutex_unlock(&(bnode[index].cnt_lck));
                } else {
                    pthread_mutex_lock(&(bnode[index].cnt_lck));
                    bnode[index].count++;
                    if(bnode[index].count == 2) {
                        pthread_cond_signal(&(bnode[index].sync));
                    }
                    if(DEBUG)   printf("level = %d, thread_id = %d waiting...\n", level, my_id);
                    while(pthread_cond_wait(&(bnode[index].release), &(bnode[index].cnt_lck)) != 0);
                    pthread_mutex_unlock(&(bnode[index].cnt_lck));
                    break;
                }
                base += thread_grp_sz / div, div *= 2;
            } while(div <= thread_grp_sz);

            div /= 2;
            while(div > 1) {
                base -= thread_grp_sz / div;
                // thread_id = (my_id - base);
                index = base + thread_id / div;
                pthread_mutex_lock(&(bnode[index].cnt_lck));
                bnode[index].count = 0;
                if(DEBUG)   printf("level = %d, thread_id = %d, releasing barrier node = %d\n", level, my_id, index);
                pthread_cond_signal(&(bnode[index].release));
                pthread_mutex_unlock(&(bnode[index].cnt_lck));
                div /= 2;
            }
            if(DEBUG)   printf("thread_id = %d exited level = %d\n", my_id, level);
        }

        for (i = ptr[my_id]; i < ptr[my_id+1]; i++) {
            if(DEBUG)   printf("level = %d, thread_id = %d, id = %d, work[id] = %d\n", level, my_id, i, work[i]);
            list[i] = work[i];
        } 

        //TODO: add barrier logic, before copying the work to list
        // thread_grp_sz = (1 << (level+1 + (level+1 == q ? 0 : 1)));
        thread_grp_sz = num_threads;
        base = (my_id / thread_grp_sz) * thread_grp_sz;
        div = 2;
        // thread_id = (my_id - base);
        thread_id = my_id;
        if(num_threads != 1) {
            do {
                index = base + thread_id / div;
                if(DEBUG)   printf("level = %d, thread_id = %d, thread_pos = %d, index = %d, base = %d, thread_grp_sz = %d\n", level, my_id, thread_id, index, base, thread_grp_sz);
                if(thread_id%div == 0) {
                    pthread_mutex_lock(&(bnode[index].cnt_lck));
                    bnode[index].count++;
                    while(bnode[index].count < 2) {
                        pthread_cond_wait(&(bnode[index].sync), &(bnode[index].cnt_lck));
                    }
                    pthread_mutex_unlock(&(bnode[index].cnt_lck));
                } else {
                    pthread_mutex_lock(&(bnode[index].cnt_lck));
                    bnode[index].count++;
                    if(bnode[index].count == 2) {
                        pthread_cond_signal(&(bnode[index].sync));
                    }
                    if(DEBUG)   printf("level = %d, thread_id = %d waiting...\n", level, my_id);
                    while(pthread_cond_wait(&(bnode[index].release), &(bnode[index].cnt_lck)) != 0);
                    pthread_mutex_unlock(&(bnode[index].cnt_lck));
                    break;
                }
                base += thread_grp_sz / div, div *= 2;
            } while(div <= thread_grp_sz);

            div /= 2;
            while(div > 1) {
                base -= thread_grp_sz / div;
                // thread_id = (my_id - base);
                index = base + thread_id / div;
                pthread_mutex_lock(&(bnode[index].cnt_lck));
                bnode[index].count = 0;
                if(DEBUG)   printf("level = %d, thread_id = %d, releasing barrier node = %d\n", level, my_id, index);
                pthread_cond_signal(&(bnode[index].release));
                pthread_mutex_unlock(&(bnode[index].cnt_lck));
                div /= 2;
            }
            if(DEBUG)   printf("thread_id = %d exited level = %d\n", my_id, level);
        }

    }
    pthread_exit(NULL);
}

// Sort list via parallel merge sort
//
// VS: ... to be parallelized using threads ...
//
void sort_list(int q) {

    int i, level, my_id; 
    int np, my_list_size; 
    int ptr[num_threads+1]; //TODO: add in main

    int my_own_blk, my_own_idx;
    int my_blk_size, my_search_blk, my_search_idx, my_search_idx_max;
    int my_write_blk, my_write_idx;
    int my_search_count; 
    int idx, i_write; 
    
    //TODO: add in main
    np = list_size/num_threads;     // Sub list size 

    //TODO: add in main
    // Initialize starting position for each sublist
    for (my_id = 0; my_id < num_threads; my_id++) {
        ptr[my_id] = my_id * np;
    }
    ptr[num_threads] = list_size;

    // Sort local lists
    for (my_id = 0; my_id < num_threads; my_id++) {
        my_list_size = ptr[my_id+1]-ptr[my_id];
        qsort(&list[ptr[my_id]], my_list_size, sizeof(int), compare_int);
    }
if (DEBUG) print_list(list, list_size); 

    // Sort list
    for (level = 0; level < q; level++) {

        // Each thread scatters its sub_list into work array
    for (my_id = 0; my_id < num_threads; my_id++) {

        my_blk_size = np * (1 << level); 

        my_own_blk = ((my_id >> level) << level);
        my_own_idx = ptr[my_own_blk];

        my_search_blk = ((my_id >> level) << level) ^ (1 << level);
        my_search_idx = ptr[my_search_blk];
        my_search_idx_max = my_search_idx+my_blk_size;

        my_write_blk = ((my_id >> (level+1)) << (level+1));
        my_write_idx = ptr[my_write_blk];

        idx = my_search_idx;
        
        my_search_count = 0;


        // Binary search for 1st element
        if (my_search_blk > my_own_blk) {
               idx = binary_search_lt(list[ptr[my_id]], list, my_search_idx, my_search_idx_max); 
        } else {
               idx = binary_search_le(list[ptr[my_id]], list, my_search_idx, my_search_idx_max); 
        }
        my_search_count = idx - my_search_idx;
        i_write = my_write_idx + my_search_count + (ptr[my_id]-my_own_idx); 
        work[i_write] = list[ptr[my_id]];

        // Linear search for 2nd element onwards
        for (i = ptr[my_id]+1; i < ptr[my_id+1]; i++) {
            if (my_search_blk > my_own_blk) {
            while ((list[i] > list[idx]) && (idx < my_search_idx_max)) {
                idx++; my_search_count++;
            }
        } else {
            while ((list[i] >= list[idx]) && (idx < my_search_idx_max)) {
                idx++; my_search_count++;
            }
        }
        i_write = my_write_idx + my_search_count + (i-my_own_idx); 
        work[i_write] = list[i];
        }
    }
        // Copy work into list for next itertion
    for (my_id = 0; my_id < num_threads; my_id++) {
        for (i = ptr[my_id]; i < ptr[my_id+1]; i++) {
            list[i] = work[i];
        } 
    }
if (DEBUG) print_list(list, list_size); 
    }
}

// Main program - set up list of random integers and use threads to sort the list
//
// Input: 
//  k = log_2(list size), therefore list_size = 2^k
//  q = log_2(num_threads), therefore num_threads = 2^q
//
int main(int argc, char *argv[]) {

    struct timespec start, stop, stop_qsort;
    double total_time, time_res, total_time_qsort;
    int k, j, error; 

    // Read input, validate
    if (argc != 3) {
    printf("Need two integers as input \n"); 
    printf("Use: <executable_name> <log_2(list_size)> <log_2(num_threads)>\n"); 
    exit(0);
    }
    k = atoi(argv[argc-2]);
    if ((list_size = (1 << k)) > MAX_LIST_SIZE) {
    printf("Maximum list size allowed: %d.\n", MAX_LIST_SIZE);
    exit(0);
    }; 
    q = atoi(argv[argc-1]);
    if ((num_threads = (1 << q)) > MAX_THREADS) {
    printf("Maximum number of threads allowed: %d.\n", MAX_THREADS);
    exit(0);
    }; 
    if (num_threads > list_size) {
    printf("Number of threads (%d) < list_size (%d) not allowed.\n", 
       num_threads, list_size);
    exit(0);
    }; 

    if(DEBUG) {
        printf("num_threads = %d\n", num_threads);
        printf("list_size = %d\n", list_size);
    }

    // Allocate list, list_orig, and work

    list = (int *) malloc(list_size * sizeof(int));
    list_orig = (int *) malloc(list_size * sizeof(int));
    work = (int *) malloc(list_size * sizeof(int));

//
// VS: ... May need to initialize mutexes, condition variables, 
// VS: ... and their attributes
//

    ptr = (int *) malloc((num_threads+1) * sizeof(int));
    thread_ids = (int *) malloc((num_threads+1) * sizeof(int));
    np = list_size / num_threads;

    for(int thread_id = 0; thread_id < num_threads; ++thread_id) {
        ptr[thread_id] = thread_id * np;
        thread_ids[thread_id] = thread_id;
    }
    ptr[num_threads] = list_size;

    mylog_init_barrier(bnode);

    // Initialize list of random integers; list will be sorted by 
    // multi-threaded parallel merge sort
    // Copy list to list_orig; list_orig will be sorted by qsort and used
    // to check correctness of multi-threaded parallel merge sort
    srand48(0);     // seed the random number generator
    if(DEBUG)   printf("Elements in the list: ");
    for (j = 0; j < list_size; j++) {
    // list[j] = (int) lrand48();
    list[j] = list_size - j;
    if(DEBUG)   printf("%d ", list[j]);
    list_orig[j] = list[j];
    }
    if(DEBUG)   printf("\n");
    // duplicate first value at last location to test for repeated values
    //TODO: remove this comment
    // list[list_size-1] = list[0]; list_orig[list_size-1] = list_orig[0];

    // Create threads; each thread executes find_minimum
    clock_gettime(CLOCK_REALTIME, &start);

//
// VS: ... may need to initialize mutexes, condition variables, and their attributes
//

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

// Serial merge sort 
// VS: ... replace this call with multi-threaded parallel routine for merge sort
// VS: ... need to create threads and execute thread routine that implements 
// VS: ... parallel merge sort

    // sort_list(q);

    int status;
    for(int i = 0; i < num_threads; ++i) {
        if(DEBUG)   printf("creating thread = %d\n", i);
        status = pthread_create(&p_threads[i], &attr, parallel_sort, (void *)&thread_ids[i]);
        if(status != 0) {
            printf("Non-zero status creating thread \n");
        }
    }
    pthread_attr_destroy(&attr);
    for(int i = 0; i < num_threads; ++i) {
        pthread_join(p_threads[i], NULL);
    }
    if(DEBUG) {
        printf("After parallel sort: ");
        for(int i = 0; i < list_size; ++i) {
            printf("%d ", list[i]);
        }
        printf("\n");
    }

    // Compute time taken
    clock_gettime(CLOCK_REALTIME, &stop);
    total_time = (stop.tv_sec-start.tv_sec)
    +0.000000001*(stop.tv_nsec-start.tv_nsec);

    // Check answer
    qsort(list_orig, list_size, sizeof(int), compare_int);
    clock_gettime(CLOCK_REALTIME, &stop_qsort);
    total_time_qsort = (stop_qsort.tv_sec-stop.tv_sec)
    +0.000000001*(stop_qsort.tv_nsec-stop.tv_nsec);

    error = 0; 
    for (j = 1; j < list_size; j++) {
    if (list[j] != list_orig[j]) error = 1; 
    }

    if (error != 0) {
    printf("Houston, we have a problem!\n"); 
    }

    // Print time taken
    printf("List Size = %d, Threads = %d, error = %d, time (sec) = %8.4f, qsort_time = %8.4f\n", 
        list_size, num_threads, error, total_time, total_time_qsort);

// VS: ... destroy mutex, condition variables, etc.

    free(list); free(work); free(list_orig); 

}

