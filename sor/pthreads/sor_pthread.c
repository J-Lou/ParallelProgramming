/*
 * A Red-Black SOR
 *
 *	using separate red and black matrices
 *	to minimize false sharing
 *
 * Solves a M+2 by 2N+2 array
 *
 *      The overall matrix is separated to 
 *      one M+2 by N+1 red matrix, and one 
 *      M+2 by N+1 black marix. And the two
 *      matrices form the overall matrix in
 *      a checkerboard pattern fashion 
 *      as illustrated:
 *      
 *      -------------------------
 *      |R(1,0) |B(1,0) |R(1,1) |
 *      -------------------------
 *      |B(2,0) |R(2,0) |B(2,1) |
 *      -------------------------
 *      |R(3,0) |B(3,0) |R(3,1) |
 *     
 *      so that each array element's four
 *      immediate neighbors are of different
 *      color.
 *
 */
#ifndef _REENTRANT
#define _REENTRANT		/* basic 3-lines for threads */
#endif
#include <pthread.h>
//#include <thread.h>

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define RESET_AFTER_ONE_ITERATION 1
/*
 * Set this flag to 1 will measure the elasping time from 
 * the time the first iteration of the SOR is finished.
 * Otherwise, the time will start from the time SOR begins.
 *
 * The purpose for it is to exclude the time for loading
 * the data into the cache, warming up the TLB, predictors,
 * etc. So with first cycle as warmup, we will see more
 * stable speedup from parallelization.
 */

struct timeval start, finish;
int iterations = 100;
int M = 4000;
int N = 500;
int verify = 1;
/*
#define M	4000
#define	N	 254
*/

float **red_;
float **black_;
int task_num = 1;

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
void
barrier (int expect)
{
    static int arrived = 0;

    pthread_mutex_lock (&mut);	//lock

    arrived++;
    if (arrived < expect)
        pthread_cond_wait (&cond, &mut);
    else {
        arrived = 0;		// reset the barrier before broadcast is important
        pthread_cond_broadcast (&cond);
    }

    pthread_mutex_unlock (&mut);	//unlock
}

/*
 * 
 * This is the main working function for SOR.
 * In each iteration, the cell is updated by averaging
 * the north, south, east, west immdediate neighbors.
 *
 * With the way black and red matrices organized,
 * in each iteration, the black cells' new values are
 * first calculated by averaging nearby red cells.
 * Then the red cells are updated using the black
 * cells' value just calculated.
 *
 * sor_odd starts from the odd row.
 */
void
sor_odd (begin, end, task_id)
     int begin;
     int end;
     int task_id;
{
    int i, j, k;

#ifdef TIME_EACH_ITERATION
    struct timeval sstart, sfinish;
#endif

    for (i = 0; i < iterations; i++) {
#ifdef TIME_EACH_ITERATION
        gettimeofday (&sstart, NULL);
#endif

        for (j = begin; j <= end; j++) {

           /* update black cells on odd rows */

            for (k = 0; k < N; k++) {

	        black_[j][k] = (red_[j - 1][k] + red_[j + 1][k] + red_[j][k] + red_[j][k + 1]) / 4.0;
            }
            if ((j += 1) > end)
	        break;

            /* update black cells on even rows */

	    for (k = 1; k <= N; k++) {

	        black_[j][k] = (red_[j - 1][k] + red_[j + 1][k] + red_[j][k - 1] + red_[j][k]) / 4.0;
            }
        }

        barrier (task_num);

        for (j = begin; j <= end; j++) {

            /* update red cells on odd rows */

            for (k = 1; k <= N; k++) {

	        red_[j][k] = (black_[j - 1][k] + black_[j + 1][k] + black_[j][k - 1] + black_[j][k]) / 4.0;
            }
            if ((j += 1) > end)
	        break;

            /* update red cells on even rows */

            for (k = 0; k < N; k++) {

	        red_[j][k] = (black_[j - 1][k] + black_[j + 1][k] + black_[j][k] + black_[j][k + 1]) / 4.0;
            }
        }

        barrier (task_num);

#ifdef	RESET_AFTER_ONE_ITERATION
        if ((i == 0) && (task_id == 0)){
            puts ("restart");
            gettimeofday (&start, NULL);
        }
#endif

#ifdef TIME_EACH_ITERATION
        gettimeofday (&sfinish, NULL);
        printf ("One Iteration Elapsed time: %.2f seconds\n",
	        (((sfinish.tv_sec * 1000000.0) + sfinish.tv_usec) -
	        ((sstart.tv_sec * 1000000.0) + sstart.tv_usec)) / 1000000.0);
#endif
    }
}

/*
 * sor_even starts from the even row.
 */
void
sor_even (begin, end, task_id)
     int begin;
     int end;
     int task_id;
{
    int i, j, k;

    for (i = 0; i < iterations; i++) {

        for (j = begin; j <= end; j++) {

            for (k = 1; k <= N; k++) {

        	black_[j][k] = (red_[j - 1][k] + red_[j + 1][k] + red_[j][k - 1] + red_[j][k]) / 4.0;
            }
            if ((j += 1) > end)
	        break;

            for (k = 0; k < N; k++) {

	        black_[j][k] = (red_[j - 1][k] + red_[j + 1][k] + red_[j][k] + red_[j][k + 1]) / 4.0;
            }
        }
        barrier (task_num);

        for (j = begin; j <= end; j++) {

            for (k = 0; k < N; k++) {

	        red_[j][k] = (black_[j - 1][k] + black_[j + 1][k] + black_[j][k] + black_[j][k + 1]) / 4.0;
            }
            if ((j += 1) > end)
        	break;

            for (k = 1; k <= N; k++) {

	        red_[j][k] = (black_[j - 1][k] + black_[j + 1][k] + black_[j][k - 1] + black_[j][k]) / 4.0;
            }
        }
        barrier (task_num);
#ifdef	RESET_AFTER_ONE_ITERATION
        if ((i == 0) && (task_id == 0)){
            gettimeofday (&start, NULL);
        }
#endif
    }
}

extern char *optarg;

void
errexit (const char *err_str)
{
    fprintf (stderr, "%s", err_str);
    exit (1);
}

void *
work_thread (void *lp)
{
    int task_id = *((int *) lp);
    int begin, end;
    struct timeval start, finish;

    /*get the divided task*/
    begin = (M * task_id) / task_num + 1;
    end = (M * (task_id + 1)) / task_num;

    fprintf (stderr, "thread %d: begin %d, end %d\n", task_id, begin, end);

    barrier (task_num);

    if(task_id==0)
        gettimeofday (&start, NULL);

    /* do the SOR work */
    if (begin & 1)
        sor_odd (begin, end, task_id);
    else
        sor_even (begin, end, task_id);

    barrier (task_num);
    gettimeofday (&finish, NULL);

    if(task_id==0)
        printf ("Elapsed time: %.2f seconds\n",
	        (((finish.tv_sec * 1000000.0) + finish.tv_usec) -
	        ((start.tv_sec * 1000000.0) + start.tv_usec)) / 1000000.0);
}

int
main (argc, argv)
     int argc;
     char *argv[];
{
    int c, i, j;
    int begin, end;
    int iTotalSize;
    pthread_attr_t attr;
    pthread_t *tid;
    int *id;

    while ((c = getopt (argc, argv, "vi:m:n:p:")) != -1)
        switch (c) {
        case 'i':
            iterations = atoi (optarg);
            break;
        case 'm':
            M = atoi (optarg);
            break;
        case 'n':
            N = atoi (optarg);
            break;
        case 'v':
            verify = 1;
            break;
        case 'p':
            task_num = atoi (optarg);
            break;
        }

    printf ("%d tasks, N=%d, M=%d, iterations=%d\n", task_num, N, M, iterations);

    // initialization
    if ((red_ = (float **) malloc ((M + 2) * sizeof (float *))) == 0)
        errexit ("out of shared memory");

    if ((black_ = (float **) malloc ((M + 2) * sizeof (float *))) == 0)
        errexit ("out of shared memory");

    for (i = 0; i <= M + 1; i++) {

        if ((red_[i] = (float *) malloc ((N + 1) * sizeof (float))) == 0)
            errexit ("out of shared memory");

        if ((black_[i] = (float *) malloc ((N + 1) * sizeof (float))) == 0)
            errexit ("out of shared memory");
    }

    Initialize (red_, black_);

    // create threads
    id = (int *) malloc (sizeof (int) * task_num);
    tid = (pthread_t *) malloc (sizeof (pthread_t) * task_num);
    if (!id || !tid)
        errexit ("out of shared memory");
    pthread_attr_init (&attr);
    pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);
    for (i = 1; i < task_num; i++) {
        id[i] = i;
        pthread_create (&tid[i], &attr, work_thread, &id[i]);
    }

    id[0]=0;
    work_thread(&id[0]);
    // wait for all threads to finish
    for (i = 1; i < task_num; i++)
        pthread_join (tid[i], NULL);

    if (verify) {
        FILE *res;

        res = fopen ("vres", "w");
        for (i = 0; i < M + 2; i++) {
            if (i & 1)
	        for (j = 0; j < N + 1; j++) {
	            fprintf (res, "[%d][%d] = %f\n", i, 2 * j, red_[i][j]);
	            fprintf (res, "[%d][%d] = %f\n", i, 2 * j + 1, black_[i][j]);
	        }
           else
	       for (j = 0; j < N + 1; j++) {
	           fprintf (res, "[%d][%d] = %f\n", i, 2 * j, black_[i][j]);
	           fprintf (res, "[%d][%d] = %f\n", i, 2 * j + 1, red_[i][j]);
	       }
        }				/* for i */
    }

    return 0;
}

/***************************************************************************\
        Initialize() intializes the array. 
        
        The black and red array are given some non-zero values as the initial
        input for the SOR test.
 
        Right now the initial values in each row repeat in an interval of 35
        and the value is calculated in an increasing exponential fashion.
        The goal here is merely to initialize some non-zero float point value,
        and the value will be calculated every iteration in SOR. So there is
        no practical meaning for this input which you need to worry about. 
        Simply check if the result is the same via output file after
        parallelized.
\***************************************************************************/

Initialize (red, black)
     float **red, **black;
{
    extern double exp ();
 
    int j, k;
 
    for (j = 0; j < M + 2; j++) {
 
            for (k = 0; k <= N; k++) {
                black[j][k] =(float) exp ((double) (k%35));
                red[j][k] = (float) exp (k%35 + (double) 3);
            }
    }

}
