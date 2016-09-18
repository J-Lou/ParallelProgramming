/* 
 * Original author:  UNKNOWN
 *
 * Modified:         Kai Shen (January 2010)
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

//#define DEBUG 

#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}

/* Solve the equation:
 *   matrix * X = R
 */

double **matrix, *X, *R;

/* Pre-set solution. */

double *X__;

/* Number of threads */

int task_num = 1;

/* Number of rows */

int nsize;

/* The current row */

int currow;

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

/** Barrier that stops all threads until each thread has reached it. */

void barrier (int expect)
{
    static int arrived = 0;

    pthread_mutex_lock (&mut);

    arrived++;
    if (arrived < expect)
        pthread_cond_wait (&cond, &mut);
    else {
        arrived = 0;		// reset the barrier before broadcast is important
        pthread_cond_broadcast (&cond);
    }

    pthread_mutex_unlock (&mut);
}

/* Initialize the Matrix */

int initMatrix(const char *fname)
{
    FILE *file;
    int l1, l2, l3;
    double d;
    int i, j;
    double *tmp;
    char buffer[1024];

    if ((file = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "The matrix file open error\n");
        exit(-1);
    }
    
    /* Parse the first line to get the matrix size. */
    fgets(buffer, 1024, file);
    sscanf(buffer, "%d %d %d", &l1, &l2, &l3);
    nsize = l1;
#ifdef DEBUG
    fprintf(stdout, "matrix size is %d\n", nsize);
#endif

    /* Initialize the space and set all elements to zero. */
    matrix = (double**)malloc(nsize*sizeof(double*));
    assert(matrix != NULL);
    tmp = (double*)malloc(nsize*nsize*sizeof(double));
    assert(tmp != NULL);    
    for (i = 0; i < nsize; i++) {
        matrix[i] = tmp;
        tmp = tmp + nsize;
    }
    for (i = 0; i < nsize; i++) {
        for (j = 0; j < nsize; j++) {
            matrix[i][j] = 0.0;
        }
    }

    /* Parse the rest of the input file to fill the matrix. */
    for (;;) {
	fgets(buffer, 1024, file);
	sscanf(buffer, "%d %d %lf", &l1, &l2, &d);
	if (l1 == 0) break;

	matrix[l1-1][l2-1] = d;
#ifdef DEBUG
	fprintf(stdout, "row %d column %d of matrix is %e\n", l1-1, l2-1, matrix[l1-1][l2-1]);
#endif
    }

    fclose(file);
    return nsize;
}

/* Initialize the right-hand-side following the pre-set solution. */

void initRHS()
{
    int i, j;

    X__ = (double*)malloc(nsize * sizeof(double));
    assert(X__ != NULL);
    for (i = 0; i < nsize; i++) {
	    X__[i] = i+1;
    }

    R = (double*)malloc(nsize * sizeof(double));
    assert(R != NULL);
    for (i = 0; i < nsize; i++) {
	    R[i] = 0.0;
	    for (j = 0; j < nsize; j++) {
	        R[i] += matrix[i][j] * X__[j];
	    }
    }
}

/* Initialize the results. */

void initResult()
{
    int i;

    X = (double*)malloc(nsize * sizeof(double));
    assert(X != NULL);
    for (i = 0; i < nsize; i++) {
	    X[i] = 0.0;
    }
}

/* Get the pivot - the element on column with largest absolute value. */

void getPivot()
{
    int i, pivotrow;

    pivotrow = currow;
    for (i = currow+1; i < nsize; i++) {
	    if (fabs(matrix[i][currow]) > fabs(matrix[pivotrow][currow])) {
	        pivotrow = i;
	    }
    }

    if (fabs(matrix[pivotrow][currow]) == 0.0) {
        fprintf(stderr, "The matrix is singular\n");
        exit(-1);
    }
    
    if (pivotrow != currow) {
#ifdef DEBUG
	fprintf(stdout, "pivot row at step %5d is %5d\n", currow, pivotrow);
#endif
        for (i = currow; i < nsize; i++) {
            SWAP(matrix[pivotrow][i],matrix[currow][i]);
        }
        SWAP(R[pivotrow],R[currow]);
    }
}

void errexit(const char *err_str)
{
    fprintf (stderr, "%s", err_str);
    exit(1);
}

void * workerThread(void *lp)
{ 
    int task_id = *((int *) lp);
    int j, k;
    double pivotval = matrix[currow][currow];
    
    if (task_id == 0) {

	    getPivot();
        
	    /* Scale the main row. */
	    if (pivotval != 1.0) {
	        matrix[currow][currow] = 1.0;
	        for (j = currow + 1; j < nsize; j++) {
	    	    matrix[currow][j] /= pivotval;
	        }
	        R[currow] /= pivotval;
	    }
    }
        
    barrier(task_num);

    int begin = currow + task_id;

	/* Factorize the rest of the matrix. */
    for (j = begin + 1; j < nsize; j = j + task_num) {
        pivotval = matrix[j][begin];
        matrix[j][begin] = 0.0;
        for (k = begin + 1; k < nsize; k++) {
            matrix[j][k] -= pivotval * matrix[begin][k];
        }
        R[j] -= pivotval * R[begin];
    }

    return NULL;
}

void initThreads()
{
    int i;
    int *id = (int *) malloc (sizeof (int) * task_num);
    pthread_t *tid = (pthread_t *) malloc (sizeof (pthread_t) * task_num);
    pthread_attr_t attr;

    if (!id || !tid)
        errexit ("out of shared memory");

    pthread_attr_init (&attr);
    pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);
    for (i = 0; i < task_num; i++) {
        id[i] = i;
        pthread_create (&tid[i], &attr, workerThread, &id[i]);
    } 
    
    for (i = 0; i < task_num; i++)
        pthread_join (tid[i], NULL);
}

/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */

void computeGauss()
{
    int i;
    for (i = 0; i < nsize; i++) {
        currow=i;
        initThreads();
    }
}

/* Solve the equation. */

void solveGauss()
{
    int i, j;

    X[nsize-1] = R[nsize-1];
    for (i = nsize - 2; i >= 0; i --) {
        X[i] = R[i];
        for (j = nsize - 1; j > i; j--) {
            X[i] -= matrix[i][j] * X[j];
        }
    }

#ifdef DEBUG
    fprintf(stdout, "X = [");
    for (i = 0; i < nsize; i++) {
        fprintf(stdout, "%.6f ", X[i]);
    }
    fprintf(stdout, "];\n");
#endif
}

int main(int argc, char *argv[])
{
    int i;
    struct timeval start, finish;
    double error;
    
    if (argc != 3) {
	    fprintf(stderr, "usage: %s <matrixfile> <# threads>\n", argv[0]);
	    exit(-1);
    }

    task_num = atoi(argv[2]);
    initMatrix(argv[1]);
    initRHS();
    initResult();

    gettimeofday(&start, 0);
    computeGauss();
    gettimeofday(&finish, 0);

    solveGauss(nsize);
    
    fprintf(stdout, "Time:  %f seconds\n", (finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec)*0.000001);

    error = 0.0;
    for (i = 0; i < nsize; i++) {
	    double error__ = (X__[i]==0.0) ? 1.0 : fabs((X[i]-X__[i])/X__[i]);
	    if (error < error__) {
	        error = error__;
	    }
    }
    fprintf(stdout, "Error: %e\n", error);

    return 0;
}
