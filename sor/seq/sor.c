/*
 * A Red-Black SOR
 *
 *      using separate red and black matrices
 *      to minimize false sharing
 *
 * Solves a M+2 by 2N+2 array
 *
 *	The overall matrix is separated to 
 *	one M+2 by N+1 red matrix, and one 
 *	M+2 by N+1 black marix. And the two
 *	matrices form the overall matrix in
 *	a checkerboard pattern fashion 
 *	as illustrated:
 *
 * 	-------------------------
 *	|R(1,0)	|B(1,0)	|R(1,1)	|
 *	-------------------------
 *	|B(2,0)	|R(2,0)	|B(2,1)	|
 *	-------------------------
 *	|R(3,0)	|B(3,0)	|R(3,1)	|
 *
 *	so that each array element's four
 *	immediate neighbors are of different
 *	color.
 *
 */
#include <stdio.h>
#include <stdlib.h>
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
int verify = 1;	/*print the result to file if set to 1 */

float **red_;
float **black_;

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
sor_odd (begin, end)
     int begin;
     int end;
{
    int i, j, k;

    for (i = 0; i < iterations; i++) {

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

#ifdef	RESET_AFTER_ONE_ITERATION
	if (i == 0)
	    gettimeofday (&start, NULL);
#endif
    }
}

extern char *optarg;

main (argc, argv)
     int argc;
     char *argv[];
{
    int c, i, j;

    while ((c = getopt (argc, argv, "vi:m:n:")) != -1)
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
	}

    if ((red_ = (float **) malloc ((M + 2) * sizeof (float *))) == 0) {
	perror ("out of shared memory");
	exit (-1);
    }

    if ((black_ = (float **) malloc ((M + 2) * sizeof (float *))) == 0) {
	perror ("out of shared memory");
	exit (-1);
    }

    for (i = 0; i <= M + 1; i++) {

	if ((red_[i] = (float *) malloc ((N + 1) * sizeof (float))) == 0) {
	    perror ("out of shared memory");
	    exit (-1);
	}

	if ((black_[i] = (float *) malloc ((N + 1) * sizeof (float))) == 0) {
	    perror ("out of shared memory");
	    exit (-1);
	}
    }

    /* initialize the matrix */
    Initialize (red_, black_);

    gettimeofday (&start, NULL);

    /* do the SOR work */
    sor_odd (1, M);

    gettimeofday (&finish, NULL);

    printf ("Elapsed time: %.2f seconds\n",
	    (((finish.tv_sec * 1000000.0) + finish.tv_usec) -
	     ((start.tv_sec * 1000000.0) + start.tv_usec)) / 1000000.0);

    /* print the full matrix out in normal order if verify is set */

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
	}			/* for i */
    }
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
