/*
 * A Red-Black SOR
 *
 *	using separate red and block matrices
 *	to minimize false sharing
 *
 * Solves a M+2 by 2N+2 array
 *
 * The overall matrix is separated to 
 * one M+2 by N+1 red matrix, and one 
 * M+2 by N+1 black marix. And the two
 * matrices form the overall matrix in
 * a checkerboard pattern fashion 
 * as illustrated:
 * 
 * -------------------------
 * |R(1,0) |B(1,0) |R(1,1) |
 * -------------------------
 * |B(2,0) |R(2,0) |B(2,1) |
 * -------------------------
 * |R(3,0) |B(3,0) |R(3,1) |
 * 
 * so that each array element's four
 * immediate neighbors are of different
 * color.
 * 
 */

/* Ported to MPI by
Cristiana Seidel 
(seidel@cos.ufrj.br) */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#define MAXPROCS 28
/*#define STA */
//define STA to collect msg stats

int statag;
#ifdef STA
int vsta[2];
int _data_send = 0;
int _msg_send = 0;
#endif

MPI_Status status; 
MPI_Status status1;
MPI_Status status2;
MPI_Request request;
MPI_Request request1;
MPI_Request request2;
int flag1, flag2;

struct	timeval	start, finish;
void sor_odd(int begin, int end);
void sor_even(int begin, int end);
void ex_black(int length);
void ex_red(int length);
void init_data(int begin, int end);
int  position;

int verify = 0;
int	iterations = 100;
int	M = 4000;
int	N = 500;	/* N.B. There are 2N columns. */
int proc_id, nprocs, mytid, tids[MAXPROCS];

float **red_;
float **black_;

extern	char	       *optarg;

/*
 * To begin statistics collection after the first iteration
 * has completed use the following #define
 *
 */
#define	RESET_AFTER_ONE_ITERATION

main(argc, argv)
     int		argc;
     char	       *argv[];
{
  int		c, i, j, NHosts;
  int		begin, end;
  int d;
  

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

  while ((c = getopt(argc, argv, "vi:m:n:")) != -1)
    switch (c) {
    case 'i':
      iterations = atoi(optarg);
      break;
    case 'm':
      M = atoi(optarg);
      break;
    case 'n':
      N = atoi(optarg);
      break;
    case 'v': verify = 1; break;
    }


  fprintf(stderr, " RED-BLACK SOR... \n");
  fprintf(stderr, "nprocs = %d\n", nprocs);
  fprintf(stderr, "%d iterations, M = %d, N = %d, %d X %d \n", iterations, M, N, M+2, 2*(N+1));
  fflush(stderr);

  begin = (M*proc_id)/nprocs + 1;
  end   = (M*(proc_id + 1))/nprocs;

     statag = 11;

#ifdef TEST_CH
    /* everybody sends to everybody else */
    for (i = 0; i < nprocs; i++)
      if (i != proc_id) {
        MPI_Isend(&d,1,MPI_INT,i,statag,MPI_COMM_WORLD,&request);
      }
    
    /* everybody recieves from everybody else */
    for (i = 1; i < nprocs; i++) {
      MPI_Recv(&d,1,MPI_INT,MPI_ANY_SOURCE,statag,MPI_COMM_WORLD,&status);
    }
#endif


  init_data(begin, end);
  MPI_Barrier(MPI_COMM_WORLD);
  gettimeofday(&start, NULL);

  if (begin & 1)
    sor_odd(begin, end);
  else 
    sor_even(begin, end);

  /*MPI_Barrier(MPI_COMM_WORLD);*/
  gettimeofday(&finish, NULL);

  printf("Elapsed time: %.2f seconds\n",
	 (((finish.tv_sec * 1000000.0) + finish.tv_usec) -
	  ((start.tv_sec * 1000000.0) + start.tv_usec)) / 1000000.0);

#ifdef STA
  printf("Data send: %d - Msg send: %d \n", _data_send, _msg_send);
#endif


  /* Some handshaking */
  /*
#ifdef STA
  statag = 11;
  if (proc_id != 0) {
    vsta[0] = _data_send;
    vsta[1] = _msg_send;
    MPI_Isend(&vsta,2,MPI_INT,0,statag,MPI_COMM_WORLD,&request);
  }
  else {
    for (i=0; i<nprocs-1; i++) {
      MPI_Recv(&vsta,2,MPI_INT,i+1,statag,MPI_COMM_WORLD,&status);
      _data_send += vsta[0];
      _msg_send += vsta[1];
    }
    printf("%d Messages, %d bytes of data sent. \n", _msg_send, _data_send);
  } 
#endif
*/

  MPI_Finalize();
  exit(0);
}


/***************************************************************************\
        init_data(begin, end)  - the mpi version of Initialize()
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


void init_data(begin, end)
int begin, end;
{
  extern double   exp();
  int i, j, k, length;

  /* The master should have the whole array */
  if (proc_id == 0)
    end = M;
  length = end - begin + 1;

  if ((red_ = (float **) malloc((length + 2)*sizeof(float *))) == 0) {
    fprintf(stderr, "out of shared memory\n");
    fflush(stderr);
    MPI_Finalize(); exit(-1);
  }
  
  if ((black_ = (float **) malloc((length + 2)*sizeof(float *))) == 0) {
    fprintf(stderr,"out of shared memory");
    fflush(stderr);
    MPI_Finalize(); exit(-1);
  }
  
  for (i=0; i <= length + 1; i++) {
    
    if ((red_[i] = (float *) malloc((N + 1)*sizeof(float))) == 0) {
      fprintf(stderr, "out of shared memory\n");
      fflush(stderr);
      MPI_Finalize(); exit(-1);
    }
    
    if ((black_[i] = (float *) malloc((N + 1)*sizeof(float))) == 0) {
      fprintf(stderr, "out of shared memory\n");
      fflush(stderr);
      MPI_Finalize(); exit(-1);
    }

    for (k = 0; k < N; k++) {
      black_[i][k] =(float) exp ((double) (k%35));
      red_[i][k] = (float) exp (k%35 + (double) 3);
    }
  }

}

    
    

/*
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
 *
 */
void	sor_odd(begin, end)
	int	begin;
	int	end;
{
  int	i, j, k, length;
  
  length = end - begin + 1;

  for (i = 0; i < iterations; i++) {
    
    for (j = 1; j <= length; j++) {
      
      for (k = 0; k < N; k++) {

	black_[j][k] = (red_[j-1][k] + red_[j+1][k] + red_[j][k] + red_[j][k+1])/4.0;
      }
      if ((j += 1) > length)
	break;
      
      for (k = 1; k <= N; k++) {
	
	black_[j][k] = (red_[j-1][k] + red_[j+1][k] + red_[j][k-1] + red_[j][k])/4.0;
      }
    }
    if (nprocs > 1)
      ex_black(length);

    for (j = 1; j <= length; j++) {
      
      for (k = 1; k <= N; k++) {

	red_[j][k] = (black_[j-1][k] + black_[j+1][k] + black_[j][k-1] + black_[j][k])/4.0;
      }
      if ((j += 1) > length)
	break;

      for (k = 0; k < N; k++) {

	red_[j][k] = (black_[j-1][k] + black_[j+1][k] + black_[j][k] + black_[j][k+1])/4.0;
      }
    }				

#ifdef      RESET_AFTER_ONE_ITERATION 
    if (i == 0) {
  	  MPI_Barrier(MPI_COMM_WORLD);
      gettimeofday(&start, NULL);
#ifdef STA
      _data_send = _msg_send = 0;
#endif
    } 
#endif

    if (nprocs > 1 && i < (iterations - 1))
      ex_red(length);
  }
}


/*
 * begin is even
 */
void	sor_even(begin, end)
     int	begin;
     int	end;
{
  int	i, j, k, length;

  length = end - begin + 1;
  for (i = 0; i < iterations; i++) {

    for (j = 1; j <= length; j++) {

      for (k = 1; k <= N; k++) {

	black_[j][k] = (red_[j-1][k] + red_[j+1][k] + red_[j][k-1] + red_[j][k])/4.0;
      }
      if ((j += 1) > length)
	break;

      for (k = 0; k < N; k++) {

	black_[j][k] = (red_[j-1][k] + red_[j+1][k] + red_[j][k] + red_[j][k+1])/4.0;
      }
    }
    if (nprocs > 1)
      ex_black(length);

    for (j = 1; j <= length; j++) {

      for (k = 0; k < N; k++) {

	red_[j][k] = (black_[j-1][k] + black_[j+1][k] + black_[j][k] + black_[j][k+1])/4.0;
      }
      if ((j += 1) > length)
	break;
      
      for (k = 1; k <= N; k++) {
	
	red_[j][k] = (black_[j-1][k] + black_[j+1][k] + black_[j][k-1] + black_[j][k])/4.0;
      }
    }		
#ifdef      RESET_AFTER_ONE_ITERATION 
    if (i ==0) {
      MPI_Barrier(MPI_COMM_WORLD);
      gettimeofday(&start, NULL);
#ifdef STA
      _data_send = _msg_send = 0;
#endif
    }
#endif
    if (nprocs > 1 && i < (iterations - 1))
      ex_red(length);

  }
}


void ex_black(length)
int length;
{
  int msgtag, r_pid, j;

  msgtag = 3;
  /* Send to (proc_id - 1) */
  if (proc_id != 0) {
       MPI_Isend(&black_[1][0],N+1,MPI_FLOAT,proc_id-1,msgtag,MPI_COMM_WORLD,&request);
#ifdef STA
    _msg_send++;
    _data_send += sizeof(int) + (N+1)*sizeof(float);
#endif
  }

  /* Send to (proc_id + 1) */
  if (proc_id != (nprocs-1)) {
       MPI_Isend(&black_[length][0],N+1,MPI_FLOAT,proc_id+1,msgtag,MPI_COMM_WORLD,&request);
#ifdef STA
    _msg_send++;
    _data_send += sizeof(int) + (N+1)*sizeof(float);
#endif
  }
  if (proc_id == 0) {
     MPI_Recv(&black_[length+1][0],N+1,MPI_FLOAT,proc_id+1,msgtag,MPI_COMM_WORLD,&status);
  }
  else {
     if (proc_id == (nprocs - 1)) {
        MPI_Recv(&black_[0][0],N+1,MPI_FLOAT,proc_id-1,msgtag,MPI_COMM_WORLD,&status);

     }
     else {
        flag1=flag2=0;
        MPI_Irecv(&black_[0][0],N+1,MPI_FLOAT,proc_id-1,msgtag,MPI_COMM_WORLD,&request1);
        MPI_Irecv(&black_[length+1][0],N+1,MPI_FLOAT,proc_id+1,msgtag,MPI_COMM_WORLD,&request2);
        while (flag1==0 &&  flag2==0) {
          MPI_Test(&request1, &flag1, &status1);
          MPI_Test(&request2, &flag2, &status2);
        } 
     }
  }
}


void ex_red(length)
int length;
{
  int msgtag, r_pid, j;

  msgtag = 4;
  /* Send to (proc_id - 1) */
  if (proc_id != 0) {
     MPI_Isend(&red_[1][0], N+1, MPI_FLOAT,proc_id-1,msgtag,MPI_COMM_WORLD,&request);
#ifdef STA
    _msg_send++;
    _data_send += sizeof(int) + (N+1)*sizeof(float);
#endif
  }

  /* Send to (proc_id + 1) */
  if (proc_id != (nprocs-1)) {
     MPI_Isend(&red_[length][0], N+1, MPI_FLOAT,proc_id+1,msgtag,MPI_COMM_WORLD,&request);
#ifdef STA
    _msg_send++;
    _data_send += sizeof(int) + (N+1)*sizeof(float);
#endif
  }

  if (proc_id == 0) {
     MPI_Recv(&red_[length+1][0],N+1,MPI_FLOAT,proc_id+1,msgtag,MPI_COMM_WORLD,&status);
  }
  else {
     if (proc_id == (nprocs - 1)) {
        MPI_Recv(&red_[0][0],N+1,MPI_FLOAT,proc_id-1,msgtag,MPI_COMM_WORLD,&status);
     }
     else {
        flag1=flag2=0;
        MPI_Irecv(&red_[0][0],N+1,MPI_FLOAT,proc_id-1,msgtag,MPI_COMM_WORLD,&request1);
        MPI_Irecv(&red_[length+1][0],N+1,MPI_FLOAT,proc_id+1,msgtag,MPI_COMM_WORLD,&request2);
        while (flag1==0 &&  flag2==0) {
          MPI_Test(&request1, &flag1, &status1);
          MPI_Test(&request2, &flag2, &status2);
        } 
     }
  }
}

