/* Blocking I/O multithreading strategy. */
#include <stdlib.h>
#include <stdio.h>
#include<assert.h>
#include<math.h>
#include<mpi.h>
#include "diffusion_common.h"

static long long int static_n; /* This will be the smaller of nx or ny */
static long long int cut_n; /* This will be the larger of nx or ny */
static long double ds;  /* Larger of 1/nx 1/ny */
static long double dc;  /* Smaller of 1/nx 1/ny */

static double kx;
static double ky;
static double ks;
static double kc;

static int nsteps;
static int wstep;

/* Interprocess communication Bookeeping */
static long long int g_end;    /* The end index of slab of diffusion medium */
static int psize;              /* Problem size along cut dimension for this proc */
static int residual;           /* residual term of problem division for load balancing*/

static int nprocs;        /* number of processes */
static int rank;          /* the rank of this process */
static int left;          /* rank of left neighbor */
static int right;         /* rank of right neighbor */
static double **medium0;  /* medium at t-1 */
static double **medium1;  /* medium at t */

void init(char* config_file){
  init_common(config_file);
  
}

void update(){
  int i,j;
  for(i=1; i<psize; i++){        /* Omit ghost points, 0 and psize+1 */
    for (j=1; j<static_n-1; j++){  /* Omit edges */
      medium1[i][j] = medium0[i][j] + 
	kc*(medium0[i+1][j] + medium0[i-1][j] - 2*medium0[i][j])+
	ks*(medium0[i][j+1] + medium0[i][j-1] - 2*medium0[i][j]);
    }
  }
  for(i=1; i<psize; i++){             /* Omit ghost points, 0 and psize+1 */
    for (j=1; j<static_n-1; j++){     /* Omit edges */
      medium0[i][j] = medium1[i][j];  /* Set medium[t-1]=medium[t] */
    }
  }
}

/* For now just implementing Blocking Strategy.  In case u forget:  */
/* int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                    int dest, int sendtag,
		    void *recvbuf, int recvcount, MPI_Datatype recvtype, 
		    int source, int recvtag,
		    MPI_Comm comm, MPI_Status *status) */
void exchange_ghost_cells(){
  MPI_Sendrecv(medium0[1], static_n, MPI_DOUBLE, left, 0,
	       medium0[psize+1], static_n, MPI_DOUBLE, right, 0,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Sendrecv(medium0[psize], static_n, MPI_DOUBLE, right, 0,
	       medium0[0], static_n, MPI_DOUBLE, left, 0,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
  
/* */
void do_work(){
  int t;
  for (t=1;t<=nsteps;t++){
    exchange_ghost_cells();
    update();
  }
}

int main(int argc,char *argv[]){
  /* Initialize MPI and assert that a config file has been received
   * with parameters 
   */
  MPI_Init(&argc, &argv);
  assert(argc==2);
  init(argv[1]);

  do_work();
  clean_up();

  MPI_Barrier(MPI_COMM_WORLD);
  print_global_data();
  MPI_Finalize();
  return 0;
}
