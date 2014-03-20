/*
 * Need a realistic implementation of the 2D heat-equation for MPI.
 */
#include <stdlib.h>
#include <stdio.h>
#include<assert.h>
#include<math.h>
#include<mpi.h>
#include "diffusion_util.h"

#define OUTPATH "output"  /*directory to write output*/

static int BLOCKING = 0;

/* Global Problem Specs */
static long long int nx;
static long long int ny;
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

MPI_Request requests[4]; /* For ghost cell exchange of nb communication */

void print_local_data();

void divide_work(){
  int residual,xBig;
  /* Make sure we divide the longer of the 2 dims amongst procs */
  xBig = (nx>ny) ? 1 : 0;
  static_n = (xBig) ? ny : nx; /* Determine which dimension to keep */
  cut_n = (xBig) ? nx : ny;    /* Determine which dimension to cut */
  ks = (xBig) ? ky : kx;
  kc = (xBig) ? kx : ky;
  
  ds=1.0/(long double)static_n;
  dc=1.0/(long double)cut_n;
  
  residual = cut_n%nprocs;
  fprintf(stderr,"residual=%d \n",residual);
  
  psize = (rank<residual) ? floor(cut_n/nprocs)+1 : floor(cut_n/nprocs);
  g_end = (rank<residual) ? psize*(rank+1) : ((psize+1)*residual + psize*((rank+1)-residual));
  
  /* Assign neighbors to each proc */
  left = (rank==0) ? MPI_PROC_NULL : rank-1;
  right = (rank==nprocs-1) ? MPI_PROC_NULL : rank+1;
  
  /* Allocate medium0 and medium1 */
  medium0 = dmatrix(static_n,psize+2);
  medium1 = dmatrix(static_n,psize+2);
}

/* Zeroed boundaries w.r.t the global border of 
 * stencil */
/* TODO: Develop function for
 *  - Periodic Boundary conditions   
 *  - An Ar/Ind Source Cond S(i,j,k) */
void init_zero_bounds(){
  long long int i;
  if (rank == 0) {
    for (i=0;i<static_n;i++) medium0[1][i]=0; /* Remember that even bounds have G.C.*/
  } else if ( rank == nprocs-1) {
    for (i=0;i<static_n;i++) medium0[psize][i]=0; /* Ditto For other end */
  }
  for (i=1;i<psize+1;i++) medium0[i][0]=medium0[i][static_n]=0;
}

/* Does book-keeping to ensure medium is smooth gausian
 * w.r.t the global decomposition */
void init_gaus_medium(){
  long long int i,j,jj;
  fprintf(stderr,"rank=%d, g_end=%lld, psize=%d, nprocs=%d \n",rank,g_end,psize,nprocs);
  
  for (i=1;i<static_n-1;i++){
    jj=g_end-psize;
    for(j=1;j<psize+1;j++){
      medium0[j][i] = (exp(-( ((5*ds*i) - 2.5)*((5*ds*i) - 2.5) ) ) )*
	(exp(-( ((5*dc*jj) - 2.5)*((5*dc*jj) - 2.5) ) ) );
      jj++;
    }
  }
}

void init(char* config_file){
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* rank = [0,n-1] */
  
  if (rank==0){
    readOldConfig(config_file,&nx,&ny,&kx,&ky,&nsteps,&wstep);
  }
  MPI_Bcast(&nx, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ny, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ky, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&wstep, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  divide_work();
  init_gaus_medium();
  init_zero_bounds();
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

/*
  int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)
  int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request)
*/
void init_ghost_cell_exchange(){
  MPI_Irecv(medium0[psize+1], static_n, MPI_DOUBLE, right, 0,
	    MPI_COMM_WORLD, &requests[0] );
  MPI_Irecv(medium0[0], static_n, MPI_DOUBLE, left, 0,
	    MPI_COMM_WORLD, &requests[1] );
  MPI_Isend(medium0[1], static_n, MPI_DOUBLE, left, 0,
	    MPI_COMM_WORLD, &requests[2] );
  MPI_Isend(medium0[psize], static_n, MPI_DOUBLE, right, 0,
	    MPI_COMM_WORLD, &requests[3] );
}

/*
  int MPI_Waitall(int count, MPI_Request array_of_requests[], 
                  MPI_Status array_of_statuses[])
*/
void fin_ghost_cell_exchange(){
  MPI_Waitall(4,requests,MPI_STATUS_IGNORE);
}

/* Omit ghost points, and points that depend on them: [0-1] & [psize-psize+1] */
/* Omit edges of grid (Consider them heat sinks) */
void update_interior(){
  int i,j;
  for(i=2; i<psize-1; i++){
    for (j=1; j<static_n-1; j++){
      medium1[i][j] = medium0[i][j] + 
	kc*(medium0[i+1][j] + medium0[i-1][j] - 2*medium0[i][j])+
	ks*(medium0[i][j+1] + medium0[i][j-1] - 2*medium0[i][j]);
    }
  }
  for(i=2; i<psize-1; i++){
    for (j=1; j<static_n-1; j++){
      /*Set medium[t-1]=medium[t]*/
      medium0[i][j] = medium1[i][j];
    }
  }
}

void update_exterior(){
  int j;
  for (j=1; j<static_n-1; j++){
    medium1[1][j] = medium0[1][j] + 
      kc*(medium0[2][j] + medium0[0][j] - 2*medium0[1][j])+
      ks*(medium0[1][j+1] + medium0[1][j-1] - 2*medium0[1][j]);
    
    medium1[psize][j] = medium0[psize][j] + 
      kc*(medium0[psize+1][j] + medium0[psize-1][j] - 2*medium0[psize][j])+
      ks*(medium0[psize][j+1] + medium0[psize][j-1] - 2*medium0[psize][j]);
  }
}

void nb_work(){
  int t;
  for (t=1;t<=nsteps;t++){
    init_ghost_cell_exchange();
    update_interior();
    fin_ghost_cell_exchange();
    update_exterior();
  }
}
  
/*
  TODO:  
  * Place option for blocking or non-blocking as a parameter here.
 */
void do_work(){
  int t;
  for (t=1;t<=nsteps;t++){
    exchange_ghost_cells();
    update();
  }
}

/* Should print rows of the globally shorter dimension */
/* TODO:
   Be sure to be able to pass the option of  which step 
   is being printed, and make a directory (somewhere else)
   for that step. */
void print_local_data(){
  FILE* fptr;
  char filename[15];
  int i,j;
  sprintf(filename,"%s/%d.out2d.dat",OUTPATH,rank);
  fptr = fopen(filename,"a");
  for(i=0;i<psize+2;i++){
    for(j=0;j<static_n;j++){
      fprintf(fptr,"%3.2f,",medium0[i][j]);
    }
    fprintf(fptr,"\n");
  }
  fclose(fptr);
}

/*  TODO::
    Be sure to add a parameter and appropriate logic 
    to indicate which timestep data is being aggregated for
 */
void print_global_data(){
  FILE* rfptr; FILE* afptr;
  char ch;   int rank;
  char filename[15];
  sprintf(filename,"%s/global_out2d.dat",OUTPATH);
  afptr = fopen(filename,"a");
  
  for(rank=0;rank<nprocs;rank++){
    char filename[15];
    sprintf(filename,"%s/%d.out2d.dat",OUTPATH,rank);
    rfptr=fopen(filename,"r");
    /* Read From Local, Write to global */
    while( (ch=getc(rfptr)) != EOF )
      putc(ch,afptr);
    fclose(rfptr);
  }
  fclose(afptr);
}

void clean_up(){
  print_local_data();
  free_dmatrix(medium0,psize);
  free_dmatrix(medium1,psize);
}

int main(int argc,char *argv[]){
  /* Initialize MPI and assert that a config file has been received
   * with parameters 
   */
  MPI_Init(&argc, &argv);
  assert(argc==2);
  init(argv[1]);
  
  if (BLOCKING)  /*Specify globale bool for*/
    do_work();
  else
    nb_work();
  
  clean_up();
  MPI_Barrier(MPI_COMM_WORLD);
  print_global_data();
  MPI_Finalize();
  return 0;
}
