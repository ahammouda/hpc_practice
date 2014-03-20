/*
   To run code:
   $ mpiexec -n <number-procs> ./a.out <dimension1-size> <dimension2-size> \
                  <parallelization-strategy>
   Ex:
   $ mpiexec -n 4 ./a.out 10000 10000 1
   
   <parallelization-strategy> can be:
   0  statically load balanced 
   1  dynamically load balanaced (master/slave)
   
   Ways to display data from the 2 strategies.
   --------------------------------------------------
   Matlab to display output data from statically load
   balanced strategy:
   --------------------------------------------------
   clear; clc;
   fp=fopen('mandelbrot.bin')
   A=fread(fp,'double');
   imagesc(reshape(A,100,100));
   
   --------------------------------------------------
   Matlab to display output from the dynamically load 
   balanced strategy.
   --------------------------------------------------
   clear; clc;
   fp=fopen('mandelbrot.bin')
   A=fread(fp,'double');
   fp2=fopen('order.out')
   O=fread(fp2,'int');
   B=zeros(100,100);
   for i=1:100
      ind=O(i,1)+1; % All this data is indexed from zero
      inds=(i-1)*100 + 1;
      inde=i*100;
      B(:,ind)=A(inds:inde,1);
   end
   imagesc(B);
*/
#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
#include "mpi.h"

static FILE* file;
static FILE* orderfp;

static int WORKTAG=1;
static int DIETAG=2;
static int MAX_ITERS=1000;

static int nx=-1;
static int ny=-1;
static int pstrategy=-1;

static int myrank=-1;
static int nworkers=-1;

typedef struct complex{
  double real;
  double imag;
} Complex;

int cal_pixel(Complex c){
  int count;
  Complex z;
  double temp, lengthsq;
  
  z.real = 0;
  z.imag = 0;
  count = 0;
  
  do{
    temp = z.real * z.real - z.imag * z.imag + c.real;
    z.imag = 2 * z.real * z.imag + c.imag;
    z.real = temp;
    lengthsq = z.real * z.real + z.imag * z.imag;
    count ++;
  }
  while ((lengthsq < 4.0) && (count < MAX_ITERS));
  return(count);
}

double* do_work(double* work, int curjob)
{
  double *data_l, *data_l_tmp;
  int tmp;
  Complex c;
  int i, j;
  /* create buffer for local work only */
  data_l = work;
  data_l_tmp = data_l;
  
  /* calc each procs coordinates and call local mandelbrot set function */
  i = curjob;
  for (j = 0; j < ny; ++j){
    c.real = i/((double) nx) * 4. - 2. ;
    c.imag = j/((double) ny) * 4. - 2. ;
    tmp = cal_pixel(c);                       /* Calculate Pixel */
    *data_l++ = (double) tmp;
  }
  data_l = data_l_tmp;
  return data_l;
}

void master()
{
  int rank;
  double *result;
  MPI_Status status;
  
  int curjob=0;
  int njobs=nx;  /* Work is doled out in rows */
  
  /* Seed the slaves; send one unit of work to each slave. */
  for (rank = 1; rank < nworkers; ++rank) {
    
    /* Send it to each rank */
    MPI_Send(&curjob,           /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
    
    /* Find the next item of work to do */
    curjob++;
  }
  printf("sent %d jobs to %d workers \n",curjob, nworkers);
  fflush(stdout);
  
  /* Loop over getting new work requests until there is no more work
     to be done */
  result = (double *) malloc( ny * sizeof(double) );
  orderfp=fopen("order.out","w");
  
  while (curjob < njobs) {
    if (curjob==(nworkers-2))
      file = fopen("mandelbrot.bin", "w");
    else
      file = fopen("mandelbrot.bin", "a");
    
    /* Receive results from a workers as they finish. */
    MPI_Recv(result,
	     ny,
	     MPI_DOUBLE, 
	     MPI_ANY_SOURCE, 
	     MPI_ANY_TAG, 
	     MPI_COMM_WORLD, 
	     &status);
    fwrite(&status.MPI_TAG,1,sizeof(int),orderfp);
    fwrite(result, ny, sizeof(double), file);
    fclose(file);
        
    printf("received result from worker number %d with %d jobs done \n",
	   status.MPI_SOURCE, curjob);
    fflush(stdout);
    
    /* Get the next unit of work to be done */
    curjob++;
    
    /* Send the slave a new work unit */
    MPI_Send(&curjob,           /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
    printf("sent job number %d to worker %d\n",curjob, status.MPI_SOURCE);
    fflush(stdout);
  }

  /* There's no more work to be done, so receive all the outstanding
     results from the slaves. */
  for (rank = 1; rank < nworkers; ++rank) {
    MPI_Recv(result, ny, MPI_DOUBLE, rank,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    printf("collecting last job from %d \n",rank);
    fflush(stdout);
    fwrite(&status.MPI_TAG,1,sizeof(int),orderfp);
    file = fopen("mandelbrot.bin", "a");
    fwrite(result, ny, sizeof(double), file);
    fclose(file);
  }
  fclose(orderfp);
  free(result);
  
  /* Tell all the slaves to exit by sending an empty message with the
     DIETAG. */
  for (rank = 1; rank < nworkers; ++rank) {
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    printf("Killing node %d ...\n", rank);
    fflush(stdout);
  }
}


void slave()
{
  int curjob;
  double* result;
  MPI_Status status;
  
  while(1){
    /* Receive a message from the master (the job count) */
    MPI_Recv(&curjob, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    printf("A slave received job number %d\n", curjob);
    fflush(stdout);
    
    /* Check the tag of the received message. */
    if (status.MPI_TAG == DIETAG) {
      printf("Killing one slave node\n");
      fflush(stdout);
      break;
    }
    /* Allocate an array of doubles for the job */
    result = (double *) malloc( 100 * sizeof(double) );
    
    /* Do the work */
    result = do_work(result, curjob);
    printf("%d finished it's work\n", curjob);
    fflush(stdout);
    
    printf("Sending over elements in our result\n");
    fflush(stdout);
    /* Send the result back */
    MPI_Send(result, ny, MPI_DOUBLE, 0, curjob, MPI_COMM_WORLD);
    
    free(result);
    printf("%d freed its memory \n", curjob);
    fflush(stdout);
  }
}

void master_slave(){
  if (myrank != 0) {
    printf("node %d calling slave() \n", myrank);
    fflush(stdout);
    slave();
    printf("node %d finished slave() \n", myrank);
    fflush(stdout);
  } else {
    master();
    printf("node %d finished master() \n", myrank);
    fflush(stdout);
  }
}

void static_loads(){
  
  int i, j;
  Complex c;
  int tmp;
  double *data_l, *data_l_tmp;
  int mystrt, myend;
  int nrows_l;
  
  MPI_Status status;
  /********************************************************************/
  /********************   DIVIDE WORK    ******************************/
  
  /* assume divides equally */
  nrows_l = nx/nworkers;
  
  /* create buffer for local work only */
  data_l = (double *) malloc(nrows_l * ny * sizeof(double));
  data_l_tmp = data_l;
  
  /* calculate each processor's region of work 
   * REMEMBER -- We've already declared the number of processors
   * we're going to use, so we guarantee here that each of these
   * get's work by limiting scope of batch with process id (myrank)
   */
  mystrt = myrank*nrows_l;
  myend  = mystrt + nrows_l - 1;
  
  /* calc each procs coordinates and call local mandelbrot set function */
  for (i = mystrt; i <= myend; ++i){
    c.real = i/((double) nx) * 4. - 2. ;
    for (j = 0; j < ny; ++j){
      c.imag = j/((double) ny) * 4. - 2. ;
      tmp = cal_pixel(c);
      *data_l++ = (double) tmp;
    }
  }
  data_l = data_l_tmp;
  
  /********************   END DIVIDE WORK    ***************************/
  /********************************************************************/
  
  printf("This is thread: %d \nStart: %d \nEnd: %d\n ",myrank, mystrt, myend);
  
  if (myrank == 0){
    /* Initial Opening and write masterpe writes his own calculations*/
    file = fopen("mandelbrot.bin", "w");
    printf("nrows_l, ny  %d %d\n", nrows_l, ny);
    fwrite(data_l, nrows_l*ny, sizeof(double), file);
    fclose(file);
    
    /* Then append the calculations of every other thread */
    for (i = 1; i < nworkers; ++i){
      /* 
	 args (in order): 
	 buffer, number of elements in buff, type in each elemenet, source, tag, 
	 comm, status_code.
      */
      MPI_Recv(data_l, nrows_l * ny, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status); 
      printf("received message from proc %d\n", i);
      file = fopen("mandelbrot.bin", "a");
      fwrite(data_l, nrows_l*ny, sizeof(double), file);
      fclose(file);
    }
  }
  
  else{
    /* 
       args (in order):
       buffer, number of elements in buff, type in each elemenet, destination, 
       tag, comm 
    */
    MPI_Send(data_l, nrows_l * ny, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }  
}

int main(int argc, char **argv)
{
  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  /* Find out my identity in the default communicator */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  /* Find out how many processes there are in the default communicator */
  MPI_Comm_size(MPI_COMM_WORLD, &nworkers);
  
  if (argc != 4){
    int err = 0;
    printf("argc %d\n", argc);
    if (myrank == 0){
      printf("usage: mandelbrot nx ny");
      MPI_Abort(MPI_COMM_WORLD,err );
    }
  }
  
  printf("Number of processes initialized:  %d \n", nworkers);
  
  if (myrank==0){
    nx=atoi(argv[1]);
    ny=atoi(argv[2]);
    pstrategy=atoi(argv[3]); /* 0 for static 1 for dynamic loads */
  }
  MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pstrategy, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  if(pstrategy){
    master_slave();
  } else {
    static_loads();
  }
  /* Shut down MPI */
  MPI_Finalize();
  return 0;
}
