/* Non-Blocking I/O multithreading strategy.
 */
#include <stdlib.h>
#include <stdio.h>
#include<assert.h>
#include<math.h>
#include<mpi.h>
#include "diffusion_common.h"

MPI_Request requests[4]; /* For ghost cell exchange of nb communication */

void init(char* config_file){
  init_common(config_file);
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

int main(int argc,char *argv[]){
  /* Initialize MPI and assert that a config file has been received
   * with parameters 
   */
  MPI_Init(&argc, &argv);
  assert(argc==2);
  init(argv[1]);
  
  nb_work();
  clean_up();
  
  MPI_Barrier(MPI_COMM_WORLD);
  print_global_data();
  
  MPI_Finalize();
  return 0;
}
