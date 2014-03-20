#define OUTPATH "output"  /*directory to write output*/

void init_common(char* config_file);

/* Divide the 2D slab amongst each process to maintain 
   appropriate load balancing */
void divide_work();

/* Initialize boundary conditions of the medium to be zero */
void init_zero_bounds();

/* Initialize the medium to have heat distributed normally */
void init_gaus_medium();

/* Prints data owned by individual process to file*/
void print_local_data();

/* Moves data for individual processes into continuous file 
   which holds all of the data. */
void print_global_data();

/* Frees malloced memory */
void clean_up();
