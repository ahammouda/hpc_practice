#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

/* Allocate 2D matrix.  Nomenclature coming from nrutil, but
 * I'm in control ! >:) */
/* ASSUMPTIONS:: */
/* More of a list of vectors--Makes ghost exchanges a bit more
 * intuitive. */
double **dmatrix(int e_static,int e_cut){
  double**  m; int i;
  
  m = (double**)malloc( e_cut*sizeof(double) );
  if(!m) fprintf(stderr,"Allocation Failure, 1st dimension");
  
  for(i=0;i<e_cut;i++) {
    m[i] = (double*)malloc( e_static*sizeof(double) );
    if (!m[i]) fprintf(stderr,"Allocation Failure, %dth row",i);
  }
  return m;
}

/* TODO: This may ignore some boundary conditions (Think about them) */
void free_dmatrix(double** dmatrix,int cut_dim){
  int i;
  for(i=0;i<cut_dim;i++)
    free(dmatrix[i]);
  free(dmatrix);
}

/* Prints usage information and exits with return status 1. */
static void quit() {
  printf("Input file must have format:\n\n");
  printf("nx = <INTEGER>\n");
  printf("k = <DOUBLE>\n");
  printf("nsteps = <INTEGER>\n");
  printf("wstep = <INTEGER>\n");
  printf("delay_dur = <DOUBLE>\n");
  printf("delay_per_mean = <DOUBLE>\n");
  printf("delay_per_stddev = <DOUBLE>\n");
  fflush(stdout);
  exit(1);
}

/* Reads a key-value pair from a file, where the value is an integer.
 * Takes as input a file pointer, a key (string), and a pointer to an
 * integer.  The text pointed to in the file must have the form "key =
 * VALUE".  This function checks that the key read matches the given
 * key; if they fail to match, an error is reported and execution
 * halts.  It then parses the integer and stores the resulting integer
 * value in the block pointed to by the given integer pointer.  Errors
 * are reported if anything does not parse.
 */
static void readint(FILE *file, char *keyword, int *ptr) {
  char buf[101];
  int value;
  int returnval;

  returnval = fscanf(file, "%100s", buf);
  if (returnval != 1) quit();
  if (strcmp(keyword, buf) != 0) quit();
  returnval = fscanf(file, "%10s", buf);
  if (returnval != 1) quit();
  if (strcmp("=", buf) != 0) quit();
  returnval = fscanf(file, "%d", ptr);
  if (returnval != 1) quit();
}

static void readlong(FILE *file, char *keyword, long long int *ptr) {
  char buf[101];
  long value;
  int returnval;

  returnval = fscanf(file, "%100s", buf);
  if (returnval != 1) quit();
  if (strcmp(keyword, buf) != 0) quit();
  returnval = fscanf(file, "%10s", buf);
  if (returnval != 1) quit();
  if (strcmp("=", buf) != 0) quit();
  returnval = fscanf(file, "%lld", ptr);
  if (returnval != 1) quit();
}


/* Read a key-value pair from a file, where the value is a double                                                                                                                 * precision floating point number.   See comments for readint.                                                                                                                   */
static void readdouble(FILE *file, char *keyword, double *ptr) {
  char buf[101];
  int value;
  int returnval;

  returnval = fscanf(file, "%100s", buf);
  if (returnval != 1) quit();
  if (strcmp(keyword, buf) != 0) quit();
  returnval = fscanf(file, "%10s", buf);
  if (returnval != 1) quit();
  if (strcmp("=", buf) != 0) quit();
  returnval = fscanf(file, "%lf", ptr);
  if (returnval != 1) quit();
}

/*
 * Note that k=(diffusivity_Const)*dt/dx^2
 */
void  readConfig(char *infilename,
		 long long *nx, long long *ny,
		 double *kx, double *ky, int *nsteps, int *wstep) {
  FILE *infile = fopen(infilename, "r");
  assert(infile);
  readlong(infile, "nx", nx);
  readlong(infile, "ny", ny);
  readdouble(infile, "kx", kx);
  readdouble(infile, "ky", ky);
  readint(infile, "nsteps", nsteps);
  readint(infile, "wstep", wstep);
  printf("diffusion1d: nx=%lld, ny=%lld, kx=%.2le, ky=%.2le, nsteps=%d, wstep=%d \n",
         *nx, *ny, *kx, *ky, *nsteps, *wstep);
  fflush(stdout);
  assert(*kx>0 && *kx<.5);
  assert(*ky>0 && *ky<.5);
  assert(*nx>=2);
  assert(*nsteps>=1);
  fclose(infile);
}
