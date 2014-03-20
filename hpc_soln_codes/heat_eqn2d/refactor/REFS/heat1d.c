#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>

void init(double *T, double *S, int n, double L, double dx);
void output(double *T,int n);

int main(int argc, char **argv){
  /* 1d array of temperatures */
  /* size 0 ... n+1 with bc's at 0 and n+1 */
  /* need both space for old and updated solution */
  double *T, *Tnew;   
  /* S is the heating rate; assume spatially dependent only */
  double *S;
  /* domain runs from 0 .. L */
  double L;
  /* diffusivity */
  double alpha;
  /*time step */
  double dt;
  /*grid spacing */
  double dx;
  /* number of grid points not including boundaries */
  int n;
  /* maximum number of iterations */
  int iter_max = 10000;
  /* print solution every iprint timesteps */
  int iprint = 5000;
  /* stuff */
  int iter, i;


  /* read input */
  n     = atoi(argv[1]);
  L     = atof(argv[2]);
  alpha = atof(argv[3]);
  dt    = atof(argv[4]);

  T    = (double *) malloc( (n+2)*sizeof(double));
  Tnew = (double *) malloc( (n+2)*sizeof(double));
  S    = (double *) malloc( (n+2)*sizeof(double));
  assert(T); assert(Tnew); assert(S);
  
  /* grid spacing */
  dx    = L/(n+1);
  init(T,S,n,L,dx);
  output(T,n);
  
  
  for (iter=1; iter<=iter_max; ++iter){
    if (!(iter % iprint)){
      output(T,n);
      exit(0);
    }
    printf("iteration: %d\n", iter);
    for (i=1;i<=n;++i){
      Tnew[i] = T[i] + alpha*dt/(dx*dx) *(T[i+1] - 2*T[i] + T[i-1]);
    }
    memcpy(T,Tnew,(n+2)*sizeof(double));
  }

}

/*
 * create initial condition
 */
void init(double *T, double *S, int n, double L, double dx){
  int i;
  double x = -5.;
  
  for (i=1; i<=n; ++i){
    x += dx*10./L;
    T[i] = exp(-x*x/2.);
    S[i] = 0.;
 }
  /* set bc's explicitly */
  T[0] = 0.;
  T[n+1] = 0.;
    
}

void output(double *T,int n){
  int i;
  FILE *outfile;
  outfile = fopen("heat1d.txt","w");
  assert(outfile);
  for (i = 0; i <= n+1; ++i){
    fprintf(outfile, "%f\n", T[i]);
  }
  fclose(outfile);
}
