
double **dmatrix(int e_static,int e_cut);

void free_dmatrix(double** dmatrix,int e_cut);

void  readConfig(char *infilename,
		 long long *nx, long long *ny, 
		 double *kx, double *ky, int *nsteps, int *wstep);
