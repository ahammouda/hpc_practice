float **fmatrix(int e_static,int e_cut);

void free_fmatrix(float** fmatrix,int cut_dim);

double **dmatrix(int e_static,int e_cut);

void free_dmatrix(double** dmatrix,int e_cut);

void readConfig(char *infilename,int* bs,
		long long *nx, long long *ny,
		double *k, int *nsteps, int *wstep,
		double *delay_dur, double *delay_per_mean,
		double *delay_per_stddev);

void  readOldConfig(char *infilename,
		    long long *nx, long long *ny,
		    double *kx, double *ky, int *nsteps, int *wstep);
