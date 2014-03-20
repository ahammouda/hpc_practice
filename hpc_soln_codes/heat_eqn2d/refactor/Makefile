par:
	mpicc -c diffusion_util.c diffusion2d_common.c diffusion2d_par.c
	mpicc -o par2dx diffusion_util.o diffusion2d_common.o diffusion2d_par.o -lm -std=c99
nb:
	mpicc -c diffusion_util.c diffusion2d_common.c diffusion2d_nb.c
	mpicc -o nb2dx diffusion_util.o diffusion2d_common.o diffusion2d_nb.o -lm -std=c99
clean-all:
	rm -f par2dx nb2dx *~ *.o output/*.dat output/*~ output/*.zip
clean-out:
	rm -f output/*.dat output/*~ output/*.zip
x:
	mpiexec -n 11 heat2d REFS/config
