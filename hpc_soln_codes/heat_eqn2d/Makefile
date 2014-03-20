all:
	mpicc -c diffusion_util.c diffusion_mpi.c
	mpicc -o heat2d diffusion_util.o diffusion_mpi.o -lm -std=c99
clean-all:
	rm -f heat2d *~ *.o output/*.dat output/*~ output/*.zip
clean-out:
	rm -f output/*.dat output/*~ output/*.zip
x:
	mpiexec -n 11 heat2d REFS/config