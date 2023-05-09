CC       = gcc
MPICC    = mpicc

#OPENMP   = -fopenmp
#CFLAGS   = -Wall -Wimplicit-function-declaration -O3 -std=c99 -march=native -Wno-maybe-uninitialized
CFLAGS   = -Wall -Wimplicit-function-declaration -O3 -std=c99 -march=native 

TARGETS = SW SW_mpi

all: $(TARGETS)

SW_mpi: SW_mpi.c
	$(MPICC) $< $(CFLAGS) $(CINCL) -o $@ $(CLIBS)

SW: SW.c
	$(CC) $< $(CFLAGS) -o $@

clean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e*
ultraclean:
	rm -rf $(TARGETS) *.o *.pe* *.po* *.o* *.e* *.txt *.ps
