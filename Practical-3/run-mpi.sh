#!/bin/bash


USAGE="\n USAGE: run-mpi.sh Bsize N nproc\n
        Bsize  -> blocking factor\n
        N      -> length of the sequences to study (< 6000)\n
        nproc  -> number of processes to use (<4)\n"

if (test $# -lt 3 || test $# -gt 3)
then
     echo -e $USAGE
     exit 0
fi

if (test $2 -gt 6000 || test $3 -gt 4)
then
	echo -e $USAGE
        exit 0
fi

PROGRAM=SW_mpi

seq1=sequences/a_500k.dat
seq2=sequences/b_500k.dat
sm=data.score
gap=-1 

make $PROGRAM


# Change the gap penalty if you wish to study the influence of this parameter

# MPIRUN=mpirun
MPIRUN=mpirun

$MPIRUN -np $3 $PROGRAM $seq1 $seq2 $sm $gap $1 $2


