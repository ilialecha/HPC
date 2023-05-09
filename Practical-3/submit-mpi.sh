#!/bin/bash

#SBATCH --job-name=submit-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-mpi.sh.o%j
#SBATCH --error=submit-mpi.sh.e%j
## #SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

USAGE="\n USAGE: sbatch submit-mpi.sh Bsize N nproc\n
        Bsize  -> blocking factor\n
        N      -> length of the sequences to study\n
        nproc  -> number of processes to use\n"

if (test $# -lt 3 || test $# -gt 3)
then
     echo -e $USAGE
     exit 0
fi

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

if (test "${HOST}" = "boada-1")
then
    echo "Use sbatch to execute this script"
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


