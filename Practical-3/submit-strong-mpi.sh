#!/bin/bash

#SBATCH --job-name=submit-strong-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-strong-mpi.sh.o%j
#SBATCH --error=submit-strong-mpi.sh.e%j
## #SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

USAGE="\n USAGE: sbatch submit-strong-mpi.sh\n
        \n
        This shell script performs a Strong Scalability analysis.\n
        It takes no parameters from the command line.\n"

if (test $# -gt 0)
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


SEQ=SW
PROG=SW_mpi

#size=32768
size=1500
size=3326400
size=369600
BS=256		# TO BE DONE: Change the value of the block size to the optimal one

np_NMIN=1
np_NMAX=12
N=2

# Make sure that all binaries exist
make $SEQ
make $PROG

#MPIRUN=mpirun
MPIRUN=mpirun

seq1=sequences/a_500k.dat
seq2=sequences/b_500k.dat
#seq1=sequences/10M.dat
#seq2=sequences/10M.dat
sm=data.score
gap=-1

HOST=$(echo $HOSTNAME | cut -f 1 -d'.')

basefn=$PROG-strong_${HOST}
psfile=$basefn.ps	# Postcript file with figures generated

out=$basefn.txt	        # File where we keep the execution results
tmp=$basefn.tmp	        # Temporary file where we save the execution results
timeout=$basefn-time.txt    # Auxiliary file used to keep output of time

outputpath1=./speedup1.txt
outputpath2=./speedup2.txt
outputpath3=./elapsed.txt
rm -rf $outputpath1 2> /dev/null
rm -rf $outputpath2 2> /dev/null
rm -rf $outputpath3 2> /dev/null

rm -rf $out 2> /dev/null

echo -n Executing $SEQ sequentially > $tmp
elapsed=0  # Acumulation of elapsed time of N executions of the program
partial=0  # Acumulation of score Matrix creation time of N executions

i=0        # Variable contador de repeticiones
while (test $i -lt $N)
do
  echo $'\n' >> $tmp
  /usr/bin/time --format=%e ./$SEQ $seq1 $seq2 $sm $gap $size >> $tmp 2>$timeout

  time1=`cat $timeout|tail -n 1`
  time2=`cat $tmp | grep "Computation of scoring matrix time"| cut -d':' -f 2 | sort | tail -1`

  elapsed=`echo $elapsed + $time1|bc`
  partial=`echo $partial + $time2|bc`

  cat $tmp >> $out
  echo "Elapsed time sequential=$time1" >> $out
  i=`expr $i + 1`
done

echo $'\n' >> $out
echo -n ELAPSED TIME AVERAGE OF $N EXECUTIONS = >> $out
sequential1=`echo $elapsed/$N|bc -l`
echo $sequential1 >> $out
echo -n COMPUTATION OF SCORING MATRIX TIME AVERAGE OF $N EXECUTIONS = >> $out
sequential2=`echo $partial/$N|bc -l`
echo $sequential2 >> $out

i=0

NP=$np_NMIN
while (test $NP -le $np_NMAX)
do
  echo $'\n' >> $out
  echo -n Executing $PROG in parallel with $NP threads >> $out
  elapsed=0  # Acumulation of elapsed time of N executions of the program
  partial=0  # Acumulation of score Matrix creation time of N executions

  while (test $i -lt $N)
  do
     echo $'\n' > $tmp
     #export OMP_NUM_THREADS=$NP
     #/usr/bin/time --format=%e ./$PROG -n $size $BS >> $tmp 2>$timeout

     /usr/bin/time --format=%e $MPIRUN -np $NP ./$PROG $seq1 $seq2 $sm $gap $BS $size >> $tmp 2>$timeout

     time1=`cat $timeout|tail -n 1`
     time2=`cat $tmp | grep "Computation of scoring matrix time"| cut -d':' -f 2 | sort | tail -1`

     elapsed=`echo $elapsed + $time1|bc`
     partial=`echo $partial + $time2|bc`

     rm -f $timeout
     cat $tmp >> $out
     echo "Elapsed time $NP processes=$time1" >> $out
     i=`expr $i + 1`
  done

  echo $'\n' >> $out
  echo -n ELAPSED TIME AVERAGE OF $N EXECUTIONS = >> $out
  average1=`echo $elapsed/$N|bc -l`
  result1=`echo $sequential1/$average1|bc -l`
  echo $average1 >> $out
  echo -n COMPUTATION OF SCORING MATRIX TIME AVERAGE OF $N EXECUTIONS = >> $out
  average2=`echo $partial/$N|bc -l`
  result2=`echo $sequential2/$average2|bc -l`
  echo $average2 >> $out

  i=0

  #output NP i average en fichero speedup1
  echo -n $NP >> $outputpath1
  echo -n "   " >> $outputpath1
  echo $result1 >> $outputpath1

  #output NP i average en fichero speedup2
  echo -n $NP >> $outputpath2
  echo -n "   " >> $outputpath2
  echo $result2 >> $outputpath2

  #output NP i average en fichero elapsed
  echo -n $NP >> $outputpath3
  echo -n "   " >> $outputpath3
  echo $average2 >> $outputpath3

  #incrementa el parametre
  NP=`expr $NP + 1`
done

jgraph -P strong-mpi.jgr >  $psfile
usuario=`whoami`
fecha=`date`
sed -i -e "s/UUU/$usuario/g" $psfile
sed -i -e "s/FFF/$fecha/g" $psfile

rm -rf ./$tmp
