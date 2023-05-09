#!/bin/csh

#SBATCH --job-name=submit-blocksize-mpi.sh
#SBATCH -D .
#SBATCH --output=submit-blocksize-mpi.sh.o%j
#SBATCH --error=submit-blocksize-mpi.sh.e%j
## #SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

echo Starting computation at `date` in node `hostname`

#setenv HOST $(echo $HOSTNAME | cut -f 1 -d'.')

if ( "$HOSTNAME" == "boada-1") then
    echo "Use sbatch to execute this script"
    exit 0
endif

setenv PROG SW_mpi

# Make sure that all binaries exist
make $PROG

#set MPIRUN = mpirun
set MPIRUN = mpirun 

setenv size 36960
#setenv size 1234
setenv seq1 sequences/a_500k.dat
setenv seq2 sequences/b_500k.dat
setenv sm data.score
set gap = -1

# PROVIDE NUMBER OF PROCESSES FROM THE COMMAND LINE
setenv NP $1

if ($NP < 1) then
  echo A parameter is expected: The number of MPI processes
  exit 1
endif

set BS_list = "8 16 32 64 128 256 512 768 1024"
#set BS_list = "16 32 "

set i = 1

set basefn = $PROG-$NP-blocksize # Basename for file names generated
set out = $basefn.txt
set tmp = $basefn.tmp
set psfile = $basefn.ps

rm -rf $out
rm -rf $tmp
rm -rf ./elapsed.txt
foreach bs ( $BS_list )
	echo "BS=$bs" > $tmp
	#$MPIRUN -np $NP -machinefile $TMPDIR/machines ./$PROG $seq1 $seq2 $sm $gap $bs $size >> $tmp
	echo "$MPIRUN -np $NP ./$PROG $seq1 $seq2 $sm $gap $bs $size" >> $tmp 
	$MPIRUN -np $NP ./$PROG $seq1 $seq2 $sm $gap $bs $size >> $tmp 
        set result = `cat $tmp  | grep "Computation of scoring matrix time"| cut -d':' -f 2 | sort | tail -1`
	echo $i >> ./elapsed.txt
	echo $result >> ./elapsed.txt
	set i = `echo $i + 1 | bc -l`
	cat $tmp >> $out
end

set i = 1
rm -rf ./hash_labels.txt
foreach bs ( $BS_list )
	echo "hash_label at " $i " : " $bs >> ./hash_labels.txt
	set i = `echo $i + 1 | bc -l`
end

jgraph -P blocksize-mpi.jgr > $psfile
set usuario=`whoami`
set fecha=`date`
sed -i -e "s/UUU/$usuario/g" $psfile
sed -i -e "s/FFF/$fecha/g" $psfile
rm -rf ./hash_labels.txt
#rm -rf ./$tmp
echo Ending computation at `date` in node `hostname`
