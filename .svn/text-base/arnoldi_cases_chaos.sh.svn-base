#!/bin/sh -e 
#PBS -N LES_cases
#PBS -l nodes=1:ppn=12:nogpu
#PBS -l walltime=600:00:00
#PBS -l mem=24gb
#PBS -j oe

# Define la lista de nodos, y guarda el nombre del 
# archivo en la variable $MACHINEFILE

MACHINEFILE=$PBS_NODEFILE
#cat $PBS_NODEFILE
# Change to the working directory
#cd $WORKDIR
# set path and goto current directory
cd $PBS_O_WORKDIR
# check ... OK!
#echo $PBS_O_WORKDIR > workdir.out

# Get the number of processes to start
NP=`wc -l $PBS_NODEFILE | awk '{print $1}'`; export NP

# Which job launcher to use
MPIEXEC="mpiexec"
#
#EXE="/home/sekimoto/fortran/shear/trunk/shear/arnoldi"
EXE="/home/sekimoto/fortran/shear/trunk/modules/arnoldi"

#caselist="icase_list"
#caselist="icase_list_Ax3y1_5_re"
#caselist="icase_list_Ax3y2_re"
caselist="icase_list_Ax3y3_re"

input_arnoldi=input_arnoldi_$caselist

mTp=5
kgmres=100

while read arg1 arg2
do 
  echo computing icase=$arg1 
  [ ! -e icase$arg1 ] && mkdir -p icase$arg1
  [ ! -e icase$arg1/ini ] &&  mkdir -p icase$arg1/ini
  # use full path name, instead of ~/
  cp -p $arg2* icase$arg1/ini/.
  outdir=icase$arg1/arnoldi/
  [ ! -e $outdir ] && mkdir $outdir
  echo $arg2 > $input_arnoldi
  echo $outdir/arnoldi_sym5 >> $input_arnoldi
  #
  #cat   arnoldi.config >> input_arnoldi
  # initial perturbation vector 
  echo 1.d-3 >> $input_arnoldi
  # m box-period to run
  echo $mTp >> $input_arnoldi
  # dimensions of GMRES iterations
  echo $kgmres >> $input_arnoldi 
  # GMRES epsilon (should be 1.d-6)
  echo 1.d-6  >> $input_arnoldi
  $MPIEXEC -np $NP -hostfile $PBS_NODEFILE $EXE < $input_arnoldi > $outdir/out_arn_icase$arg1_$PBS_JOBID
done < $caselist
