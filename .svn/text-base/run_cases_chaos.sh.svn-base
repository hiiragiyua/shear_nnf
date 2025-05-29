#!/bin/sh -e 
#PBS -N LES_cases_run
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
#!/bin/sh
#
EXE="/home/sekimoto/fortran/shear/trunk/modules/disk3d_ini1"
# 
#caselist="icase_list"
#caselist="icase_list_Ax3y1_5"
#caselist="icase_list_Ax3y2"
caselist="icase_list_Ax3y3"
ofile='var_sym5'
tmpfile='input_run_'$caselist
while read arg1 arg2
do 
  echo computing icase=$arg1 
  [ ! -e icase$arg1 ] && mkdir -p icase$arg1
  #mkdir icase$arg1/ini
  # use full path name, instead of ~/
  #cp -p $arg2* icase$arg1/ini/.
  outdir=icase$arg1/pool_run/
  [ ! -e $outdir ] && mkdir -p $outdir
  echo $arg2 > $tmpfile
  echo $outdir/$ofile >> $tmpfile
  #
  $MPIEXEC -np $NP -hostfile $PBS_NODEFILE $EXE < $tmpfile > $outdir/output_pool_run_icase$arg1
done < $caselist
