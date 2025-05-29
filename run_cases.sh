#!/bin/sh
#
EXE="/home/sekimoto/fortran/shear/trunk/shear/disk3d_ini1"
#mpirun -np 4 $EXE 
caselist="icase_list"
ofile='var_sym5'

while read arg1 arg2
do 
  echo computing icase=$arg1 
  mkdir icase$arg1
  #mkdir icase$arg1/ini
  # use full path name, instead of ~/
  #cp -p $arg2* icase$arg1/ini/.
  outdir=icase$arg1/pool_run/
  mkdir $outdir
  echo $arg2 > input_run
  echo $outdir/$ofile >> input_run
  #
  mpirun -np 4 $EXE < input_run > $outdir/output_pool_run_icase$arg1
done < $caselist
