#!/bin/sh
#
EXE="/home/sekimoto/fortran/shear/trunk/shear/arnoldi"

caselist="icase_list"
mTp=5
kgmres=100

while read arg1 arg2
do 
  echo computing icase=$arg1 
  mkdir icase$arg1
  mkdir icase$arg1/ini
  # use full path name, instead of ~/
  cp -p $arg2* icase$arg1/ini/.
  outdir=icase$arg1/arnoldi/
  mkdir $outdir
  echo $arg2 > input_arnoldi
  echo $outdir/arnoldi_sym5 >> input_arnoldi
  #
  #cat   arnoldi.config >> input_arnoldi
  # initial perturbation vector 
  echo 1.d-3 >> input_arnoldi
  # m box-period to run
  echo $mTp >> input_arnoldi
  # dimensions of GMRES iterations
  echo $kgmres >> input_arnoldi 
  # GMRES epsilon (should be 1.d-6)
  echo 1.d-6  >> input_arnoldi
  mpirun -np 4 $EXE < input_arnoldi > $outdir/output_arnoldi_icase$arg1
done < $caselist
