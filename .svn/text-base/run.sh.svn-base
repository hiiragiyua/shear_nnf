#!/bin/sh -e 
#PBS -N shear_test
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2:00:00
#PBS -l mem=4gb
#PBS -j oe

# switch ppn=12:nogpu or ppn12:gpus=1
# Define la lista de nodos, y guarda el nombre del 
# archivo en la variable $MACHINEFILE

MACHINEFILE=$PBS_NODEFILE

#cat $PBS_NODEFILE
# Change to the working directory
#WORKDIR=/wamba11/sekimoto/shear/vulcano/test1
#cd $WORKDIR
# set path and goto current directory
cd $PBS_O_WORKDIR
# check ... OK!
#echo $PBS_O_WORKDIR > workdir.out

# Get the number of processes to start
NP=`wc -l $PBS_NODEFILE | awk '{print $1}'`; export NP

# Which job launcher to use
MPIEXEC="mpirun"

# Executable binary
EXE="../disk3d"
#cat $PBS_NODEFILE > nodes_$PBS_JOBID.out
#echo $NP > nproc_$NP
$MPIEXEC -np $NP -hostfile $PBS_NODEFILE $EXE > output_$PBS_JOBID
