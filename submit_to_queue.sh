#!/bin/csh -fx
#PBS -N Iceberg_Weddell_Test3
#PBS -l walltime=0:40:00
#PBS -l size=32
#PBS -S /bin/tcsh
#PBS -r n
#PBS -m ae
#PBS -j oe
#PBS -E
#PBS -A gfdl_o

#set compiler = intel


module load totalview
env
pwd
time aprun -n 32 ../../build/intel/ice_ocean_SIS2/repro/MOM6
#date

#time aprun ...

#date
