#!/bin/bash
#PBS -N FLP_kazeroon_n1
#PBS -o /public/home/leisichao/FLP_n1/FLP_kazeroon_n1.out
#PBS -e /public/home/leisichao/FLP_n1/FLP_kazeroon_n1.err
#PBS -q queue2
#PBS -l nodes=1:ppn=1
#PBS -M 1120878054@qq.com
#PBS -m abe
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
/opt/MATLAB/R2018b/bin/matlab -nodisplay -nosplash -nojvm -r "addpath('/public/home/leisichao/FLP_n1/');main_runtimes_kazeroon_n1;exit;"

