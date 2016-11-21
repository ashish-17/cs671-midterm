#!/bin/bash

#SBATCH -J AJ523_MIDTERM
#SBATCH -o AJ523_MIDTERM.%J.stdout
#SBATCH -e AJ523_MIDTERM.%J.stderr
#SBATCH -p main
#SBATCH -N 2
#SBATCH -t 00:10:00

#Uncomment following lines before running on caliburn
#cd $HOME/cs671-midterm
#module load openmpi
#sleep 3

set -e

make clean

find . -type f -name '*.csv' -delete

make


size=1000
max=10000
while [ "$size" -lt "$max" ] 
do
    ./main 0 20 $size >> "stats_serial.csv"
    ./main 1 20 $size 2 >> "stats_omp_2.csv"
    ./main 1 20 $size 4 >> "stats_omp_4.csv"
    ./main 1 20 $size 8 >> "stats_omp_8.csv"
    ./main 1 20 $size 16 >> "stats_omp_16.csv"
    ./main 1 20 $size 24 >> "stats_omp_24.csv"
    ./main 1 20 $size 32 >> "stats_omp_32.csv"
    ./main 1 20 $size 64 >> "stats_omp_64.csv"
    size=$(($size+1000))
done
:
