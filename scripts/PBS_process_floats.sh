#!/bin/bash
#
#PBS -W group_list=yetiphyso
#PBS -N process_floats
#PBS -l nodes=1:ppn=1,walltime=12:00:00,mem=4000mb
#PBS -m abe

cd $PBS_O_WORKDIR
module add anaconda
floats_to_tables floats.epac.16deg.h5

