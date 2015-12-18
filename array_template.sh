#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -j oe
#PBS -J 0-{scriptcount}
module load python-uon
cd {path}
echo "Running on `hostname`"
echo  "The TASK ID of *this* task is : $PBS_ARRAY_INDEX"
./{scriptcommand} --index=$PBS_ARRAY_INDEX
echo "Finished job now"
