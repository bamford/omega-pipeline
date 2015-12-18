#!/usr/bin/env python
import numpy as np
from os.path import join
from glob import glob

project_path = '/home/ppzsb1/projects/OMEGA'

script_path = join(project_path, 'mcmc_runs')
script_template = 'run_mcmc_evidence.py --aper={aper} --field={field}'

data_path = join(project_path, 'plot_aper/all_spectra_hst')

with open(os.path.join(script_path, 'array_template.sh')) as f:
    template = f.readlines()

#make array_mcmc_pt_aper5_field*.sh
for field in fieldlist:
    info = dict(path=script_path, scriptcommand=scriptcommand, scriptcount=scriptcount)
    outlines = [line.format(info) for line in template]
    outfilename = 'array_mcmc_pt_{}_field{}.sh'.format(aper=aper, field=field)
    with open(outfilename, 'w') as outfile:
        outfile.writelines(outlines)

#make array_mcmc_pt_total_aper_field*.sh
for key in fieldDict.keys():
    fout = open('array_mcmc_pt_total_aper_field%d.sh'%(key),'w')
    fout.write('#!/bin/bash\n')

    fout.write('# This an example task array script\n')

    fout.write('# OPTIONS FOR PBS PRO ==============================================================\n')
    fout.write('#\n')
    fout.write('#PBS -l walltime=01:00:00\n')
    fout.write('# This specifies the job will run within 1 hour \'real\' time\n')

    fout.write('#PBS -l select=1:ncpus=1:mem=1gb\n')
    fout.write('# This specifies the resources needed for each task\n')

    fout.write('#PBS -j oe\n')
    fout.write('# Join up output and error into one file\n')

    fout.write('#PBS -J 0-%d\n'%fieldDict[key][1])
    fout.write('# This is the job array request for 100 jobs labelled 1 to 100\n')

    fout.write('# OPTIONS FOR PBS PRO=================================================================\n')
    fout.write('# Change to correct directory\n')

    fout.write('module load python-uon\n')
    fout.write('cd /home/ppztw1/mcmc_runs\n')

    fout.write('# Here we just use Unix command to run our program\n')

    fout.write('echo "Running on `hostname`"\n')

    fout.write('# This script just gives the value of the PBS_ARRAY_INDEX variable\n')

    fout.write('echo  "The TASK ID of *this* task is : $PBS_ARRAY_INDEX"\n')

    fout.write('echo "If I named my input file INPUT.$PBS_ARRAY_INDEX, I could run my program on this file in parallel with all the other input files"\n')
    fout.write('./script_total_aper_field%d.$PBS_ARRAY_INDEX\n'%key)
    fout.write('sleep 20\n')


    fout.write('echo "Finished job now"\n')

    fout.close()
