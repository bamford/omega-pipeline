#!/usr/bin/env python
import numpy as np

#key = field, val[0] = number for aper5, val[1] = number for aper_tot
#fieldDict = {2:[535,534],3:[518,514],5:[587,587],6:[745,744],7:[685,685],8:[617,617],9:[582,582],10:[718,717],11:[794,794],12:[685,685],13:[650,649],14:[645,643],15:[702,702],16:[641,641],17:[662,662],18:[732,732],19:[613,613],20:[614,614],21:[727,725],22:[701,701]}

all_fields, ids5, idstotal = np.loadtxt('../plot_aper/galaxies_per_field_hst.txt',dtype='i',unpack=True)

fieldDict={}
for i in range(len(all_fields)):
    key = all_fields[i]
    fieldDict[key] = [ids5[i], idstotal[i]]

print fieldDict

#make array_mcmc_pt_aper5_spb_field*.sh
for key in fieldDict.keys():
    fout = open('array_mcmc_pt_aper5_spb_field%d.sh'%(key),'w')
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

    fout.write('#PBS -J 0-%d\n'%fieldDict[key][0])
    fout.write('# This is the job array request for 100 jobs labelled 1 to 100\n')

    fout.write('# OPTIONS FOR PBS PRO=================================================================\n')
    fout.write('# Change to correct directory\n')

    fout.write('module load python-uon\n')
    fout.write('cd /home/ppzsb1/projects/OMEGA/mcmc_runs\n')

    fout.write('# Here we just use Unix command to run our program\n')

    fout.write('echo "Running on `hostname`"\n')

    fout.write('# This script just gives the value of the PBS_ARRAY_INDEX variable\n')

    fout.write('echo  "The TASK ID of *this* task is : $PBS_ARRAY_INDEX"\n')

    fout.write('echo "If I named my input file INPUT.$PBS_ARRAY_INDEX, I could run my program on this file in parallel with all the other input files"\n')
    fout.write('./script_aper5_spb_field%d.$PBS_ARRAY_INDEX\n'%key)
    fout.write('sleep 20\n')

    fout.write('echo "Finished job now"\n')

    fout.close()

#make array_mcmc_pt_total_aper_spb_field*.sh
for key in fieldDict.keys():
    fout = open('array_mcmc_pt_total_aper_spb_field%d.sh'%(key),'w')
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
    fout.write('cd /home/ppzsb1/projects/OMEGA/mcmc_runs\n')

    fout.write('# Here we just use Unix command to run our program\n')

    fout.write('echo "Running on `hostname`"\n')

    fout.write('# This script just gives the value of the PBS_ARRAY_INDEX variable\n')

    fout.write('echo  "The TASK ID of *this* task is : $PBS_ARRAY_INDEX"\n')

    fout.write('echo "If I named my input file INPUT.$PBS_ARRAY_INDEX, I could run my program on this file in parallel with all the other input files"\n')
    fout.write('./script_total_aper_spb_field%d.$PBS_ARRAY_INDEX\n'%key)
    fout.write('sleep 20\n')


    fout.write('echo "Finished job now"\n')

    fout.close()
