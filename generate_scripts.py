#!/usr/bin/env python
import numpy as np
#all_fields, aper5, total_aper =
all_fields, ids5, idstotal = np.loadtxt('../../plot_aper/galaxies_per_field_hst.txt',dtype='i',unpack=True)

for ii in range(0,len(all_fields)):
    for kk in range(0,ids5[ii]):
        
        f1 = open('script_aper5_field'+str(int(all_fields[ii]))+'.'+str(kk)+'','w+')
       
        f1.write('%1s%20s\n'%('#','!/home/ppztw1/py-tim/bin/python')) 
        f1.write('%6s %5s % 2s %2s\n'%('import', 'numpy', 'as', 'np'))
        f1.write('%10s%1s%23s%i%4s%1s\n'%('filename0=',"'",'../plot_aper/aper5_ids_hst_F',all_fields[ii],'.txt',"'"))
        f1.write('%30s %1s%1s%1s%1s\n'%('f0=np.loadtxt(filename0,dtype=',"'",'i',"'",')'))  
        f1.write('%7s%i%1s\n'%('glx=f0[',kk,']'))
        f1.write('%5s %3s\n'%('print', 'glx'))
        f1.write('%6s %20s\n'%('import','run_mcmc_pt_aper5_all'))  
        f1.write('%31s%i%1s\n'%('run_mcmc_pt_aper5_all.run_glx(glx,',all_fields[ii],')'))
        f1.close()
    
    for jj in range(0,idstotal[ii]):
        f1 = open('script_total_aper_field'+str(int(all_fields[ii]))+'.'+str(jj)+'','w+')
    
        f1.write('%1s%20s\n'%('#','!/home/ppztw1/py-tim/bin/python')) 
        f1.write('%6s %5s % 2s %2s\n'%('import', 'numpy', 'as', 'np'))
        f1.write('%10s%1s%23s%i%4s%1s\n'%('filename0=',"'",'../plot_aper/total_aper_ids_hst_F',all_fields[ii],'.txt',"'"))
        f1.write('%30s %1s%1s%1s%1s\n'%('f0=np.loadtxt(filename0,dtype=',"'",'i',"'",')'))  
        f1.write('%7s%i%1s\n'%('glx=f0[',jj,']'))
        f1.write('%5s %3s\n'%('print', 'glx'))
        f1.write('%6s %20s\n'%('import','run_mcmc_pt_total_aper_all'))  
        f1.write('%39s%i%1s\n'%('run_mcmc_pt_total_aper_all.run_glx(glx,',all_fields[ii],')'))
        f1.close()
