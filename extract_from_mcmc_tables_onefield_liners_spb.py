#!/usr/bin/env /home/project/OMEGA/py-omega/bin/python
### In this script we take the information from the tables generated in the MCMC runs.
### A catalogue will be created with the median, maximum-likelihood values and errors of the
### four parameters used in our model. In the cases where there is no measurement usign the
### total aperture we copy the values obtained for the aperture of five pixels.

import pyfits
import numpy as np
from scipy import stats
import os
import math

def open_fits(aperture,glx,field):

    ### Open table aper 5
    hdu_list = pyfits.open('galaxy_pages_'+str(aperture)+'/F'+str(field)+'/table_'+str(glx)+'.fits')
    #hdu_list = pyfits.open('galaxy_pages_'+str(aperture)+'/plot_aper_web/table_'+str(glx)+'.fits')


    ### Accesing the data table:
    table_hdu = hdu_list[1]
    ### Header of the table:
    table_header = table_hdu.header
    ### To see what columns are in the table:
    # table_hdu.columns
    ### Extracting the data from the table:
    table_data = table_hdu.data
    ### The first column of the table would be:
    par0 = table_data.field(0)
    ### or:
    continuum = np.array(table_data.field('continuum'))
    redshift = np.array(table_data.field('redshift'))
    flux_ha = np.array(table_data.field('flux Ha'))
    flux_nii = np.array(table_data.field('flux NII'))
    lnprob = np.array(table_data.field('Probability'))

    max_prob = np.where(lnprob == max(lnprob))[0][0]

    prob_posi_ew = stats.percentileofscore(flux_ha/continuum,0,'weak')

    prob_agn_lineratio = stats.percentileofscore(flux_nii/flux_ha,10**(-0.4),'weak')
    prob_agn_halpha = stats.percentileofscore(flux_ha/continuum,6,'weak')

    prob_sf = 0.01*(100-stats.percentileofscore(flux_ha/continuum,3,'weak'))*stats.percentileofscore(flux_nii/flux_ha,10**(-0.4),'weak')
    prob_strong_agn = 0.01*(100-stats.percentileofscore(flux_ha/continuum,6,'weak'))*(100-stats.percentileofscore(flux_nii/flux_ha,10**(-0.4),'weak'))
    prob_retired = 0.01*(stats.percentileofscore(flux_ha/continuum,3,'weak'))

    prob_weak_agn = 0.01*0.01*(100-stats.percentileofscore(flux_ha/continuum,3,'weak'))*(stats.percentileofscore(flux_ha/continuum,6,'weak'))*(100-stats.percentileofscore(flux_nii/flux_ha,10**(-0.4),'weak'))

    cont_lowper = stats.scoreatpercentile(continuum,16)
    cont_highper = stats.scoreatpercentile(continuum,84)
    cont_maxprob = continuum[max_prob]
    cont_median = np.median(continuum)

    z_lowper = stats.scoreatpercentile(redshift,16)
    z_highper = stats.scoreatpercentile(redshift,84)
    z_maxprob = redshift[max_prob]
    z_median = np.median(redshift)

    fha_lowper = stats.scoreatpercentile(flux_ha,16)
    fha_highper = stats.scoreatpercentile(flux_ha,84)
    fha_maxprob = flux_ha[max_prob]
    fha_median = np.median(flux_ha)

    lha_median = fha_median*4*math.pi*(2.44*10**27)**2

    fnii_lowper = stats.scoreatpercentile(flux_nii,16)
    fnii_highper = stats.scoreatpercentile(flux_nii,84)
    fnii_maxprob = flux_nii[max_prob]
    fnii_median = np.median(flux_nii)

    rel_error_cont_median = (cont_highper - cont_lowper)/(2*cont_median)
    rel_error_z_median = (z_highper - z_lowper)/2
    rel_error_fha_median = (fha_highper - fha_lowper)/(2*fha_median)
    rel_error_fnii_median = (fnii_highper - fnii_lowper)/(2*fnii_median)

    rel_error_cont_maxlik = (cont_highper - cont_lowper)/(2*cont_maxprob)
    rel_error_z_maxlik = (z_highper - z_lowper)/2
    rel_error_fha_maxlik = (fha_highper - fha_lowper)/(2*fha_maxprob)
    rel_error_fnii_maxlik = (fnii_highper - fnii_lowper)/(2*fnii_maxprob)

    ewha_median = fha_median/cont_median
    error_ewha_median = math.sqrt((rel_error_fha_median*fha_median/cont_median)**2 + (fha_median*rel_error_cont_median*cont_median/(cont_median**2))**2)

    lineratio_median = fnii_median/fha_median
    error_lineratio_median = math.sqrt((rel_error_fnii_median*fnii_median/fha_median)**2 + (fnii_median*rel_error_fha_median*fha_median/(fha_median**2))**2)

    ewha_maxlik = fha_maxprob/cont_maxprob
    error_ewha_maxlik = math.sqrt((rel_error_fha_maxlik*fha_maxprob/cont_maxprob)**2 + (fha_maxprob*rel_error_cont_maxlik*cont_maxprob/(cont_maxprob**2))**2)

    lineratio_maxlik = fnii_maxprob/fha_maxprob
    error_lineratio_maxlik = math.sqrt((rel_error_fnii_maxlik*fnii_maxprob/fha_maxprob)**2 + (fnii_maxprob*rel_error_fha_maxlik*fha_maxprob/(fha_maxprob**2))**2)


    ### Number of points to calculate reduced chi-square
    #filename='../plot_aper/F'+str(field)+'/F'+str(field)+'_spectrum_gal_'+str(glx)+'_'+str(aperture)+'.txt'
    filename='../plot_aper/all_spectra_hst/F'+str(field)+'_spectrum_gal_'+str(glx)+'_'+str(aperture)+'_hst.txt'

    #filename='../plot_aper/plot_aper_web/spectrum_gal_'+str(glx)+'_'+str(aperture)+'.txt'

    f = np.loadtxt(filename, dtype='f')
    #xmeans, ymeans, error_flux_phot, error_flux_zeropointfit, error_flux = f[:,:5].T
    xmeans, ymeans, error_flux = f[:,:3].T
    lambdamin = xmeans.min()
    lambdamax = xmeans.max()
    ok = np.logical_not(np.isnan(ymeans) | np.isnan(error_flux))
    degs = len(ok)

    ##chi-squared
    chisq =  (-2*lnprob[max_prob])/degs

    ## are the two lines inside the wavelength range?
    if ((1+z_median)*6583.5 < lambdamax - 3.75)&((1+z_median)*6562.8 > lambdamin +3.75):
       twolines = 1
    else:
       twolines = 0
    if ((1+z_median)*6562.8 < lambdamax - 3.75)&((1+z_median)*6562.8 > lambdamin +3.75):
       oneline = 1
    else:
       oneline = 0





    return  int(glx), 100 - prob_posi_ew,prob_strong_agn,prob_sf,prob_retired,prob_weak_agn, cont_maxprob, cont_median,rel_error_cont_median,rel_error_cont_maxlik, z_maxprob, z_median,rel_error_z_median, fha_maxprob, fha_median,rel_error_fha_median, rel_error_fha_maxlik,lha_median, fnii_maxprob, fnii_median,rel_error_fnii_median,rel_error_fnii_maxlik, ewha_median,ewha_maxlik,error_ewha_median,error_ewha_maxlik,lineratio_median,lineratio_maxlik,error_lineratio_median,error_lineratio_maxlik, chisq,twolines,oneline
formato = ('%5i %.5f %.5f %.5f %.5f %.5f %.5g %.5g %.5g %.5g %.5g %.5g %.5g  %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g  %.5g %.5g %5i %.5f %.5f %.5f %.5f %.5f %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g  %.5g %.5g')
#print len(formato)


def in_field(field):

#def all():
    array = np.loadtxt('../plot_aper/aper5_ids_hst_F'+str(field)+'.txt', dtype='i')
    #array = np.loadtxt('../plot_aper/all_galaxies_comboflag3_aper5_python.txt', dtype='i')

    #for ii in range(0,1):
    #    each = array[ii]
    for each in array:
        ### We try to open the tables of aper5 and total_aperture.
        try:
            f5 = open_fits('aper5_spb',each,field)
            print 'opened fits aper5_spb'
            two = np.ones([1,2*len(f5)])
            ### We set to '-99' all the values for which we don't have a measurement.
            two = -99*two
            try:
                open_fits('total_aper_spb',each,field)
                print 'opened fits total_aper_spb'
                ft = open_fits('total_aper_spb',each,field)
                two[:,0:len(f5)] = f5
                two[:,len(f5):]=ft
                f_handle = file('liners_catalogues_spb/F'+str(field)+'_catalogue_hst_liners.txt','a')
                np.savetxt(f_handle,two,fmt=formato)
                f_handle.close()
            except IOError:
                print 'IOError: Not in total aper',each
                two[:,0:len(f5)] = f5
                ### Because there is no measurement using total aperture, we copy the values from the 5pixel aperture
                two[:,len(f5):]=f5
                f_handle = file('liners_catalogues_spb/F'+str(field)+'_catalogue_hst_liners.txt','a')
                np.savetxt(f_handle,two,fmt=formato)
                f_handle.close()
            except IndexError:
                print 'IndexError: Not in total aper',each
                two[:,0:len(f5)] = f5
                ### Because there is no measurement using total aperture, we copy the values from the 5pixel aperture
                two[:,len(f5):]=f5
                f_handle = file('liners_catalogues_spb/F'+str(field)+'_catalogue_hst_liners.txt','a')
                np.savetxt(f_handle,two,fmt=formato)
                f_handle.close()

        except IOError:
            print 'IOError: Not in aper 5',each

            try:
                open_fits('total_aper_spb',each,field)
                print 'opened fits total_aper_spb'
                ft = open_fits('total_aper_spb',each,field)
                two = np.ones([1,2*len(ft)])
                two = -99*two
                two[:,len(ft):]=ft
                f_handle = file('liners_catalogues_spb/F'+str(field)+'_catalogue_hst_liners.txt','a')
                np.savetxt(f_handle,two,fmt=formato)
                f_handle.close()

            except IOError:
                print 'IOError: Not in any aperture', each
            except IndexError:
                print 'IndexError: Not in any aperture', each
