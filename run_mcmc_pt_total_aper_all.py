### Running the mcmc with parallel tempered

import numpy as np
import emcee
from math import sqrt, pi,log10
from matplotlib.ticker import MaxNLocator
from scipy.integrate import quad
import math
from scipy import stats
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from emcee import PTSampler
#from astropy.io import fits
import pyfits
import triangle

wlHa=6562.8
wlNIIa =6548.1
wlNIIb=6583.5

par = ['continuum', 'redshift', 'flux Ha', 'flux NII']
ndim_fixha = 4

nwalkers = 300
nburn = 500
nsamp = 500


def gaussian(x, ampl, centre, sigma):
    # with unit maximum
    return (ampl/(sigma*sqrt(2*pi))) * np.exp(-(x-centre)**2/(2.0*sigma**2))

def fixha_model(x, p):
    model = p[0] + gaussian(x, p[2], (1+p[1])*wlHa, 7.5)
    model += gaussian(x, p[3]/3.06, (1+p[1])*wlNIIa, 7.5) + gaussian(x, p[3], (1+p[1])*wlNIIb, 7.5)
    return model
  
def lnprob_fixha(p, y, x, icov, max_fha, zmin, zmax):
    # Here we set the hard priors:
    # Prior limiting the absorption in Ha
    if (p[3] < 0) or (p[1] > zmax) or (p[1] < zmin) or ((p[3]/p[2]) > 3) or (p[3] > (3*p[2] + 3*p[0])):

    # No prior allowing a 3A absorption in Ha 
    #if (p[3] < 0) or (p[1] > zmax) or (p[1] < zmin) or ((p[3]/p[2]) > 3):
        return -np.inf
    model = fixha_model(x, p)
    chisq = ((model - y)**2 * icov).sum()
    return -chisq/2.0
  
def get_data(glx,field):

    #filename='../plot_aper/plot_aper_web/spectrum_gal_'+str(glx)+'_total_aper.txt'
    #filename='../plot_aper/F'+str(field)+'/F'+str(field)+'_spectrum_gal_'+str(glx)+'_total_aper.txt'
    filename='../plot_aper/all_spectra_hst/F'+str(field)+'_spectrum_gal_'+str(glx)+'_total_aper_hst.txt'

    f = np.loadtxt(filename, dtype='g')
    #xmeans, ymeans, error_flux_phot, error_flux_zeropointfit, error_flux = f[:,:5].T
    xmeans, ymeans, error_flux = f[:,:3].T

    ## Getting rid of the nan values...
    ok = np.logical_not(np.isnan(ymeans) | np.isnan(error_flux))
    
    error_flux[error_flux == 0] = 1.0  # ARBITRARY?!?
        
    xmeans = xmeans[ok]
    ymeans = ymeans[ok]
    error_flux = error_flux[ok]
    xmin = np.min(xmeans)
    xmax = np.max(xmeans)
    ymin = np.min(ymeans)
    ymax = np.max(ymeans)
    icov = 1/error_flux**2
    zmax = (xmax+20)/wlHa - 1
    zmin = (xmin-20)/wlHa - 1
    return xmeans, ymeans, error_flux, icov, xmin, xmax, ymin, ymax, zmin, zmax   
  
def init_p0_fixha(ymin, ymax,zmin, zmax):
    # very broad uniform initial distributions
    np.random.seed(1234)### Establish the same random numbers every time we run it!!
    c0_init = np.random.uniform(ymin, ymax, nwalkers)
    z_init = np.random.uniform(zmin, zmax, nwalkers)
    fha_init = (np.random.uniform(ymin, ymax, nwalkers) - ymin) * sqrt(2*pi) * 7.5
    fnii_init = (np.random.uniform(ymin, ymax, nwalkers) - ymin) * sqrt(2*pi) * 7.5
    return np.transpose([c0_init, z_init, fha_init, fnii_init])

def flatten_without_burn(sampler, nburn):
    c = sampler.chain
    if c.ndim == 4:
        c = c[0]
    c = c[:, nburn:]
    return c.reshape((np.product(c.shape[:-1]), c.shape[-1]))


def weight_without_burn(sampler, nburn):
    c = sampler.lnprobability
    if c.ndim == 3:
        c = c[0]
    c = c[:, nburn:]
    w = np.exp(c.reshape(np.product(c.shape)))
    return w / w.sum() 


  
def run_glx(glx,field): 

    xmeans, ymeans, error_flux, icov, xmin, xmax, ymin, ymax, zmin, zmax = get_data(glx,field)

    p0_fixha = init_p0_fixha(ymin, ymax, zmin, zmax)

    # value doesn't actually matter here
    max_fha = ymax*sqrt(2*pi)*7

    def logl(x):

        return lnprob_fixha(x, ymeans, xmeans, icov, max_fha, zmin, zmax)

    # Use a flat prior
    def logp(x):
        return 0.0

    ntemps = 10  
    sampler_fixha_pt = PTSampler(ntemps, nwalkers, ndim_fixha, logl, logp)
    r_pt = sampler_fixha_pt.run_mcmc(p0_fixha*np.ones((ntemps, nwalkers, ndim_fixha)), nburn+nsamp)
    
######## AUTOCORRELATION CHECKS!!!!!! #####################
    a_exp = sampler_fixha_pt.acor[0]
    print a_exp
    a_int = np.max([emcee.autocorr.integrated_time(sampler_fixha_pt.chain[0, i, nburn:]) for i in range(len(sampler_fixha_pt.chain[0]))], 0)
    print a_int
    a_exp = max(a_exp)
    a_int = max(a_int)
    print('A reasonable burn-in should be around {:d} steps'.format(int(10*a_exp)))
    print('After burn-in, each chain produces one independent sample per {:d} steps'.format(int(a_int)))

    g = [emcee.autocorr.function(sampler_fixha_pt.chain[0, i, :nburn])
         for i in range(3)]

    ax = plt.figure().add_subplot(111)
    for i in range(3):
        ax.plot(g[i][:200], "k")
    ax.axhline(0, color="k")
    ax.set_xlim(0, 200)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"Autocorrelation")
    #ax.figure.savefig("acor.png", transparent=True, dpi=300)
    plt.savefig('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/F'+str(field)+'_mcmc_burnin_acor_'+str(glx)+'.pdf')
    
   ### Autocorrelation for AFTER the burn-in

    f = [emcee.autocorr.function(sampler_fixha_pt.chain[0, i, nburn:])
         for i in range(3)]

    ax = plt.figure().add_subplot(111)
    for i in range(3):
        ax.plot(f[i][:200], "k")
    ax.axhline(0, color="k")
    ax.set_xlim(0, 200)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"Autocorrelation")
    #ax.figure.savefig("acor.png", transparent=True, dpi=300)
    plt.savefig('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/F'+str(field)+'_mcmc_acor_'+str(glx)+'.pdf')

    ### gk_i es la suma de todos los valores anteriores del parametro gk en la cadena i dividido entre el numero de valores. Vamos, un valor medio.
    
    ### gk_twodots es la suma de todos los valores anteriores de gk en todas las cadenas dividido entre numero total de cadenas multiplicado por el numero de iterations.

    all_R_k = []
    gk_all = []    
    B_over_q_each = [] 
    W_each = []

    for qq in range(0,nsamp):
        for ll in range(0,nwalkers):

            gk_ll = sum(sampler_fixha_pt.chain[0,ll,:qq,0])/nsamp
            gk_twodots = sum(sum(sampler_fixha_pt.chain[0,:,:qq,0]))/(nsamp*nwalkers)  
	    
	    gk_ll_qq = sampler_fixha_pt.chain[0,ll,qq,0]
            W_each.append((gk_ll_qq - gk_ll)**2)
            
            B_over_q_each.append(((gk_ll - gk_twodots)**2))
        
        B_over_q_total = sum(B_over_q_each)/(nwalkers-1)        
        W_total = sum(W_each)/(nwalkers*(qq - 1))
        sigma_2_plus = (qq - 1)*W_total/qq + B_over_q_total
        
        all_R_k.append((sigma_2_plus + B_over_q_total/nwalkers)/W_total)

    f1 = open('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/F'+str(field)+'_gal_'+str(glx)+'_convergence_test.txt','a+')
    for hh in range(0,len(all_R_k)):
        f1.write('%5i %.3f\n'%(hh,all_R_k[hh]))
    f1.close()
    
    ax = plt.figure().add_subplot(111)
    ax.plot(all_R_k)
    plt.savefig('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/F'+str(field)+'_gal_'+str(glx)+'_convergence_test.pdf')


##################### END OF AUTOCORRELATION CHECKS!!!!!!! #################################    

    f1 = open('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/F'+str(field)+'_acceptance_fractions.txt','a+')
    f1.write('%5i %.3f\n'%(glx,np.mean(sampler_fixha_pt.acceptance_fraction)))
    f1.close()
      
    flattable = flatten_without_burn(sampler_fixha_pt, nburn)

    prob = sampler_fixha_pt.lnprobability[0,:, nburn:]
        
    prob_reshaped = prob.reshape((np.product(prob.shape[:]), 1))

    likel = sampler_fixha_pt.lnlikelihood[0,:, nburn:]
    likel_reshaped = likel.reshape((np.product(likel.shape[:]), 1))

   # print 'prob', -2*prob_reshaped

    
    a0 = flattable[:,0]
    a1 = flattable[:,1]
    a2 = flattable[:,2]
    a3 = flattable[:,3]
    a4 = prob_reshaped[:,0]

    
    ### Creating the FITS file with the samplers and the probability
    col1 = pyfits.Column(name='continuum', format='E', array=a0)
    col2 = pyfits.Column(name='redshift', format='E', array=a1)
    col3 = pyfits.Column(name='flux Ha', format='E', array=a2)
    col4 = pyfits.Column(name='flux NII', format='E', array=a3)
    col5 = pyfits.Column(name='Probability', format='E', array=a4)
    
    cols = pyfits.ColDefs([col1, col2, col3, col4, col5])
    
    tbhdu = pyfits.new_table(cols)
    
    phdu = pyfits.PrimaryHDU()
    
    hdulist = pyfits.HDUList([phdu, tbhdu])
    
    hdulist.writeto('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/table_'+str(glx)+'.fits')
    #hdulist.writeto('table_'+str(glx)+'.fits')
    '''

    plt.figure(figsize=[20,10])
    for i, p in enumerate(par):
        plt.subplot(2,2,i+1)
        for t in [0]: #range(ntemps):
            for w in range(nwalkers):
                plt.plot(np.arange(len(sampler_fixha_pt.chain[t,0,:,0])), sampler_fixha_pt.chain[t,w,:,i], 'r-', alpha=0.1)
        plt.xlabel(p)
        aymin, aymax = plt.ylim()
        plt.vlines(nburn, aymin, aymax, linestyle=':')
        plt.ylim(aymin, aymax)
        plt.savefig('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/mcmc_runs'+str(glx)+'.pdf')
        #plt.savefig('mcmc_runs'+str(glx)+'.pdf')
    
    plt.figure(figsize=[20,10])
    for i, p in enumerate(par):
        plt.subplot(2,2,i+1)
        plt.hist(flatten_without_burn(sampler_fixha_pt, nburn)[:,i], bins=100, histtype='stepfilled', label='PT', alpha=0.75)
        plt.xlabel(p)
    plt.legend()    
    plt.savefig('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/mcmc_runs_histograms_'+str(glx)+'.pdf') 
    plt.close()
    #plt.savefig('mcmc_runs_histograms_'+str(glx)+'.pdf')    
    
    xchain = np.arange(xmin, xmax, 1.0)
    #ychain = [fixha_model(xchain, p) for p in flatten_without_burn(sampler_fixha, nburn)]
    ychain_pt = [fixha_model(xchain, p) for p in flatten_without_burn(sampler_fixha_pt, nburn)]
    
    
    plt.figure(figsize=[20,10])
    plt.errorbar(xmeans, ymeans, error_flux, None, 'o')
    for i, y in enumerate(ychain_pt[::100]):
        plt.plot(xchain, y, 'g-', alpha=0.1)
    plt.savefig('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/mcmc_runs_fits_'+str(glx)+'.pdf')    
    #plt.savefig('mcmc_runs_fits_'+str(glx)+'.pdf')    

    
    #plt.figure(figsize=[20,10])
    triangle.corner(flatten_without_burn(sampler_fixha_pt, nburn), labels=par)
    gs1 = gridspec.GridSpec(12, 12)
    gs1.update(left=0.08, right=0.90,wspace=0.15)
    plt.subplots_adjust(wspace = 0.1,hspace=0.1,top=0.96,bottom=0.1,left=0.1,right=0.95)

    ax = plt.subplot2grid((12,12), (0,7), colspan=5,rowspan=5)
    for y in ychain_pt[::100]:
        plt.plot(xchain, y, 'r-', alpha=0.01)
    plt.errorbar(xmeans, ymeans, error_flux, None, 'o')
    plt.savefig('../mcmc_fits/galaxy_pages_total_aper/F'+str(field)+'/mcmc_runs_triangle_'+str(glx)+'.pdf')    
        '''

    
   
    

 
