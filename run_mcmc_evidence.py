### Running the mcmc with parallel tempered

import numpy as np
from astropy.io import fits as pyfits

from omega_mcmc.convenience import (flatten_without_burn,
                                    lnprobability_without_burn)
from omega_mcmc.omega_priors import create_priors
from omega_mcmc.omega_init_pars import create_init_pars
from omega_mcmc.run_fit import run_fixha_sigma, run_line_sigma, run_flat_sigma


run_emcee_conf = dict(ntemps=0,  # 0 to use cunning ladder
                      nwalkers=50,
                      nsamples=2000,
                      nupdates=10,
                      nburn=1000,
                      minlogbeta=None,
                      xlabel='wavelength',
                      ylabel='flux')

output_path = '../mcmc_fits/galaxy_pages_{aper}_spb/F{field}/'

fluxunit = 1e-16


def get_data(glx, field, aper):
    filetemplate = ('../plot_aper/all_spectra_hst/'
                    'F{field}_spectrum_gal_{glx}_{aper}_hst.txt')
    filename = filetemplate.format(glx=glx, field=field, aper=aper)

    f = np.loadtxt(filename, dtype='g')
    xmeans, ymeans, yerror = f[:, :3].T

    # Getting rid of the nan values...
    ok = np.logical_not(np.isnan(ymeans) | np.isnan(yerror))

    yerror[yerror == 0] = 1.0  # ARBITRARY?!?

    xmeans = xmeans[ok].astype(np.float64)
    ymeans = ymeans[ok].astype(np.float64)
    yerror = yerror[ok].astype(np.float64)

    ymeans /= fluxunit
    yerror /= fluxunit
    return xmeans, ymeans, yerror


def run_glx(glx, field, aper, reprocess=False):
    x, y, yerror = get_data(glx, field, aper)

    lnpriors, ranges = create_priors(x, y)
    init_pars = create_init_pars(x, y)

    sampler_fixha, stats_fixha = run_fixha_sigma(
        x, y, yerror, lnpriors, init_pars,
        run_emcee_conf=run_emcee_conf,
        label=(output_path + '{glx}').format(
            aper=aper, field=field, glx=glx),
        reprocess=reprocess)

    run_emcee_conf_line = run_emcee_conf.copy()
    run_emcee_conf_line['nburn'] /= 2
    run_emcee_conf_line['nsamples'] = run_emcee_conf_line['nburn'] + 1000
    sampler_line, stats_line = run_line_sigma(
        x, y, yerror, lnpriors, init_pars,
        run_emcee_conf=run_emcee_conf_line,
        label=(output_path + '{glx}').format(
            aper=aper, field=field, glx=glx),
        reprocess=reprocess)

    run_emcee_conf_flat = run_emcee_conf.copy()
    run_emcee_conf_flat['nburn'] /= 4
    run_emcee_conf_flat['nsamples'] = run_emcee_conf_flat['nburn'] + 1000
    sampler_flat, stats_flat = run_flat_sigma(
        x, y, yerror, lnpriors, init_pars,
        run_emcee_conf=run_emcee_conf_flat,
        label=(output_path + '{glx}').format(
            aper=aper, field=field, glx=glx),
        reprocess=reprocess)

    nburn = run_emcee_conf['nburn']
    flattable = flatten_without_burn(sampler_fixha, nburn)
    prob = lnprobability_without_burn(sampler_fixha, nburn)

    continuum = flattable[:, 0]
    redshift = flattable[:, 1]
    fluxHa = flattable[:, 2] / flattable[:, 3]
    fluxNII = flattable[:, 2]
    NIIHa = flattable[:, 3]
    absEWHa = flattable[:, 4]
    sconst = flattable[:, 5]
    sfactor = flattable[:, 6]
    badfrac = flattable[:, 7]

    continuum *= fluxunit
    fluxHa *= fluxunit
    fluxNII *= fluxunit

    # Creating the FITS file with the samplers and the probability
    cols = [pyfits.Column(name='continuum', format='E', array=continuum),
            pyfits.Column(name='redshift', format='E', array=redshift),
            pyfits.Column(name='flux Ha', format='E', array=fluxHa),
            pyfits.Column(name='flux NII', format='E', array=fluxNII),
            pyfits.Column(name='NIIHa', format='E', array=NIIHa),
            pyfits.Column(name='absEWHa', format='E', array=absEWHa),
            pyfits.Column(name='sconst', format='E', array=sconst),
            pyfits.Column(name='sfactor', format='E', array=sfactor),
            pyfits.Column(name='badfrac', format='E', array=badfrac),
            pyfits.Column(name='Probability', format='E', array=prob)]

    tbhdu = pyfits.new_table(pyfits.ColDefs(cols))
    phdu = pyfits.PrimaryHDU()
    phdu.header['A_EXP'] = stats_fixha['a_exp']
    phdu.header['A_INT'] = stats_fixha['a_int']
    phdu.header['ACC_FRAC'] = stats_fixha['acc_frac']
    phdu.header['Z_FIXHA'] = stats_fixha['logz']
    phdu.header['ZE_FIXHA'] = stats_fixha['logzerr']
    phdu.header['Z_LINE'] = stats_line['logz']
    phdu.header['ZE_LINE'] = stats_line['logzerr']
    phdu.header['Z_FLAT'] = stats_flat['logz']
    phdu.header['ZE_FLAT'] = stats_flat['logzerr']
    phdu.header['ITEMP50'] = stats_fixha['itemp50']
    phdu.header['ITEMP90'] = stats_fixha['itemp90']
    hdulist = pyfits.HDUList([phdu, tbhdu])
    hdulist.writeto((output_path + 'table_{glx}.fits').format(
        aper=aper, field=field, glx=glx))
