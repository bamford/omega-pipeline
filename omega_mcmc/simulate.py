import numpy as np
from matplotlib import pyplot as plt

from .convenience import LogLikelihood, run_emcee
from .plots import print_emcee
from .omega_models import *
from .omega_priors import *
from .omega_init_pars import create_init_pars


ntemps = 50
nwalkers = 50
nsamples = 2500
nburn = 500


def sim_line(continuum, redshift, fluxHa, fluxNII, xmin, xmax, nx=50, SN=10):
    np.random.seed(54321)
    x = np.linspace(xmin, xmax, nx)
    y_true = fixha_model(x, (continuum, redshift, fluxHa, fluxNII))
    sigma_true = continuum / float(SN)
    y_noisy = np.random.normal(y_true, sigma_true)
    return x, y_true, y_noisy, sigma_true


def fixha(x, y, yerror, lnpriors, init_pars, truths):
    sampler_fixha = run_emcee(LogLikelihood(fixha_model,
                                            x, y, yerror),
                              lnpriors['fixha_normal'],
                              init_pars['fixha'],
                              ntemps=ntemps, nwalkers=nwalkers,
                              nsamples=nsamples,
                              minlogbeta=-6, nupdates=10,
                              outfile='sim_sampler_fixha')
    print_emcee(sampler_fixha, par_fixha, fixha_model, x, y, yerror,
                nburn=nburn, truths=truths, outfile='sim_fixha.pdf')


def fixha_better(x, y, yerror, lnpriors, init_pars):
    sampler_fixha_better = run_emcee(LogLikelihood(fixha_better_model,
                                                   x, y, yerror),
                                     lnpriors['fixha_better'],
                                     init_pars['fixha_better'],
                                     ntemps=ntemps, nwalkers=nwalkers,
                                     nsamples=nsamples,
                                     minlogbeta=-6, nupdates=10,
                                     outfile='sim_sampler_fixha_better')
    print_emcee(sampler_fixha_better, par_fixha_better, fixha_better_model,
                x, y, yerror, nburn=nburn, outfile='sim_fixha_better.pdf')


def line(x, y, yerror, lnpriors, init_pars):
    sampler_line = run_emcee(LogLikelihood(line_model, x, y, yerror),
                             lnpriors['line_normal'],
                             init_pars['line'],
                             ntemps=ntemps, nwalkers=nwalkers,
                             nsamples=nsamples,
                             minlogbeta=-6, nupdates=10,
                             outfile='sim_sampler_line')
    print_emcee(sampler_line, par_line, line_model, x, y, yerror,
                nburn=nburn, outfile='sim_line.pdf')


def flat(x, y, yerror, lnpriors, init_pars):
    sampler_flat = run_emcee(LogLikelihood(flat_model, x, y, yerror),
                             lnpriors['flat_normal'],
                             init_pars['flat'],
                             ntemps=ntemps, nwalkers=nwalkers,
                             nsamples=nsamples,
                             minlogbeta=-6, nupdates=10,
                             outfile='sim_sampler_flat')
    print_emcee(sampler_flat, par_flat, flat_model, x, y, yerror,
                nburn=nburn, outfile='sim_flat.pdf')


def run():
    truths = (1, 0.15, 30, 10)
    x, y_true, y, yerror = sim_line(*truths, SN=10,
                                    xmin=1.14 * wlHa, xmax=1.16 * wlHa)
    plt.plot(x, y_true, '-')
    plt.errorbar(x, y, yerror, fmt='o')
    plt.savefig('sim_data.pdf')

    lnpriors, ranges = create_priors(x, y)
    init_pars = create_init_pars(x, y)

    fixha(x, y, yerror, lnpriors, init_pars, truths)
    fixha_better(x, y, yerror, lnpriors, init_pars)
    line(x, y, yerror, lnpriors, init_pars)
    flat(x, y, yerror, lnpriors, init_pars)


if __name__ == '__main__':
    run()
