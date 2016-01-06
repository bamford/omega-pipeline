import numpy as np
from matplotlib import pyplot as plt

from .omega_models import fixha_model, wlHa
from .omega_priors import create_priors
from .omega_init_pars import create_init_pars
from .run_fit import run_fixha, run_line, run_flat, run_fixha_poor


run_emcee_conf = dict(ntemps=0,  # use cunning ladder
                      nwalkers=50,
                      nsamples=200,
                      nupdates=10,
                      nburn=100,
                      minlogbeta=None,
                      xlabel='wavelength',
                      ylabel='flux')


# It would be nice to produce a set of simulations covering reasonable
# parameter space to estimate biases, uncertainties and incompleteness.
# However, not enough time at the moment.


def sim_line(continuum, redshift, fluxHa, fluxNII, xmin, xmax, nx=62, SN=10):
    np.random.seed(54321)
    x = np.linspace(xmin, xmax, nx)
    y_true = fixha_model(x, (continuum, redshift, fluxNII,
                             fluxNII / fluxHa, 0.0))
    sigma_true = continuum / float(SN)
    y_noisy = np.random.normal(y_true, sigma_true)
    return x, y_true, y_noisy, sigma_true


def run():
    truths = (3.0, 0.145, 30.0, 20.0)
    x, y_true, y, yerror = sim_line(*truths, SN=10,
                                    xmin=1.14 * wlHa, xmax=1.16 * wlHa)

    #truths = (1.0, 0.15, 30.0, 10.0)
    #x, y_true, y, yerror = sim_line(*truths, SN=10,
    #                                xmin=1.14 * wlHa, xmax=1.16 * wlHa)

    plt.plot(x, y_true, '-')
    plt.errorbar(x, y, yerror, fmt='o')
    plt.savefig('sim_data.pdf')

    lnpriors, ranges = create_priors(x, y)
    init_pars = create_init_pars(x, y)

    run_fixha_poor(x, y, yerror, lnpriors, init_pars, truths=truths,
                   run_emcee_conf=run_emcee_conf)
    truths = (truths[0], truths[1], truths[3], truths[3] / truths[2], 0.0)
    run_fixha(x, y, yerror, lnpriors, init_pars, truths=truths,
              run_emcee_conf=run_emcee_conf)
    run_line(x, y, yerror, lnpriors, init_pars,
             run_emcee_conf=run_emcee_conf)
    run_flat(x, y, yerror, lnpriors, init_pars,
             run_emcee_conf=run_emcee_conf)


if __name__ == '__main__':
    run()