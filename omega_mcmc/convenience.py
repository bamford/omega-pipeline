from __future__ import print_function, division
import sys
import time
import numpy as np
from numpy import ma
from scipy import integrate
from emcee import PTSampler, autocorr
import nestle
from multiprocess import Pool

from .pdfs import norm_logpdf


def flatten_without_burn(sampler, nburn, itemp=0):
    c = sampler.chain
    if c.ndim == 4:
        c = c[itemp]
    c = c[:, nburn:]
    return c.reshape((np.product(c.shape[:-1]), c.shape[-1]))


def weight_without_burn(sampler, nburn, itemp=0):
    c = sampler.lnprobability
    if c.ndim == 3:
        c = c[itemp]
    c = c[:, nburn:]
    w = np.exp(c.reshape(np.product(c.shape)))
    return w / w.sum()


def lnprobability_without_burn(sampler, nburn, itemp=0):
    c = sampler.lnprobability
    if c.ndim == 3:
        c = c[itemp]
    c = c[:, nburn:]
    w = c.reshape(np.product(c.shape))
    return w


def lnlikelihood_without_burn(sampler, nburn, itemp=0):
    c = sampler.lnlikelihood
    if c.ndim == 3:
        c = c[itemp]
    c = c[:, nburn:]
    w = np.exp(c.reshape(np.product(c.shape)))
    return w / w.sum()


def expand_ranges(ranges, fraction):
    ranges = ranges.copy()
    range_widths = ranges[:, 1] - ranges[:, 0]
    ranges[:, 0] -= fraction * range_widths
    ranges[:, 1] += fraction * range_widths
    return ranges


def autocor_checks(sampler, nburn, itemp=0):
    a_exp = sampler.acor[0]
    a_int = np.max([autocorr.integrated_time(sampler.chain[itemp, i, nburn:])
                    for i in range(sampler.chain.shape[1])], 0)
    a_exp = max(a_exp)
    a_int = max(a_int)
    print('A reasonable burn-in should be around '
          '{:d} steps'.format(int(10 * a_exp)))
    print('After burn-in, each chain produces one independent '
          'sample per {:d} steps'.format(int(a_int)))
    return a_exp, a_int


def round_sig(x, sig=1):
    d = sig - int(np.floor(np.log10(x))) - 1
    d = max(0, d)
    return np.round(x, d), d


def summary(samples, par, truths=None):
    mean = samples.mean(0)
    sigma = samples.std(0)
    for i, p in enumerate(par):
        err, dp = round_sig(sigma[i], 1)
        val = round(mean[i], dp)
        post = 'e}' if np.abs(np.log10(np.abs(val))) > 3 else 'f}'
        dps = str(dp) + post
        outstr = ('{:16s} = {:8.' + dps +
                  ' +- {:<8.' + dps).format(p, val, err)
        if truths is not None:
            dpt = str(dp + 1) + post
            outstr += ('   ({:8.' + dpt + ')').format(truths[i])
        print(outstr)
    return mean, sigma


def log_evidence(sampler, nburn):
    logls = sampler.lnlikelihood[:, :, nburn:]
    logls = ma.masked_array(logls, mask=logls == -np.inf)
    mean_logls = logls.mean(axis=-1).mean(axis=-1)
    logZ = -integrate.simps(mean_logls, sampler.betas)
    logZ2 = -integrate.simps(mean_logls[::2], sampler.betas[::2])
    logZerr = abs(logZ2 - logZ)
    print('estimated evidence = {} +- {}'.format(logZ, logZerr))
    return logZ, logZerr


class LogLikelihood:
    def __init__(self, model, x, y, sigma):
        self.model = model
        self.y = y
        self.x = x
        self.sigma = sigma

    def __call__(self, p):
        yM = np.atleast_2d(self.model(self.x, p))
        lnL = norm_logpdf(self.y, yM, self.sigma).sum(axis=1)
        return lnL.squeeze()


def check_init_pars(logl, logp, p0):
    ok = True
    prior0bad = (logp(p0) == -np.inf).sum() / np.product(p0.shape)
    if prior0bad > 0:
        ok = False
        print('Warning: {:.2f}% of initial parameters '
              'have zero prior'.format(prior0bad * 100))
    lnfunc0bad = (logl(p0) == -np.inf).sum() / np.product(p0.shape)
    if lnfunc0bad > 0:
        ok = False
        print('Warning: {:.2f}% of initial parameters '
              'have zero likelihood'.format(lnfunc0bad * 100))
    return ok


def run_emcee(logl, logp, p0func,
              ntemps=50, nwalkers=50, nsamples=2500,
              minlogbeta=-6, nupdates=10,
              threads=1, outfile=None):
    betas = np.logspace(0, minlogbeta, ntemps)
    pos = p0func((ntemps, nwalkers))
    if not check_init_pars(logl, logp, pos):
        return None
    ndim = pos.shape[-1]
    if threads > 1:
        pool = Pool(threads)
    else:
        pool = None
    sampler = PTSampler(ntemps, nwalkers, ndim, logl, logp,
                        betas=betas, pool=pool)
    if nupdates > 0:
        start = time.clock()
        print('Steps:', end='')
        sys.stdout.flush()
        nsteps = nsamples // nupdates
        for i in range(nupdates):
            pos, lnprob, rstate = sampler.run_mcmc(pos, nsteps)
            print(' {}'.format((i + 1) * nsteps), end='')
            sys.stdout.flush()
        print('\nTime taken = {:.2f} secs'.format(time.clock() - start))
    nsteps = nsamples - sampler.chain.shape[-2]
    if nsteps > 0:
        pos, lnprob, rstate = sampler.run_mcmc(pos, nsteps)
    if outfile is None:
        outfile = 'emcee_sampler.npz'
    np.savez(outfile,
             acceptance_fraction=sampler.acceptance_fraction,
             acor=sampler.acor, beta=sampler.betas, chain=sampler.chain,
             lnlikelihood=sampler.lnlikelihood,
             lnprobability=sampler.lnprobability,
             tswap_acceptance_fraction=sampler.tswap_acceptance_fraction)
    return sampler


def run_nestle(logl, logp, ranges, npoints=1000, method='single'):
    def logprob(p):
        return logl(p) + logp(p)

    def prior_transform(p):
        return (ranges[:, 1] - ranges[:, 0]) * p + ranges[:, 0]

    return nestle.sample(logprob, prior_transform, len(ranges),
                         method=method, npoints=npoints)
