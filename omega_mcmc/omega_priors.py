import numpy as np
from scipy import stats

from .pdfs import norm_logpdf, gamma_logpdf, beta_logpdf
from .omega_models import wlHa, linewidth, get_ranges


# Purposefully-informative priors
# These are defined independently for each parameter, so they can be combined
# while maintaining normalisation.


def lnprior_continuum(continuum, ymin, ymax):
    # proper prior on continuum
    # same as flat model normal prior
    return norm_logpdf(continuum, 0.5 * ymin, 0.25 * (ymax - ymin))


def lnprior_fluxNII(fluxNII, fmax):
    # proper prior on flux of NII line: flat with cutoff at high fluxes
    # flux > fmax disfavoured by lnL = -1
    # flux > 10*fmax disfavoured by lnL = -10
    # and positive definite
    return gamma_logpdf(fluxNII, a=1.0, scale=fmax)


def lnprior_NIIHa(NIIHa):
    # proper prior on NII/Ha ratio: peaked at ~0.5
    # fairly flat, with cutoffs
    # NII/Ha < 0.1 and > 3 disfavoured by lnL = -2.7
    # and positive definite
    return gamma_logpdf(NIIHa, a=2.7, scale=0.5)


def lnprior_absEWHa(absEWHa):
    # proper prior on EW of Ha absorption: peaked at ~0.5
    # fairly flat, with cutoffs
    # absEWHa > 3 disfavoured by lnL = -3
    # and positive definite
    return gamma_logpdf(absEWHa, a=1.0, scale=1.0)


def lnprior_redshift(redshift, zmin, zmax):
    # Set edges of redshift box (prior > 0.9) such that
    # peak of Halpha is in the wavelength range
    # Set cutoffwidth such that lnL ~ -2 by point at which
    # Halpha is out of wavelength range by 3 sigma
    c = 1.4
    return beta_logpdf(redshift, c, c, loc=zmin, scale=(zmax - zmin))


class LogPriorFixHa:
    def __init__(self, ymin, ymax, zmin, zmax, fmax):
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.fmax = fmax

    def __call__(self, p):
        p = np.atleast_2d(p)
        continuum, redshift, fluxNII, NIIHa, absEWHa = p.T
        lnprior = lnprior_continuum(continuum, self.ymin, self.ymax)
        lnprior += lnprior_redshift(redshift, self.zmin, self.zmax)
        lnprior += lnprior_fluxNII(fluxNII, self.fmax)
        lnprior += lnprior_NIIHa(NIIHa)
        lnprior += lnprior_absEWHa(absEWHa)
        return lnprior.squeeze()


# Less-informative priors
# Simpler priors for comparison and testing


class LogPriorUniform:
    def __init__(self, ranges):
        self.rv = stats.uniform(ranges[:, 0], ranges[:, 1] - ranges[:, 0])

    def __call__(self, x):
        x = np.atleast_2d(x)
        lnP = self.rv.logpdf(x).sum(axis=1)
        return lnP.squeeze()


class LogPriorNormal:
    def __init__(self, ranges, trunc=False):
            mean = 0.5 * (ranges[:, 1] + ranges[:, 0])
            sigma = 0.25 * (ranges[:, 1] - ranges[:, 0])
            if trunc:
                self.rv = stats.truncnorm(-2, 2, mean, sigma)
            else:
                self.rv = stats.norm(mean, sigma)

    def __call__(self, x):
        x = np.atleast_2d(x)
        lnP = self.rv.logpdf(x).sum(axis=1)
        return lnP.squeeze()


def create_priors(x, y):
    xmin, xmax, ymin, ymax, zmin, zmax, fmin, fmax = get_ranges(x, y)

    lnpriors, ranges = {}, {}

    lnpriors['fixha'] = LogPriorFixHa(ymin, ymax, zmin, zmax, fmax)

    # ranges for fixha are just for nestle
    ranges['fixha'] = np.array([[ymin, ymax], [zmin, zmax],
                                       [fmin, fmax], [-1, 10], [-1, 10]])

    ranges['fixha_poor'] = np.array([[ymin, ymax], [zmin, zmax],
                             [fmin, fmax], [fmin, fmax]])
    lnpriors['fixha_poor_uniform'] = LogPriorUniform(ranges['fixha_poor'])
    lnpriors['fixha_poor_normal'] = LogPriorNormal(ranges['fixha_poor'])

    ranges['line'] = ranges['fixha_poor'][0:3]
    lnpriors['line_uniform'] = LogPriorUniform(ranges['line'])
    lnpriors['line_normal'] = LogPriorNormal(ranges['line'])
    lnpriors['line_truncnormal'] = LogPriorNormal(ranges['line'], trunc=True)

    ranges['flat'] = ranges['fixha_poor'][0:1]
    lnpriors['flat_uniform'] = LogPriorUniform(ranges['flat'])
    lnpriors['flat_normal'] = LogPriorNormal(ranges['flat'])

    return lnpriors, ranges
