import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

from .pdfs import norm_logpdf, gamma_logpdf, beta_logpdf
from .plots import plot_prior
from .omega_models import get_ranges


# Purposefully-informative priors
# These are defined independently for each parameter, so they can be combined
# while maintaining normalisation.


def lnprior_continuum(continuum, ymin, ymax):
    # proper prior on continuum
    # same as flat model normal prior
    return norm_logpdf(continuum, 0.5 * (ymax - ymin), 0.5 * (ymax - ymin))


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
    # proper prior on EW of Ha absorption:
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


def lnprior_sconst(sconst):
    # proper prior on constant error
    # sconst > x disfavoured by lnL = -x
    # and positive definite
    return gamma_logpdf(sconst, a=1.0, scale=1.0)


def lnprior_sfactor(sfactor):
    # proper prior on error factor: peaked at ~0.5
    # fairly flat, with cutoffs
    # sfactor < 0.1 and > 3 disfavoured by lnL = -2.7
    # and positive definite
    return gamma_logpdf(sfactor, a=2.7, scale=0.5)


def lnprior_badfrac(badfrac):
    # proper prior on bad fraction
    # Constrained to range [0, 1],
    # but skewed to low values
    return beta_logpdf(badfrac, 0.5, 4.0)


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


class LogPriorLine:
    def __init__(self, ymin, ymax, zmin, zmax, fmax):
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.fmax = fmax

    def __call__(self, p):
        p = np.atleast_2d(p)
        continuum, redshift, flux = p.T
        lnprior = lnprior_continuum(continuum, self.ymin, self.ymax)
        lnprior += lnprior_redshift(redshift, self.zmin, self.zmax)
        lnprior += lnprior_fluxNII(flux, self.fmax)
        return lnprior.squeeze()


class LogPriorFlat:
    def __init__(self, ymin, ymax):
        self.ymin = ymin
        self.ymax = ymax

    def __call__(self, p):
        p = np.atleast_2d(p)
        continuum = p.T
        lnprior = lnprior_continuum(continuum, self.ymin, self.ymax)
        return lnprior.squeeze()


# Priors incorporating doubt about errors

class LogPriorSigma():
     def __call__(self, p):
        p = np.atleast_2d(p)
        sconst, sfactor, badfrac = p.T
        lnprior = lnprior_sconst(sconst)
        lnprior += lnprior_sfactor(sfactor)
        lnprior += lnprior_badfrac(badfrac)
        return lnprior.squeeze()


class LogPriorFixHaSigma(LogPriorFixHa, LogPriorSigma):
    def __call__(self, p):
        # last three parameters relate to sigma
        s = p[..., -3:]
        lnprior = LogPriorSigma.__call__(self, s)
        q = p[..., :-3]
        lnprior += LogPriorFixHa.__call__(self, q)
        return lnprior


class LogPriorLineSigma(LogPriorLine, LogPriorSigma):
    def __call__(self, p):
        # last three parameters relate to sigma
        s = p[..., -3:]
        lnprior = LogPriorSigma.__call__(self, s)
        q = p[..., :-3]
        lnprior += LogPriorLine.__call__(self, q)
        return lnprior


class LogPriorFlatSigma(LogPriorFlat, LogPriorSigma):
    def __call__(self, p):
        # last three parameters relate to sigma
        s = p[..., -3:]
        lnprior = LogPriorSigma.__call__(self, s)
        q = p[..., :-3]
        lnprior += LogPriorFlat.__call__(self, q)
        return lnprior


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
    lnpriors['line'] = LogPriorLine(ymin, ymax, zmin, zmax, fmax)
    lnpriors['flat'] = LogPriorFlat(ymin, ymax)

    lnpriors['fixha_sigma'] = LogPriorFixHaSigma(ymin, ymax, zmin, zmax, fmax)
    lnpriors['line_sigma'] = LogPriorLineSigma(ymin, ymax, zmin, zmax, fmax)
    lnpriors['flat_sigma'] = LogPriorFlatSigma(ymin, ymax)

    # ranges for mostly just for nestle
    # uniform and normal priors only used for testing

    ranges['fixha'] = np.array([[ymin, ymax], [zmin, zmax],
                                [fmin, fmax], [-1, 10], [-1, 10]])

    ranges['fixha_sigma'] = np.array([[ymin, ymax], [zmin, zmax],
                                      [fmin, fmax], [-1, 10], [-1, 10],
                                      [-1, 10], [-1, 10], [-0.1, 1.1]])

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


def plot_priors(x, y, outfile=None):
    xmin, xmax, ymin, ymax, zmin, zmax, fmin, fmax = get_ranges(x, y)
    fig, axarr = plt.subplots(2, 4, figsize=(15, 10))
    axarr = axarr.ravel()
    plot_prior(axarr[0], lnprior_continuum, [ymin, ymax],
               'continuum', args=[ymin, ymax])
    plot_prior(axarr[1], lnprior_redshift, [zmin - 0.01, zmax + 0.01], 'redshift',
               args=[zmin, zmax])
    plot_prior(axarr[2], lnprior_fluxNII, [-fmax, 5*fmax], 'fluxNII',
               args=[fmax])
    plot_prior(axarr[3], lnprior_NIIHa, [-1, 5], 'NIIHa')
    plot_prior(axarr[4], lnprior_absEWHa, [-1, 10], 'absEWHa')
    plot_prior(axarr[5], lnprior_sconst, [-1, 10], 'sconst')
    plot_prior(axarr[6], lnprior_sfactor, [-1, 10], 'sfactor')
    plot_prior(axarr[7], lnprior_badfrac, [-0.1, 1.1], 'badfrac')
    plt.tight_layout()
    if outfile is not None:
        plt.savefig(outfile)
