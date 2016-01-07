from __future__ import print_function, division
import numpy as np

from .pdfs import norm_logpdf, norm_pdf, _norm_pdf_C


wlHa = 6562.8
wlNIIa = 6548.1
wlNIIb = 6583.5
NIIratio = 3.06
linewidth = 7.5


def fixha_model(x, p):
    p = np.atleast_2d(p)
    continuum, redshift, fluxNII, NIIHa, absEWHa = p.T[:5]
    fluxHa = fluxNII / NIIHa - absEWHa * continuum
    zfactor = 1 + redshift
    model = continuum + norm_pdf(x, zfactor * wlHa, linewidth) * fluxHa
    model += (norm_pdf(x, zfactor * wlNIIa, linewidth) * fluxNII / NIIratio +
              norm_pdf(x, zfactor * wlNIIb, linewidth) * fluxNII)
    return model.squeeze()

par_fixha = ['continuum', 'redshift', 'flux NII', 'NII/Ha', 'abs EW Ha']


def fixha_poor_model(x, p):
    p = np.atleast_2d(p)
    continuum, redshift, fluxHa, fluxNII = p.T[:4]
    zfactor = 1 + redshift
    model = continuum + norm_pdf(x, zfactor * wlHa, linewidth) * fluxHa
    model += (norm_pdf(x, zfactor * wlNIIa, linewidth) * fluxNII / NIIratio +
              norm_pdf(x, zfactor * wlNIIb, linewidth) * fluxNII)
    return model.squeeze()

par_fixha_poor = ['continuum', 'redshift', 'flux Ha', 'flux NII']


def flat_model(x, p):
    p = np.atleast_2d(p)
    return x * 0 + p.T[0]

par_flat = ['continuum']


def line_model(x, p):
    p = np.atleast_2d(p)
    continuum, redshift, flux = p.T[:3]
    zfactor = 1 + redshift
    model = continuum + norm_pdf(x, zfactor * wlHa, linewidth) * flux
    return model.squeeze()

par_line = ['continuum', 'redshift', 'flux']


models = dict(fixha=fixha_model,
              fixha_poor=fixha_poor_model,
              line=line_model,
              flat=flat_model)

pars = dict(fixha=par_fixha,
            fixha_poor=par_fixha_poor,
            line=par_line,
            flat=par_flat)


def get_ranges(x, y):
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    yrange = ymax - ymin
    ymin -= 2 * yrange
    ymax += 2 * yrange
    zmin = xmin / wlHa - 1
    zmax = xmax / wlHa - 1
    zcutoff = 3 * linewidth / wlHa
    zmin -= zcutoff
    zmax += zcutoff
    fmax = 10 * (ymax - ymin) * _norm_pdf_C * linewidth
    fmin = -0.1 * fmax
    return xmin, xmax, ymin, ymax, zmin, zmax, fmin, fmax


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


class LogLikelihoodSigma:
    def __init__(self, model, x, y, sigma):
        self.model = model
        self.y = y
        self.x = x
        self.sigma = sigma

    def __call__(self, p):
        # last three parameters relate to sigma
        # a constant error (independent of that given in the data)
        # a uniform scaling of the given errors
        # a fraction of bad data points, with much larger errors
        q = p[..., :-3]
        yM = np.atleast_2d(self.model(self.x, q))
        s = p[..., -3:]
        sconst, sfactor, badfrac = s.T
        sigma = np.sqrt(sconst**2 + (sfactor * self.sigma)**2)
        badsigma = 100 * sigma
        lnL_good = norm_logpdf(self.y, yM, sigma)
        lnL_good += np.log(1 - badfrac)
        lnL_bad = norm_logpdf(self.y, yM, badsigma)
        lnL_bad += np.log(badfrac)
        # using np.logaddexp helps maintain numerical precision
        lnL = np.logaddexp(lnL_good, lnL_bad).sum(axis=1)
        return lnL.squeeze()


par_sigma = ['sconst', 'sfactor', 'badfrac']

for key in pars.keys():
    pars[key + '_sigma'] = pars[key] + par_sigma
