from __future__ import print_function, division
import numpy as np

from .pdfs import norm_pdf, _norm_pdf_C

wlHa = 6562.8
wlNIIa = 6548.1
wlNIIb = 6583.5
NIIratio = 3.06
linewidth = 7.5


def fixha_better_model(x, p):
    p = np.atleast_2d(p)
    continuum, redshift, fluxNII, NIIHa, absEWHa = p.T
    fluxHa = fluxNII / NIIHa - absEWHa * continuum
    zfactor = 1 + redshift
    model = continuum + norm_pdf(x, zfactor * wlHa, linewidth) * fluxHa
    model += (norm_pdf(x, zfactor * wlNIIa, linewidth) * fluxNII / NIIratio +
              norm_pdf(x, zfactor * wlNIIb, linewidth) * fluxNII)
    return model.squeeze()

par_fixha_better = ['continuum', 'redshift', 'flux NII', 'NII/Ha', 'abs EW Ha']


def fixha_model(x, p):
    p = np.atleast_2d(p)
    continuum, redshift, fluxHa, fluxNII = p.T
    zfactor = 1 + redshift
    model = continuum + norm_pdf(x, zfactor * wlHa, linewidth) * fluxHa
    model += (norm_pdf(x, zfactor * wlNIIa, linewidth) * fluxNII / NIIratio +
              norm_pdf(x, zfactor * wlNIIb, linewidth) * fluxNII)
    return model.squeeze()

par_fixha = ['continuum', 'redshift', 'flux Ha', 'flux NII']


def flat_model(x, p):
    p = np.atleast_2d(p)
    return x * 0 + p.T[0]

par_flat = ['continuum']


def line_model(x, p):
    p = np.atleast_2d(p)
    continuum, redshift, flux = p.T
    zfactor = 1 + redshift
    model = continuum + norm_pdf(x, zfactor * wlHa, linewidth) * flux
    return model.squeeze()

par_line = ['continuum', 'redshift', 'flux']


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
