import numpy as np
from scipy import stats
from .pdfs import _norm_pdf_C
from .omega_models import wlHa, linewidth, get_ranges


class InitPars:
    def __init__(self, ranges):
        ranges = np.asarray(ranges)
        self.npar = len(ranges)
        self.rv = stats.uniform(ranges[:, 0], ranges[:, 1] - ranges[:, 0])

    def __call__(self, shape):
        shape =  list(shape) + [self.npar]
        return self.rv.rvs(size=shape)


def create_init_pars(x, y):
    xmin, xmax, ymin, ymax, zmin, zmax, fmin, fmax = get_ranges(x, y)

    init_pars = {}
    init_pars['fixha_better'] = InitPars([[ymin, ymax],
                                          [zmin, zmax],
                                          [0, fmax],
                                          [0, 3],
                                          [0, 5]])

    init_pars['fixha'] = InitPars([[ymin, ymax],
                                   [zmin, zmax],
                                   [0, fmax],
                                   [0, fmax]])

    init_pars['line'] = InitPars([[ymin, ymax],
                                  [zmin, zmax],
                                  [0, fmax]])

    init_pars['flat'] = InitPars([[ymin, ymax]])

    return init_pars
