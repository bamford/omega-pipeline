import numpy as np
from scipy import stats
from .pdfs import _norm_pdf_C
from .omega_models import wlHa, linewidth, get_ranges


class InitPars:
    def __init__(self, ranges):
        self.ranges = np.asarray(ranges)
        self.npar = len(self.ranges)
        self.rv = stats.uniform(self.ranges[:, 0],
                                self.ranges[:, 1] - self.ranges[:, 0])

    def __call__(self, shape):
        shape = list(shape) + [self.npar]
        return self.rv.rvs(size=shape)


def create_init_pars(x, y):
    xmin, xmax, ymin, ymax, zmin, zmax, fmin, fmax = get_ranges(x, y)

    init_pars = {}
    init_pars['fixha'] = InitPars([[ymin, ymax],
                                   [zmin, zmax],
                                   [0, fmax],
                                   [0, 3],
                                   [0, 5]])

    init_pars['fixha_poor'] = InitPars([[ymin, ymax],
                                   [zmin, zmax],
                                   [0, fmax],
                                   [0, fmax]])

    init_pars['line'] = InitPars([[ymin, ymax],
                                  [zmin, zmax],
                                  [0, fmax]])

    init_pars['flat'] = InitPars([[ymin, ymax]])

    for key in init_pars.keys():
        ranges = init_pars[key].ranges
        ranges = np.concatenate((ranges, [[0, 1],
                                          [0.5, 2.0],
                                          [0, 0.1]]))
        init_pars[key + '_sigma'] = InitPars(ranges)

    return init_pars
