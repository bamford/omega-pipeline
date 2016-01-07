from .convenience import run_emcee
from .plots import print_emcee
from .omega_models import models, pars, LogLikelihood, LogLikelihoodSigma


def run(modelname, sigma_model,
        x, y, yerror,
        lnpriors, init_pars,
        run_emcee_conf, truths=None,
        label='sim', output=True):
    if sigma_model:
        LogLikeClass = LogLikelihoodSigma
        suffix = '_sigma'
    else:
        LogLikeClass = LogLikelihood
        suffix = ''
    sampler = run_emcee(
        LogLikeClass(models[modelname], x, y, yerror),
        lnpriors[modelname + suffix],
        init_pars[modelname + suffix],
        outfilename='{}_sampler_{}'.format(label, modelname),
        **run_emcee_conf)
    if output:
        stats = print_emcee(
            sampler, pars[modelname + suffix],
            models[modelname], x, y, yerror,
            outfilename='{}_{}'.format(label, modelname + suffix),
            **run_emcee_conf)
        return sampler, stats


def run_fixha_poor(*args, **kwargs):
    return run('fixha_poor', False, *args, **kwargs)


def run_fixha(*args, **kwargs):
    return run('fixha', False, *args, **kwargs)


def run_line(*args, **kwargs):
    return run('line', False, *args, **kwargs)


def run_flat(*args, **kwargs):
    return run('flat', False, *args, **kwargs)


def run_fixha_sigma(*args, **kwargs):
    return run('fixha', True, *args, **kwargs)


def run_line_sigma(*args, **kwargs):
    return run('line', True, *args, **kwargs)


def run_flat_sigma(*args, **kwargs):
    return run('flat', True, *args, **kwargs)
