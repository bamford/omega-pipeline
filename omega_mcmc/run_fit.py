from .convenience import LogLikelihood, run_emcee
from .plots import print_emcee
from .omega_models import *
from .omega_priors import *


def run_fixha_poor(x, y, yerror, lnpriors, init_pars,
                   run_emcee_conf, truths=None,
                   label='sim', output=True):
    sampler = run_emcee(
        LogLikelihood(fixha_poor_model,
                      x, y, yerror),
        lnpriors['fixha_poor_normal'],
        init_pars['fixha_poor'],
        outfilename='{}_sampler_fixha_poor'.format(label),
        **run_emcee_conf)
    if output:
        stats = print_emcee(sampler, par_fixha_poor,
                            fixha_poor_model, x, y, yerror,
                            outfilename='{}_fixha_poor'.format(label),
                            **run_emcee_conf)
        return sampler, stats


def run_fixha(x, y, yerror, lnpriors, init_pars,
              run_emcee_conf, truths=None,
              label='sim', output=True):
    sampler = run_emcee(
        LogLikelihood(fixha_model,
                      x, y, yerror),
        lnpriors['fixha'],
        init_pars['fixha'],
        outfilename='{}_sampler_fixha'.format(label),
        **run_emcee_conf)
    if output:
        stats = print_emcee(sampler, par_fixha, fixha_model,
                            x, y, yerror, outfilename='{}_fixha'.format(label),
                            **run_emcee_conf)
        return sampler, stats


def run_line(x, y, yerror, lnpriors, init_pars,
             run_emcee_conf, truths=None,
             label='sim', output=True):
    sampler = run_emcee(
        LogLikelihood(line_model, x, y, yerror),
        lnpriors['line_normal'],
        init_pars['line'],
        outfilename='{}_sampler_line'.format(label),
        **run_emcee_conf)
    if output:
        stats = print_emcee(sampler, par_line, line_model, x, y, yerror,
                            outfilename='{}_line'.format(label),
                            **run_emcee_conf)
        return sampler, stats


def run_flat(x, y, yerror, lnpriors, init_pars,
             run_emcee_conf, truths=None,
             label='sim', output=True):
    sampler = run_emcee(
        LogLikelihood(flat_model, x, y, yerror),
        lnpriors['flat_normal'],
        init_pars['flat'],
        outfilename='{}_sampler_flat'.format(label),
        **run_emcee_conf)
    if output:
        stats = print_emcee(sampler, par_flat, flat_model, x, y, yerror,
                            outfilename='{}_flat'.format(label),
                            **run_emcee_conf)
        return sampler, stats
