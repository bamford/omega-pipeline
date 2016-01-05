from __future__ import print_function, division
import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from palettable.colorbrewer.qualitative import Set3_12 as palette
import corner
from emcee import autocorr
import nestle

from .convenience import flatten_without_burn, summary, log_evidence


# better-looking plots
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.figsize'] = (10, 5)
plt.rcParams['font.size'] = 10


def plot_chain(sampler, par, nburn=None, itemp=0, outfile=None):
    nwalkers = sampler.chain.shape[1]
    if nburn is None:
        nburn = sampler.chain.shape[2] // 10
    nrow = int(np.ceil(len(par) / 2.0))
    for i, p in enumerate(par):
        plt.subplot(nrow, 2, i + 1)
        for w in range(nwalkers):
            plt.plot(np.arange(len(sampler.chain[itemp, 0, :, 0])),
                     sampler.chain[itemp, w, :, i], 'r-', alpha=1.0 / nwalkers)
        plt.xlabel(p)
        aymin, aymax = plt.ylim()
        plt.vlines(nburn, aymin, aymax, linestyle=':')
        plt.ylim(aymin, aymax)
    if outfile is not None:
        plt.savefig(outfile)


def plot_hist(sampler, par, nburn=None, weights=None, outfile=None):
    if nburn is None:
        nburn = sampler.chain.shape[2] // 10
    nrow = int(np.ceil(len(par) / 2.0))
    for i, p in enumerate(par):
        plt.subplot(nrow, 2, i + 1)
        plt.hist(flatten_without_burn(sampler, nburn)[:, i], weights=weights,
                 bins=100, histtype='stepfilled', alpha=0.75)
        plt.xlabel(p)
    if outfile is not None:
        plt.savefig(outfile)


def plot_func(sampler, model, xmin, xmax, xdata, ydata, yerror,
              nburn=None, outfile=None, model_pars=[]):
    if nburn is None:
        nburn = sampler.chain.shape[2] // 10
    xchain = np.arange(xmin, xmax, 1.0)
    ychain = [model(xchain, p, *model_pars)
              for p in flatten_without_burn(sampler, nburn)]
    plt.errorbar(xdata, ydata, yerror, None, 'o')
    for i, y in enumerate(ychain[::100]):
        plt.plot(xchain, y, '-', alpha=0.1, c='green')
    plt.xlim(xmin, xmax)
    if outfile is not None:
        plt.savefig(outfile)


def plot_triangle(samples, par, model, xdata, ydata, yerror,
                  itemp=0, outfile=None, model_pars=[]):
    xchain = np.linspace(xdata.min(), xdata.max(), 100)
    ychain = [model(xchain, p, *model_pars) for p in samples]
    if len(par) == 1:
        plt.subplot(1, 2, 1)
        plt.hist(samples[:, 0], bins=100, histtype='stepfilled', alpha=0.75)
        plt.subplot(1, 2, 2)
    else:
        corner.corner(samples, labels=par)
        if len(par) % 2 == 0:
            n = 2
        else:
            n = int(np.ceil(len(par) / 2.0)) + 1
        plt.subplot(n, n, n)
    for y in ychain[::100]:
        plt.plot(xchain, y, 'r-', alpha=1000.0 / len(ychain))
    plt.errorbar(xdata, ydata, yerror, None, 'o')
    plt.subplots_adjust(wspace=0.2, hspace=0.2)
    if outfile is not None:
        plt.savefig(outfile)


def check_betas(sampler, nburn=500):
    logls = sampler.lnlikelihood[:, :, nburn:]
    logls = ma.masked_array(logls, mask=logls == -np.inf)
    mean_logls = logls.mean(axis=-1).mean(axis=-1)
    fig, ax = plt.subplots()
    plt.semilogx(sampler.betas, mean_logls, "-o")


def plot_autocorr(sampler, nburn, itemp=0, outfile=None):
    nwalkers = sampler.chain.shape[1]

    samples_before = sampler.chain[itemp, :, :nburn]
    samples_after = sampler.chain[itemp, :, nburn:]

    a_before = [autocorr.function(samples_before[i]) for i in range(nwalkers)]
    a_int_before = max(np.max([autocorr.integrated_time(samples_before[i])
                               for i in range(nwalkers)], 0))

    fig, [ax1, ax2] = plt.subplots(2)
    for a in a_before:
        ax1.plot(a[:200], "k", alpha=0.1)
    ax1.axhline(0, color="k")
    ax1.set_xlim(0, 200)
    ax1.set_xlabel(r"$\tau$")
    ax1.set_ylabel(r"Autocorrelation during burn-in")
    ax1.text(0.9, 0.9, '{}'.format(a_int_before), horizontalalignment='right',
             verticalalignment='top', transform=ax1.transAxes)

    a_after = [autocorr.function(samples_after[i]) for i in range(nwalkers)]
    a_int_after = max(np.max([autocorr.integrated_time(samples_after[i])
                              for i in range(nwalkers)], 0))

    for a in a_after:
        ax2.plot(a[:200], "k", alpha=0.1)
    ax2.axhline(0, color="k")
    ax2.set_xlim(0, 200)
    ax2.set_xlabel(r"$\tau$")
    ax2.set_ylabel(r"Autocorrelation after burn-in")
    ax2.text(0.9, 0.9, '{}'.format(a_int_after), horizontalalignment='right',
             verticalalignment='top', transform=ax2.transAxes)
    if outfile is not None:
        plt.savefig(outfile)


def plot_convergence(sampler, itemp=0, outfile=None):
    niter = sampler.chain.shape[2]
    iterno = np.arange(1, niter + 1)
    mean_over_walkers = sampler.chain[itemp].mean(axis=0)
    mean_vs_iteration = (np.cumsum(mean_over_walkers, axis=0).T / iterno).T
    stdev_over_walkers = sampler.chain[itemp].std(axis=0)
    stdev_vs_iteration = (np.cumsum(stdev_over_walkers, axis=0).T / iterno).T
    mean_track = ((mean_vs_iteration - mean_vs_iteration[-1]) /
                  stdev_vs_iteration[-1])
    stdev_track = (stdev_vs_iteration / stdev_vs_iteration[-1]) - 1
    fig, ax = plt.subplots()
    ax.set_color_cycle(palette.mpl_colors)
    ax.plot(iterno, mean_track, linestyle='-')
    ax.set_color_cycle(palette.mpl_colors)
    ax.plot(iterno, stdev_track, linestyle='--')
    if outfile is not None:
        plt.savefig(outfile)


def print_emcee(sampler, par, model, x, y, yerror, nburn, truths=None,
                outfile=None):
    mean, sigma = summary(flatten_without_burn(sampler, nburn),
                          par, truths=truths)
    samples = flatten_without_burn(sampler, nburn)
    logz, logzerr = log_evidence(sampler, nburn)
    ntemp = sampler.chain.shape[0]
    print()
    if outfile is not None:
        pdf = PdfPages(outfile)
        def page(title):
            plt.suptitle(title)
            pdf.savefig()
            plt.close()
    else:
        def page(title):
            plt.suptitle(title)
    plot_triangle(samples, par, model, x, y, yerror)
    page('triangle')
    check_betas(sampler, nburn)
    page('check_betas')
    for itemp in (0, ntemp // 2, ntemp - 1):
        plot_chain(sampler, par, nburn, itemp=itemp)
        page('chain itemp={}'.format(itemp))
        plot_autocorr(sampler, nburn, itemp=itemp)
        page('autocorr itemp={}'.format(itemp))
        plot_convergence(sampler, itemp=itemp)
        page('convergence itemp={}'.format(itemp))
    pdf.close()


def print_nestle(res, par, model, x, y, yerror, outfile=None):
    print(res.summary)
    print(nestle.mean_and_cov(res.samples, res.weights))
    plot_triangle(res.samples, par, model, x, y, yerror, weights=res.weights)
    if outfile is not None:
        plt.savefig(outfile)
