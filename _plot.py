#!/usr/bin/python

# Common
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import itertools as it
import os
import sys

_fs = (18,8)
_dmargin = dict(left=0.06, right=0.96, bottom=0.06, top=0.96,
                wspace=0.2, hspace=0.2)

eco_vars = ['omega','lambda','d','N','pi']
gdp_vars = ['Y','Y0','g','i']
emi_vars = ['Eind','Eland','sigma','gsigma','n']
clim_vars = ['T','T0','CO2at','CO2up','CO2lo']
pcar_vars = ['pC','pbs']
dam_vars = ['A','DY','DK','D']

def plot_basic(lGemInt,
               vars = eco_vars,
               fs=None, dmargin=None, draw=True, fignumber=1):
    """ Very basic plotting function to visualize the time traces of all parameters
    
    fs = figure size in inches
    draw = bool, useful if heavy figure or embedded in GUI that require not
           drawing immediately
    """

    # Default inputs
    if fs is None:
        fs = _fs
    if dmargin is None:
        dmargin = _dmargin

    # Prepare data
    if type(lGemInt) is not list:
        lGemInt = [lGemInt]
    nvars = len(vars)

    params = lGemInt[0].params.lparams

    # Prepare figure and axes grid
    fig = plt.figure(fignumber,figsize=fs)
    plt.clf()
    
    nrows = int(np.sqrt(nvars))
    ncols = int(np.ceil(nvars/nrows))
    axarr = fig.add_gridspec(nrows,ncols, **dmargin)
    # Create and format axes
    lax = []
    daxprop = {'frameon':True}#, 'facecolor':'w'}
    for i in range(0,nvars):
        ii = searchdict(params,vars[i])
        if i==0:
            lax.append(fig.add_subplot(axarr[i], **daxprop))
        else:
            # Make sure a zoom on an axes is replicated in all (for t)
            # lax.append(fig.add_subplot(axarr[i], sharex=lax[0], **daxprop))
            lax.append(fig.add_subplot(axarr[i], **daxprop))
        tit = params[ii]['meaning']
        lab = "%s (%s)"%(params[ii]['symbol'], params[ii]['units'])
        lax[i].set_title(tit)
        lax[i].set_ylabel(lab)
        string = params[ii]['ylim']
        if string!='':
            lax[i].set_ylim(np.array([string.split(',')[i] for i in range(2)]).astype(float))
        if i<(nrows-1)*ncols:
            lax[i].set_xticklabels([])
        else:
            lax[i].set_xlabel("t (years)") 

    # Add data to the axes
    for gg in lGemInt:
        for i in range(0,nvars):
            ii = searchdict(params,vars[i])            
            t = gg.sol['t']
            u = gg.sol['U'][:,ii]
            ti,tf = t[0],t[-1]
            if params[ii]['plot']=='log':
                lax[i].semilogy(t,u,ls='-', lw=2., label=gg.name)
            else:
                lax[i].plot(t,u,ls='-', lw=2., label=gg.name)
            lax[i].grid(True)
            lax[i].set_xlim([ti,tf])

    if draw:
        fig.canvas.draw()
    return lax

def searchdict(ldict,varname):
    for i in range(len(ldict)):
        if ldict[i]['var']==varname:
            return i

