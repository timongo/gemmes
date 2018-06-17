#!/usr/bin/python

# Common
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

_fs = (18,8)
_dmargin = dict(left=0.06, right=0.98, bottom=0.04, top=0.96,
                wspace=0.2, hspace=0.2)

def plot_basic(lGemInt,
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
    params = lGemInt[0].params
    nparams = lGemInt[0].sol['U'].shape[1]

    # Prepare figure and axes grid
    fig = plt.figure(fignumber,figsize=fs)
    axarr = GridSpec(3,5, **dmargin)
    # Create and format axes
    lax = []
    daxprop = {'frameon':True}#, 'facecolor':'w'}
    for ii in range(0,nparams):
        if ii==0:
            lax.append(fig.add_subplot(axarr[ii], **daxprop))
        else:
            # Make sure a zoom on an axes is replicated in all (for t)
            lax.append(fig.add_subplot(axarr[ii], sharex=lax[0], **daxprop))
        tit = params[ii]['meaning']
        lab = "%s (%s)"%(params[ii]['symbol'], params[ii]['units'])
        lax[ii].set_title(tit)
        lax[ii].set_ylabel(lab)
        if ii<10:
            lax[ii].set_xticklabels([])
        else:
            lax[ii].set_xlabel("t (years)") 

    # Add data to the axes
    for gg in lGemInt:
        for ii in range(0,nparams):
            lax[ii].plot(gg.sol['t'], gg.sol['U'][:,ii],
                         ls='-', lw=1., label=gg.name)

    if draw:
        fig.canvas.draw()
    return lax

