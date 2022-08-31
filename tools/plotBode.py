# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 14:56:43 2022

@author: Vern Chen
"""

import numpy as np
import control.matlab as ctrl
import matplotlib.pyplot as plt

# Define line style for draw multi-line in the same plot.
def linestyle_generator():
    linestyle = ['-', '--', '-.', ':']
    lineID = 0
    while True:
        yield linestyle[lineID]
        lineID = (lineID + 1) % len(linestyle)

# Draw bode plot based on frequency fs instead of omega (2*pi*fs)
def bodeplot_set(fig_ax, *args):
    fig_ax[0].grid(which="both", ls='-')
    fig_ax[0].set_ylabel('Magnitude [dB]')

    fig_ax[1].grid(which="both", ls='-')
    fig_ax[1].set_xlabel('Frequency [Hz]')
    fig_ax[1].set_ylabel('Phase [deg]')
    
    if len(args) > 0:
        fig_ax[0].legend(loc=args[0], fontsize=16, ncol=2, shadow=True)
    if len(args) > 1:
        fig_ax[1].legend(loc=args[1], fontsize=16, shadow=True)

# Draw bode plot according to transfer function.
def Tf2Bodeplot(Tf, Gstring):
    
    LS = linestyle_generator()
    
    fig, ax = plt.subplots(2, 1, figsize = (11, 8))
    
    for i in range(len(Tf)):
        # Get the Gain and Phase and angular frequency 1Hz ~ 10MHz.
        gain, phase, omega = ctrl.bode(Tf[i], ctrl.logspace(1, 8), plot = False)
        
        pltargs = {'ls': next(LS), 'label': (Gstring[i])}
        ax[0].semilogx(omega/(2*np.pi), 20*np.log10(gain), **pltargs)
        ax[1].semilogx(omega/(2*np.pi), phase*180/np.pi, **pltargs)
        
    bodeplot_set(ax, 3, 3)
    ax[1].set_ylim(-200,10)
    ax[1].set_yticks([-180, -135, -90, -45, 0])
    fig.tight_layout()
