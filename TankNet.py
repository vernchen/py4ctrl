# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:12:29 2022

Resonant Converters

The tank networks can be written in the form:
                —————
    ———————————| Xs |—————————————————————
               —————       |             |
                           |             |
                         —————         —————
                        | Xp |        | Re |
                        —————         —————
                           |            |
                           |            |
    —————————————————————————————————————

Zi0 = Xs for Re = 0;
Zi∞ = (Xs + Xp) for Re = ∞;

The unloaded tank transfer function is:
                Xp
    H∞(ω) = —————————
             Xs + Xp

The matched-load impedance is:
               XsXp
    Zo0(ω) = ————————— = Xs * H∞(ω)
              Xs + Xp

Matched-load resistance occurs at Re = Ro0 = ||Zo0||.

The critical load resistance at the boundary between ZVS and ZCS is:
                        ————————            —————————————
                       /   Zi∞             /      Xs
    Rcrit = ||Zo0|| * / - —————  = |Xp| * / - —————————
                     V     Zi0           V     Xs + Xp
                     
The frequency f = fm, where ||Zi∞|| = ||Zi0||, can be shown to occur at the frequency
where Xs = -Xp / 2.

If we define:
    conversion ratio M = Vout / Vin
    Mormalized load current J = Iout * Ro / Vin
    Effective quality factor as Qe = R0 / Re
Then the elliptical output characteristic can be written:
    
      M  2      J  2
    (———)  +  (———)  =  1
      a         b
      
and the control characteristic can be written:
    
                  1
    M = ————————————————————
           —————————————————
          / 1             2
         / ———— + (Qe / b)
        V  a^2
        
Where the parameters a and b are given by:
    
                        |Xp|
    a = ||H∞(ω)|| = ———————————
                     |Xs + Xp|
    
        ||H∞(ω)|| * R0    R0
    b = —————————————— = ————
          ||Zo0(ω)||     |Xs|
    

Vern's Noteboot Power Electronics --> Resonant Converters

@author: Vern Chen

"""

import math as math
import numpy as np
import control.matlab as ctrl
import matplotlib.pyplot as plt
import ResurrectionF as resf

def Tank_Network_TF(Xs, Xp, Re, f0):
    
    # Transfer function of major load condition, Re as full load.
    Gt100 = resf.Parallel_Impedance(Xp, Re) / (resf.Parallel_Impedance(Xp, Re) + Xs)
    R = Re * 2      # 50% load, half load.
    Gt50 = resf.Parallel_Impedance(Xp, R) / (resf.Parallel_Impedance(Xp, R) + Xs)
    R = Re * 5      # 20% load.
    Gt20 = resf.Parallel_Impedance(Xp, R) / (resf.Parallel_Impedance(Xp, R) + Xs)
    R = Re * 10     # 10% load.
    Gt10 = resf.Parallel_Impedance(Xp, R) / (resf.Parallel_Impedance(Xp, R) + Xs)
    print(Gt100)
    
    # Drawing based on normalized frequency.
    fig, ax = plt.subplots(2, 1, figsize = (11, 8))
    
    ax[0].grid(which="both", ls=':')
    ax[1].grid(which="both", ls=':')
    ax[1].set_xlabel('f$_n$, Normalized Frequency')
    ax[0].set_ylabel('M$_g$, Voltage Gain [dB]')
    ax[1].set_ylabel('Phase [deg]')
    ax[0].set_xlim(0.1, 10)
    ax[1].set_xlim(0.1, 10)
    ax[0].set_ylim(-20, 20)
    ax[1].set_yticks([-90, -45, 0, 45, 90, 135])
    
    Gtank = (Gt10, Gt20, Gt50, Gt100)
    Gstring = ('10% load', '20% load', '50% load', '100% load')
        
    for i, value in enumerate(Gtank):
        # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
        gain, phase, omega = ctrl.bode(value, ctrl.logspace(4, 7), plot = False)
        
        ax[0].semilogx(omega/(2*np.pi*f0), 20*np.log10(gain), label = Gstring[i])
        ax[1].semilogx(omega/(2*np.pi*f0), phase*180/np.pi + 360, label = Gstring[i])
        
    ax[0].legend(loc=1, fontsize=16, shadow=True)
    ax[1].legend(loc=1, fontsize=16, shadow=True)
    fig.tight_layout()

def LLC_Tank_Critical_load(Cr, Lr, Lm, f0, fp, Re):
    
    fs = np.linspace(fp+1, f0-1, 1000)
    W = 2*np.pi*fs
    
    # Rf = np.linspace(Re, Re, 1000)
    '''
    Input impedance with load Re is:
        Zi = (jωLm)||Re + jωLr + 1/jωCr
        
    R critical occurs when imaginary part of Zi, Im(Zi) = 0
    Solve equation Im(Zi) = 0, and get Rcrit from fp->f0.
    '''
    Rcrit = np.sqrt((W**2)*(Lm**2)*(1/(W*Cr) - W*Lr) / (W*Lm + W*Lr - 1/(W*Cr)))
    
    fig, ax = plt.subplots(1, 1, figsize = (11, 8))
    
    ax.grid(which="both", ls=':')
    ax.set_xlabel('f$_s$, Switching Frequency [Hz]')
    ax.set_ylabel('||R|| , Load impedance [dB]')
    
    ax.semilogx(fs, 20*np.log10(Rcrit), label = 'R$_{crit}$')
    # ax.semilogx(fs, 20*np.log10(Rf), label = 'R$_{full}$', linestyle = '-.')
    
    # Get the frequency at full load ZVS/ZCS boundary.
    y = 20*np.log10(Re)
    for i, value in enumerate(Rcrit):
        if value < Re:
            x = fs[i]
            break
    
    ax.axhline(y, color='darkorange', lw=1, ls = '--')
    ax.axvline(x, color='darkorange', lw=1, ls = '--')
    ax.scatter(x, y, color='k', lw=2)
    text = '  fs   = ' + str(round(x/1000, 2)) + \
        ' kHz\n ||R|| = ' + str(round(y, 2)) + ' dB\n'
    ax.annotate(text, xy=(x, y), xytext=(x-2e4, y-1e1), arrowprops = dict(facecolor ='g'))
    ax.text(6e4, 5e1, 'ZVS', color='g', size=39)
    ax.text(4e4, 0e1, 'ZCS', color='salmon', size=39)
    ax.legend(loc=3, fontsize=16, shadow=True)
    
def LLC_Tank_Input_Impedance(Cr, Lr, Lm):
    # Switching frequency range from 1e2 to 1e7.
    fs = np.linspace(1e2, 1e7, 10000)
    
    # Convert to angular frequency.
    W = 2*np.pi*fs
    
    Xs = 1j*W*Lr + 1 / (1j*W*Cr)
    Xp = 1j*W*Lm
    
    # Zis represent input impedance with output short circuit.
    Zis = abs(Xs)
    # Zio represent input impedance with output open circuit.
    Zio = abs(Xs + Xp)
    
    '''
    # fm defined at frequency where Zis = Zio:
        abs(jwLr + 1 / jwCr) = abs(jw(Lr + Lm) + 1 / jwCr)
        abs(jwLr - j / wCr) = abs(jw(Lr + Lm) - j / wCr)
        abs(wLr - 1 / wCr) = abs(w(Lr + Lm) - 1 / wCr)
        1 / wCr - wLr = w(Lr + Lm) - 1 / wCr
        2 / wCr = w(2Lr + Lm)
    '''
    fm = math.sqrt(2 / (Cr*(2*Lr + Lm))) / (2*math.pi)
    Zm = abs(1j*2*math.pi*fm*Lr + 1 / (1j*2*math.pi*fm*Cr))
    print('fm = ', fm, '\nZm = ', Zm)
    
    # Drawing Tank input impedance.
    fig, ax = plt.subplots(1, 1, figsize = (11, 8))

    ax.grid(which="both", ls=':')
    ax.set_xlabel('f$_s$, Switching Frequency [Hz]')
    ax.set_ylabel('||Z$_i$|| , Input Impedance [dB]')
    
    ax.semilogx(fs, 20*np.log10(Zis), label = 'Z$_{is}$')
    ax.semilogx(fs, 20*np.log10(Zio), label = 'Z$_{io}$')
    
    x = fm
    y = 20*np.log10(Zm)
    
    ax.scatter(x, y, color='k', lw=2)
    text = '  fs   = ' + str(round(x/1000, 2)) + \
        ' kHz\n ||Zi|| = ' + str(round(y, 2)) + ' dB\n'
    ax.annotate(text, xy=(x, y), xytext=(x+5e4, y+3e1), arrowprops = \
                dict(arrowstyle = "->", connectionstyle = \
                     "angle, angleA = 0, angleB = 90, rad = 10"))
    
    ax.legend(loc=3, fontsize=16, shadow=True)
    
    
def LLC_Tank_Network(Vi, Vo, Po, Cr, Lr, Lm, n):
    
    LLC_Tank_Input_Impedance(Cr, Lr, Lm)
    
    # f0 is peak resonant frequency at load short circuit.
    # fp is peak resonant frequency at load open circuit (no load).
    f0 = 1 / (2*math.pi*math.sqrt(Lr*Cr))
    fp = 1 / (2*math.pi*math.sqrt((Lr+Lm)*Cr))
    
    # Inductance ratio, Normalized inductance.
    Ln = Lm / Lr
    
    '''
    The AC equivalent load resistance:
              Voe     8 * n^2     Vo     8 * n^2
        Re = ————— = ————————— * ———— = ————————— * Rl
              Ioe      pi^2       Io      pi^2 
    '''
    Re = (8 * n**2/math.pi**2) * Vo**2 / Po
    
    LLC_Tank_Critical_load(Cr, Lr, Lm, f0, fp, Re)
    
    # Standarized tank resonant convertor.
    s = ctrl.tf('s')
    Xs = s*Lr + 1 / (s*Cr)
    Xp = s*Lm
    
    # Draw tank transfer fucntion and bode plot.
    Tank_Network_TF(Xs, Xp, Re, f0)
    
    '''
    The quality factor of the series resonant circuit is defined as:
               _________
              V Lr / Cr
        Qe = ————————————
                  Re

    Increasing Qe is due to an increased load.
    Short circuit, Re -> 0, Qe -> ∞.
    Open circuit, Re -> ∞, Qe -> 0.
    '''
    Qe = math.sqrt(Lr/Cr) / Re
    
    fs = np.linspace(0.1*f0, 10*f0, 1000)
    
    # Normalized frequency, series resonant frequency (f0) can be selected as 
    # the base for normalization.
    fn = fs / f0
    
    # Mg_∞ is a special value of Mg that presents the gain value at no load 
    # when fn is approaching infinity. In other words, at no load, the gain 
    # curve (Mg) is approaching an asymptotic horizontal line.
    Mg_inf = Ln / (Ln + 1)
    
    '''
    Voltage gain Function, normalized from:
              Voe    |       (jωLm) || Re          |
        Mg = ————— = | ——————————————————————————— |
              Vge    | (jωLm)||Re + jωLr + 1/jωCr  |
    '''
    Mg = abs(Ln * fn**2 / ((Ln + 1) * fn**2 - 1 + ((fn**2 - 1) * fn * Qe * Ln)*1j))
    
    print('f0 = ', f0, '\nfp = ', fp, '\nLn = ', Ln, '\nRe = ', Re, '\nQe = ', Qe)
    print('Mg_inf = ', Mg_inf)
    
    fig, ax = plt.subplots(1, 1, figsize = (11, 8))

    ax.grid(which="both", ls=':')
    ax.set_xlabel('f$_n$, Normalized Frequency')
    ax.set_ylabel('M$_g$, Voltage Gain')
    # ax.set_xlim(0.8, 1.2)
    ax.set_ylim(0.9, 1.1)
    
    Qe = [0, 0.1, 0.2, 0.38, 0.5, 0.8, 1, 2, 5, 8, 10]
    
    for i, value in enumerate(Qe):
        Mg = abs(Ln * fn**2 / ((Ln + 1) * fn**2 - 1 + ((fn**2 - 1) * fn * value * Ln)*1j))
        ax.semilogx(fn, Mg, label = 'Qe = ' + str(value))
        ax.legend(loc=1, fontsize=16, shadow=True)
        
        
# Example, LLC Tank parameters:
Vi = 426            # Bulk voltage = 426V
Vo = 12.5           # Output voltage = 12.5V
Po = 4200/2         # Full load power = 4200W, divided by 2 as interleaved.
Cr = 220e-9         # 220nF
Lr = 9.7e-6         # 9.7uH
Lm = 145e-6         # 145uH
n = 17              # Transformer turns ratio.

LLC_Tank_Network(Vi, Vo, Po, Cr, Lr, Lm, n)

