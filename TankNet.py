# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:12:29 2022

Resonant Converters

The tank networks can be written in the form:
                ——————
    ———————————| jXs |—————————————————————
               ——————      |              |
                           |              |
                         ——————         —————
                        | jXp |        | Re |
                        ——————         —————
                           |             |
                           |             |
    ——————————————————————————————————————

Zi0 = jXs for Re = 0;
Zi∞ = j(Xs + Xp) for Re = ∞;

The unloaded tank transfer function is:
                Xp
    H∞(ω) = —————————
             Xs + Xp

The matched-load impedance is:
               jXsXp
    Zo0(ω) = ————————— = jXs * H∞(ω)
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

def Tank_Network(Xs, Xp, Re, f0):
    
    # Rcrit = abs(Xp) * math.sqrt(-Xs / (Xs + Xp))
    # print(Rcrit)
    
    # Transfer function.
    Gtank = resf.Parallel_Impedance(Xp, Re) / (resf.Parallel_Impedance(Xp, Re) + Xs)
    
    print(Gtank)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
    gain, phase, omega = ctrl.bode(Gtank, ctrl.logspace(4, 7), plot = False)
    
    # Drawing based on normalized frequency.
    fig, ax = plt.subplots(2, 1, figsize = (11, 8))
    
    ax[0].grid(which="both", ls=':')
    ax[1].grid(which="both", ls=':')
    ax[1].set_xlabel('f$_n$, Normalized Frequency')
    ax[0].set_ylabel('M$_g$, Voltage Gain [dB]')
    ax[1].set_ylabel('Phase [deg]')
    ax[0].set_xlim(0.1, 10)
    ax[1].set_xlim(0.1, 10)
    ax[1].set_yticks([-450, -405, -360, -315, -270, -225])
    
    ax[0].semilogx(omega/(2*np.pi*f0), 20*np.log10(gain))
    ax[1].semilogx(omega/(2*np.pi*f0), phase*180/np.pi)
    
    fig.tight_layout()
    
    
def LLC_Tank_Network(Vi, Vo, Po, Cr, Lr, Lm, n):
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
    
    # Standarized tank resonant convertor.
    s = ctrl.tf('s')
    Xs = s*Lr + 1 / (s*Cr)
    Xp = s*Lm
    
    # Draw tank transfer fucntion and bode plot.
    Tank_Network(Xs, Xp, Re, f0)
    
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
    
    Qe = [0.1, 0.2, 0.5, 0.8, 1, 2, 5, 8, 10]
    
    for i, value in enumerate(Qe):
        Mg = abs(Ln * fn**2 / ((Ln + 1) * fn**2 - 1 + ((fn**2 - 1) * fn * value * Ln)*1j))
        ax.semilogx(fn, Mg, label = 'Qe = ' + str(value))
        ax.legend(loc=1, fontsize=16, shadow=True)

    
# Example.
# LLC_Tank_Network(27.3e-9, 60e-6, 210e-6, 90)

# LLC Tank parameters:
Vi = 426            # Bulk voltage = 426V
Vo = 12.5           # Output voltage = 12.5V
Po = 4200           # Full load power = 4200W
Cr = 220e-9         # 220nF
Lr = 9.7e-6         # 9.7uH
Lm = 145e-6         # 145uH
n = 17              # Transformer turns ratio.

LLC_Tank_Network(Vi, Vo, Po, Cr, Lr, Lm, n)

