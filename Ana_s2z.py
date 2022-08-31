# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 17:27:56 2022

@author: Vern Chen

"""

import numpy as np
import matplotlib.pyplot as plt

'''
This function shows that the magnitude response of the discrete-time integrator 
approximates very well with magnitude response of the continuous-time integrator 
at frequencies sufficiently low compared to the sampling frequency.

To evaluate frequency response, let  s = jw

The s-domain continuous-time integral compensator transfer function
              w0
Gc(jw) = -j * ——
              w
              
The z-domain discrete-time integral compensator transfer function
               w0 * Ts   cos(w*Ts/2)
Gcd(jw) = -j * ——————— * ———————————
                  2      sin(w*Ts/2)

The phase responses of G_cd  and G_c  are exactly the same at all frequencies. 
Both transfer functions exhibit −90◦ phase at all frequencies.

Vern's OneNote Signals and Systems --> Introduction to Discrete-Time Systems
'''
def Magnitude_domain_z_VS_s():
    
    fs = 100 * 10**3                  # Hz, Switching frequency fs of convertor.
    Ts = 1 / fs                       # Second, Switching period.
    f = np.linspace(1, 10*fs, 1000)
    
    omega = 2*np.pi*f
    omega_0 = 2*np.pi*100
    
    Gain_Gcd = abs((np.cos(omega*Ts/2)/np.sin(omega*Ts/2)) * omega_0 * Ts / 2)
    Gain_Gc = omega_0/omega
    
    fig, ax = plt.subplots(figsize = (11, 8))
    ax.grid(which="both", ls=':')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Magnitude [dB]')
    ax.set_ylim(-100, 60)
    
    ax.semilogx(omega/(2*np.pi), 20*np.log10(Gain_Gcd), label = 'G$_{cd}$')
    ax.semilogx(omega/(2*np.pi), 20*np.log10(Gain_Gc), label = 'G$_c$')
    ax.legend(loc=3, fontsize=16, shadow=True)
    fig.tight_layout()
    
'''
The continuous-time to discrete-time mapping approaches.

As presented by its name, z-domain transfer and s-domain transfer:
    z^-1 repersent delay in time domain;
    s^-1 repersent integral in time domain;

With unit delay deduction, the frequency response of a z-domain transfer function 
can be found by replacing z^(−1)  with e^(−sTs), and then s with jω

    z = e^(jwTs) = cos(wTs) + j*sin(wTs)

With integral deduction by bilinear mapping of Continuous Time to Discrete Time
         2    z - 1
    s = —— * ———————
        Ts    z + 1
        
To solving for z in terms of s:
         1 + jwTs/2
    z = —————————————
         1 - jwTs/2
    
The Gain responses of above two approaches are exactly the same at all frequencies. 
Both transfer functions exhibit unit (0dB) gain at all frequencies.

The definition e^(−sTs), or the trapezoidal approximation, leads to a direct 
connection between the Z-transform of the digital domain and the Laplace transform 
of the analog domain.

Vern's OneNote Signals and Systems --> Introduction to Discrete-Time Systems

'''
def Mapping_of_z_transfer():
    
    fs = 100 * 10**3                  # Hz, Switching frequency fs of convertor.
    Ts = 1 / fs                       # Second, Switching period.
    f = np.linspace(1, 10*fs, 1000)   # Frequency sweeping range.
    omega = 2*np.pi*f                 # Sweeping range in rad/sec.
    
    # (cosx)^2 + (sinx)^2 = 1
    Gain_z_delay = np.sqrt((np.cos(omega*Ts))**2 + (np.sin(omega*Ts))**2)
    Phase_z_delay = np.arctan(np.sin(omega*Ts) / np.cos(omega*Ts))
    
    Gain_z_integ = np.sqrt(1 + (omega*Ts/2)**2) / np.sqrt(1 + (omega*Ts/2)**2)
    # Phase_z_integ = np.arctan((omega*Ts/2)) - np.arctan((-omega*Ts/2))
    # Calcuation above doesn't within +/-90 degree, use below one instead.
    Phase_z_integ = np.arctan(omega*Ts / (1 - (omega*Ts)**2 / 4))
    
    # Seems there is bug exist with np.log10, using Gain directly without dB.
    # Gain_z_delay = 20*np.log10(Gain_z_delay)
    # Gain_z_integ = 20*np.log10(Gain_z_integ)
    Phase_z_delay = 180*Phase_z_delay / np.pi
    Phase_z_integ = 180*Phase_z_integ / np.pi
    
    fig, ax = plt.subplots(2, 1, figsize = (11, 8))
    
    ax[0].grid(which="both", ls=':')
    ax[0].set_ylabel('Gain')
    
    ax[1].grid(which="both", ls=':')
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].set_ylabel('Phase [deg]')
    
    ax[0].semilogx(f, (Gain_z_delay), label = 'By Delay (z$^{-1}$)')
    ax[0].semilogx(f, (Gain_z_integ), label = 'Integral (s$^{-1}$)', linestyle = '-.')
    ax[1].semilogx(f, (Phase_z_delay), label = 'By Delay (z$^{-1}$)')
    ax[1].semilogx(f, (Phase_z_integ), label = 'Integral (s$^{-1}$)', linestyle = ':')
    
    ax[0].legend(loc=3, fontsize=16, shadow=True)
    ax[1].legend(loc=3, fontsize=16, shadow=True)
    fig.tight_layout()
    
Magnitude_domain_z_VS_s()
Mapping_of_z_transfer()
