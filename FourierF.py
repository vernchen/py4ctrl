# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 10:45:41 2022

Fourier Series and Fourier Transform
Square Wave Study

	• For periodic signals, the spectral coefficients have amplitudes C_n  and 
      occur at discrete set of harmonically related frequencies (nω_0), where, 
      (n = 0, ±1, ±2…).
	• For non-periodic signals, the complex exponentials occur at continuous 
      frequencies with magnitude [X(ω)dω/2π].

Vern's Noteboot Signals and Systems --> Fourier Series and Fourier Transform

@author: Vern Chen

"""

import numpy as np
import matplotlib.pyplot as plt

# Amplitude in time domain.
A = 1
# The fundamental frequency.
T = 1

'''
The exponential form of Fourier series of a continuous-time periodic signal x(t).

The set of coefficients [Cn] is called the set of the Fourier series 
coefficients or the spectral coefficients of signal x(t). The complex 
coefficients [C_n] measure the portion of the signal x(t), that is at 
each harmonic of the fundamantal component.

'''
def Periodic_Square(A, T):
    
    # For t axis, define length of range to include 10 period of square wave.
    p_c = 10
    
    # At time domain, define t to include 10 period of square wave.
    t = np.linspace(-p_c*T/2, p_c*T/2, 100*p_c+1)
    x_t = np.empty(100*p_c+1)
    
    # Generate square wave function.
    for i, value in enumerate(t):
        if t[i] % T < T/2:
            x_t[i] = A
        else:
            x_t[i] = -A
    
    # At frequency domain.
    # Total Harmonic numbers to be plotted.
    h_n = 25
    n = np.linspace(-h_n, h_n, 2*h_n+1)
    An = np.empty(2*h_n+1)
    Bn = np.empty(2*h_n+1)
    
    # The Fourier series coefficients of the function of x(t) are discrete in 
    # nature and hance we obtain a discrete spectrum.
    # Cn = An + j*Bn, An = 0 for this case.
    #            A 
    # Cn = -j * ———— * (1 - cos(n*pi))
    #           n*pi
    for i, value in enumerate(n):
        if n[i] == 0:
            An[i] = 0
            Bn[i] = 0
        else:
            An[i] = 0
            Bn[i] = -(1-np.cos(value*np.pi))*(A/value*np.pi)
    
    fig, ax = plt.subplots(2, 1, figsize = (11, 8))
    
    ax[0].set_xlabel('time [s]')
    ax[0].set_ylabel('x(t)')
    
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].set_ylabel('B$_n$')
    
    ax[0].plot(t, x_t)
    ax[1].stem(n/T, Bn)
    
'''
The Fourier Transform of a continuous-time non-periodic signal x(t) defined as X(w).

The Fourier transform X(ω)  of an aperiodic signal x(t) is called the spectrum 
of the signal x(t).

The equation of inverse Fourier transform plays a role of non-periodic signals 
similar to the equation of Fourier series for periodic signals. Because both 
the equations represent the linear combination of complex exponentials.

'''
def Non_Periodic_Square(A, T):
    
    # As non periodic doesn't have period, keep same number as periodic square.
    p_c = 10
    
    # At time domain.
    t = np.linspace(-p_c*T/2, p_c*T/2, 100*p_c+1)
    x_t = np.empty(100*p_c+1)
    
    # Generate square wave function.
    for i, value in enumerate(t):
        if -T < t[i] < T:
            x_t[i] = A
        else:
            x_t[i] = 0
    
    # At frequency domain.
    # Harmonic numbers.
    f = np.linspace(-10/T, 10/T, 1000)
    
    # The Fourier transfer.
    #         A * sin(2*pi*f*T)
    # X(jw) = —————————————————
    #              pi * f
    X_w = A * np.sin(2*np.pi*f*T) / (np.pi*f)
    
    fig, ax = plt.subplots(2, 1, figsize = (11, 8))
    
    ax[0].set_xlabel('time [s]')
    ax[0].set_ylabel('x(t)')
    
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].set_ylabel('X(j$\omega$)')
    
    ax[0].plot(t, x_t)
    ax[1].plot(f, X_w)
    
'''
Averaging is an artifice that facilitates the derivation of tractable equations 
describing the low-frequency dynamics of the switching converter. It removes the 
waveform components at the switching frequency and its harmonics, while preserving
the magnitude and phase of the waveform low-frequency components.

The averaging operator is a transformation that effectively performs a lowpass
function to remove the switching ripple. We can take the Laplace transformation:
    <X(s)>Ts = Gav(s) * X(s)
    
              e^(s*Ts/2) - e^(-s*Ts/2)
    Gav(s) = ——————————————————————————
                      s*Ts
                      
compute the effect of the averaging operator on a sinusoid of angular frequency 
ω by letting s = jω in the above equation:
               sin(w*Ts/2)
    Gav(jw) = —————————————
                 w*Ts/2
    
Vern's Noteboot Math --> Laplace Transfer of Average Approximation

'''
def Average_Approximation():
    # Consider frequency from 0.1 to 10 times of switching frequency.
    f = np.linspace(0.1/T, 10/T, 1000)
    fs = 1/T
    x = f / fs
    
    # Calculate magnitude of Gav.
    Gav = abs(np.sin(2*np.pi*f*T / 2) / (2*np.pi*f*T / 2))
    
    fig, ax = plt.subplots(1, 1, figsize = (11, 8))
    
    ax.grid(which="both")
    ax.set_xlim(0.1, 10)
    ax.set_ylim(-50, 0)
    
    ax.set_xlabel('f / f$_s$')
    ax.set_ylabel('Magnitude [dB]')
    
    ax.semilogx(x, 20*np.log10(Gav))

Periodic_Square(A, T)
Non_Periodic_Square(A, T)
Average_Approximation()
