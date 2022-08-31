# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 08:56:54 2022

@author: Vern Chen

Discrete-Time Compensator Design.

This script was converted from Fundamental of Power Electronics 19.3.2
Design Example of Matlab script.

"""

import math as math
import control.matlab as ctrl
import plotBode as pB

# Synchronous Buck converter parameters
Vg = 5; Vref = 1.8; D = Vref/Vg       # Input and reference voltages, duty cycle
L = 1e-6; RL = 30e-3                  # Inductance and series resistance
C = 200e-6; Resr = 0.8e-3             # Capacitance and capacitor ESR
fo = 1/(2*math.pi*math.sqrt(L*C))     # Pole frequency
R = 1000                              # Load resistance
fs = 1e6; Ts = 1/fs                   # Switching frequency and period

s = ctrl.tf('s')
z = ctrl.tf('z')                      # Define s and z

# Open-loop control to output transfer function
Gvd = Vg*(Resr+1/s/C)/(Resr + 1/s/C + s*L + RL)
fp2 = 1e6; H = 1/(1 + s/2/math.pi/fp2)# Sensor transfer function
Tu = H * Gvd                          # Uncompensated loop gain, no delay

# Analog PID compensator
fc = 100e3                            # Cross-over frequency
fL = 8e3; fz = 33e3; fp1 = 300e3;     # Corner frequencies
Gcm = math.sqrt(fz/fp1)*(fc/fo)**2/Vg # Mid-frequency gain

#Analog compensator transfer function
Gc = Gcm*(1 + 2*math.pi*fL/s)*(1 + s/2/math.pi/fz)/(1+s/2/math.pi/fp1);
T = Gc*Tu                             # Loop gain with analog compensator

# Uncompensated loop gain, including delay
td = D*Ts                             # Delay in the digital control loop
Tu.IODelay = td                       # Delay
Tud = ctrl.c2d(Tu,Ts,'impulse')       # Mapping of Tu with delay

# Analog PID compensator redesigned for digital implementation
fL = 8e3; fz = 22e3; fp1 = 450e3;     # Corner frequencies
Gcm = math.sqrt(fz/fp1)*(fc/fo)**2/Vg # Mid-frequency gain
Gca = Gcm*(1 + 2*math.pi*fL/s)*(1 + s/2/math.pi/fz)/(1+s/2/math.pi/fp1)

# Digital compensator transfer function
Gcd = ctrl.c2d(Gca, Ts, method='bilinear', prewarp_frequency=2*math.pi*fc)
Td = Tud*Gcd                          # Loop gain with digital compensator

# Compare magnitude and phase responses of T and Td
pB.Tf2Bodeplot((T, Td), ('Analog', 'Digital'))

'''
ctrl.bode(Gca, 
          ctrl.logspace(0, 7),
          dB = True, 
          Hz = True,
          deg = True, 
          plot = True)                # Bode plot of T

ctrl.bode(Gcd,
          ctrl.logspace(0, 7),
          dB = True, 
          Hz = True,
          deg = True, 
          plot = True)                # Bode plot of Td
'''
