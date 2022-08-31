# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:56:25 2022

The file name 'Resurrection F' is coming from Dragon Ball movie in 2015.
F stands for Frieza in the movie while here stands for transfer function in control.

@author: Vern Chen
"""

# To support local include style as C language.
# import glob    # File management and filtering.
import os      # Helps manage and create specific paths.
import sys     # System-specific parameters and functions.

# Get current working directory.
cwd = os.getcwd()
# print(cwd)

# Add relative path for import.
sys.path.append(cwd + '/tools/')

import math as math
import numpy as np
import sympy
import control.matlab as ctrl
import matplotlib.pyplot as plt

"""
Single Pole Response:
    
Consider the simple R–C low-pass filter
  ______[]______
 |      R      |   +
 |             |
(+)          C =  v2(s)
(-) v1(s)      |
 |             |
 |_____________|   -
 
The transfer function is given by the voltage divider ratio
G(s) = v2(s) / v1(s)
     = (1/sC) / (1/sC + R)
     = 1 / (1 + sRC)
     = 1 / (1 + s/w0)
     
w0 = 1 / RC

"""
def Single_Pole(R, C):
    # Circuit parameters init.
    # R = 1000      # 1k Ohm
    # C = 0.00001   # 10u F
    
    # Corner Frequency.
    omega_0 = 1 / (R * C)
    
    # Generate transfer function base on the circuit.
    single_p = ctrl.tf([1], [1/omega_0, 1])
    print(single_p)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 100kHz.
    gain, phase, omega = ctrl.bode(single_p, 
                                   ctrl.logspace(0, 5), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return single_p

"""
Single Zero Response:

A single zero response contains a root in the numerator of the transfer function
G(s) = 1 + s/w0

"""
def Single_Zero(omega_0):
    single_z = ctrl.tf([1/omega_0, 1], [1])
    print(single_z)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 100kHz.
    gain, phase, omega = ctrl.bode(single_z, 
                                   ctrl.logspace(0, 5), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return single_z

"""
Right Half-Plane Zero (RHP Zero):
 
Right half-plane zeroes are often encountered in the small-signal transfer functions 
of switching converters. These terms have the following normalized form:
G(s) = 1 - s/w0

"""
def RHP_Zero(omega_0):
    RHP_z = ctrl.tf([-1/omega_0, 1], [1])
    print(RHP_z)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 100kHz.
    gain, phase, omega = ctrl.bode(RHP_z, 
                                   ctrl.logspace(0, 5), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return RHP_z

"""
Inverted Pole:

The inverted pole has the transfer function
G(s) = 1 / (1 + w0/s)
     = (s/w0)/(1 + s/w0)

"""
def Inverted_Pole(omega_0):
    inverted_p = ctrl.tf([1/omega_0, 0], [1/omega_0, 1])
    print(inverted_p)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 100kHz.
    gain, phase, omega = ctrl.bode(inverted_p, 
                                   ctrl.logspace(0, 5), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return inverted_p

"""
Inverted Zero:

The inverted zero has the transfer function
G(s) = (1 + w0/s)
     = (1 + s/w0) / (s/w0)

"""
def Inverted_Zero(omega_0):
    inverted_z = ctrl.tf([1/omega_0, 1], [1/omega_0, 0])
    print(inverted_z)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 100kHz.
    gain, phase, omega = ctrl.bode(inverted_z, 
                                   ctrl.logspace(0, 5), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return inverted_z

"""
Quadratic Pole Response: Resonance:
    
Consider next the transfer function G(s) of the two-pole low-pass filter:
  _____ooo_____________
 |      L      |      |    +
 |             |      |
(+)          C =   R []    v2(s)
(-) v1(s)      |      |
 |             |      |
 |_____________|______|    -

The buck converter contains a filter of this type. When manipulated into 
canonical form, the models of the boost and buck–boost also contain similar
filters. One can show that the transfer function of this network is The 
transfer function is given by the voltage divider ratio:
G(s) = v2(s) / v1(s)
                1
     = --------------------
       (1 + L/R*s + LC*s^2)
       
                    1
     = ----------------------------
       1 + 2*zeta*(s/w0) + (s/w0)^2
       
                   1
     = -------------------------
       1 + (s/(Q*w0)) + (s/w0)^2

where:
    Q = 1 / (2*zeta)
The parameter Q is called the quality factor of the circuit, and is a measure 
of the dissipation in the system. A more general definition of Q, for sinusoidal 
excitation of a passive element or network, is
                  peak stored energy
    Q = 2*pi*-----------------------------
              energy dissipated per cycle
                                   ___ 
    f0 = w0 / (2*pi) = 1 / (2*pi* V LC)
           ____
    Q = R*V C/L
    
"""
def Quadratic_Pole(R, L, C):
    # R = 10        # 10 ohm
    # L = 0.00016   # 160u H
    # C = 0.00016   # 160u F
    
    omega_0 = (1 / math.sqrt(L*C))
    Q = R * math.sqrt(C/L)
    quadratic_p = ctrl.tf([1], [1/(omega_0**2), 1 / (Q*omega_0), 1])
    print(quadratic_p)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 100kHz.
    gain, phase, omega = ctrl.bode(quadratic_p, 
                                   ctrl.logspace(0, 5), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return quadratic_p

"""
Series Impedances: Addition of Asymptotes

Consider the series-connected R–C network
  ______[]______
        R      |
               |
Z(s)           = C
-->            |
               |
  _____________|
  
To construct the magnitude asymptotes of the total series impedance Z(s):
Z(s) = R + 1/sC
     = (sRC + 1) / sC
     = (s/w0 + 1) / (s/(R*w0))

where:
    w0 = 1/RC
"""
def Series_Impedances(R, C):
    # R = 10        # 10 ohm
    # C = 0.000001   # 1u F
    
    # The capacitor has an impedance magnitude of 1/ωC. 
    # This quantity varies inversely with ω, and hence its magnitude Bode plot 
    # is a line with slope −20 dB/decade. The line passes through 1Ω ⇒ 0dBΩ 
    # at the angular frequency ω where 1/ωC = 1 ohm
    omega = 1 / (1*C)
    f = omega / (2 * np.pi)
    print(f)
    
    # The corner frequency f0, where the asymptotes intersect, can now be 
    # easily deduced. At angular frequency ω0 = 2πf0, the two asymptotes are 
    # equal in value:
    omega_0 = 1 / (R*C)
    f0 = omega_0 / (2 * np.pi)
    print(f0)
    
    series_z = ctrl.tf([1/omega_0, 1], [1/(R*omega_0), 0])
    print(series_z)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(series_z, 
                                   ctrl.logspace(0, 6), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return series_z
    
"""
Series Resonant:
To construct the magnitude asymptotes for the series R–L–C circuit. 
  _____________
              |
              [] R
              |
Z(s)          Q  L
-->           Q 
              |
              =  C
  ____________|
  
The series impedance Z(s) is:
Z(s) = R + sL + 1/sC

"""
def Series_Resonant(R, L, C):
    # R = 10          # 10 ohm
    # L = 0.001       # 1m H
    # C = 0.0000001   # 0.1u F
    
    # The series impedance Z(s) is dominated by the capacitor at low frequency, 
    # by the resistor at mid frequencies, and by the inductor at high frequencies
    
    # The impedance Z(s) contains a zero at angular frequency ω1, where the 
    # capacitor and resistor asymptotes intersect. By equating the expressions 
    # for the resistor and capacitor asymptotes, we can find ω1
    omega_1 = 1 / (R*C)
    
    # A second zero occurs at angular frequency ω2, where the inductor and 
    # resistor asymptotes intersect. Upon equating the expressions for the 
    # resistor and inductor asymptotes at ω2, we obtain
    omega_2 = R / L
    
    # As R decreased, the approximate corner frequencies ω1 and ω2 move closer 
    # together until reducing R further in value causes the asymptotes to become 
    # independent of the value of R, The ||Z|| asymptotes now switch directly 
    # from ωL to 1/ωCAt corner frequency ω0, the inductor and capacitor 
    # asymptotes are equal in value ω0L = 1 / (ω0C) = R0
    omega_0 = 1 / math.sqrt(L*C)
    R0 = omega_0 * L
    
    # At ω = ω0, the inductor and capacitor impedances are equal in magnitude 
    # but opposite in phase. Hence, they exactly cancel out in the series 
    # impedance, and we are left with Z(s) = R
    # Z(jw0) = R + jw0L + 1 / (jW0C) = R + jR0 + R0/j = R + jR0 - jR0
    # The actual curve in the vicinity of the resonance at ω = ω0 can deviate
    # significantly from the asymptotes, because its value is determined by R 
    # rather than ωL or 1/ωC, the deviation of the actual curve from the 
    # asymptotes at ω = ω0 is equal to Q,
    Q = R0 / R
    
    print('w1 =', omega_1, 'w2 =', omega_2, 'w0 =', omega_0, '\nR0 =', R0, 'Q =', Q)
    
    series_r = ctrl.tf([L*C, R*C, 1], [C, 0])
    print(series_r)
    
    # Get the Gain and Phase and angular frequency 100Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(series_r, 
                                   ctrl.logspace(2, 6), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return series_r

"""
Parallel Resonant: Inverse Addition of Asymptotes

To construct the magnitude asymptotes for the parallel R–L–C circuit. 
  ____________________
       |      |      |
       |      |      |
Z(s)   Q L    = C   [] R
-->    Q      |      |
       |      |      |
  _____|______|______|

The parallel impedance Z(s) is:
Z(s) = 1 / (1/R + 1/sL + sC)
     = sL / (1 + sL/R + LC*s^2)
     
"""
def Parallel_Resonant(R, L, C):
    # R = 1000        # 1k ohm
    # L = 0.001       # 1m H
    # C = 0.0000001   # 0.1u F
    
    # The series impedance Z(s) is dominated by the capacitor at low frequency, 
    # by the resistor at mid frequencies, and by the inductor at high frequencies
    
    # The impedance Z(s) contains a zero at angular frequency ω1, where the 
    # capacitor and resistor asymptotes intersect. By equating the expressions 
    # for the resistor and capacitor asymptotes, we can find ω1
    omega_1 = R / L
    
    # A second zero occurs at angular frequency ω2, where the inductor and 
    # resistor asymptotes intersect. Upon equating the expressions for the 
    # resistor and inductor asymptotes at ω2, we obtain
    omega_2 = 1 / (R*C)
    
    # As R decreased, the approximate corner frequencies ω1 and ω2 move closer 
    # together until reducing R further in value causes the asymptotes to become 
    # independent of the value of R, The ||Z|| asymptotes now switch directly 
    # from ωL to 1/ωCAt corner frequency ω0, the inductor and capacitor 
    # asymptotes are equal in value ω0L = 1 / (ω0C) = R0
    omega_0 = 1 / math.sqrt(L*C)
    R0 = omega_0 * L
    
    # At ω = ω0, the inductor and capacitor impedances are equal in magnitude 
    # but opposite in phase. Hence, they exactly cancel out in the series 
    # impedance, and we are left with Z(s) = R
    # Z(jw0) = R + jw0L + 1 / (jW0C) = R + jR0 + R0/j = R + jR0 - jR0
    # The actual curve in the vicinity of the resonance at ω = ω0 can deviate
    # significantly from the asymptotes, because its value is determined by R 
    # rather than ωL or 1/ωC, the deviation of the actual curve from the 
    # asymptotes at ω = ω0 is equal to Q,
    Q = R / R0
    
    print('w1 =', omega_1, 'w2 =', omega_2, 'w0 =', omega_0, '\nR0 =', R0, 'Q =', Q)
    
    parallel_r = ctrl.tf([L, 0], [L*C, L/R, 1])
    # Converter Transfer Functions of Vern's notebook for dividor transfer function.
    # dividor_tf = ctrl.tf([1], [L*C, L/R, 1])
    print(parallel_r)
    
    # Get the Gain and Phase and angular frequency 100Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(parallel_r, 
                                   ctrl.logspace(2, 6), 
                                   dB = True, 
                                   Hz = True,
                                   deg = True, 
                                   plot = True)
    return parallel_r

"""
Chose pole and zero frequency based on phase lead of compensator requirment at 
frequency fc.
"""
def Zero_Frequency(fc, Theta):
    return fc * math.sqrt((1 - math.sin(Theta)) / (1 + math.sin(Theta)))

def Pole_Frequency(fc, Theta):
    return fc * math.sqrt((1 + math.sin(Theta)) / (1 - math.sin(Theta)))

"""
Lead (PD) compensator

This type of compensator transfer function is used to improve the phase margin. 
A zero is added to the loop gain, at a frequency fz sufficiently far below the 
crossover frequency fc.

The lead compensator is also called a proportional-plus-derivative, PD controller.

The transfer function of the lead compensator therefore contains a low-frequency 
zero and several high-frequency poles.

A side effect of the zero is that it causes the compensator gain to increase 
with frequency, with a +20 dB/decade slope.

Of particular concern are the switching frequency harmonics present in the 
output voltage and feedback signals. If the compensator gain at the switching 
frequency is too great, then these switching harmonics are amplified by the 
compensator, and can disrupt the operation of the pulse-width modulator. 
So the compensator network should contain poles at a frequency less than the 
switching frequency.

These considerations typically restrict the crossover frequency fc to be less 
than approximately 10% of the converter switching frequency fs.

Gc(s) = Gc0 * (1 + s/wz) / (1 + s/wp)

The maximum phase occurs at a frequency fmax given by the geometrical mean of 
pole and zero frequencies
        _____
fmax = V fzfp

To obtain the maximum improvement in phase margin, we should design our 
compensator so that the frequency fϕmax coincides with the loop gain crossover 
frequency fc

"""
def Lead_Compensator(fz, fp):
    omega_z = 2 * np.pi * fz
    omega_p = 2 * np.pi * fp
    
    # Gc0 = math.sqrt(fz / fp)    # Gain value at maximum phase.
    # phase_max = math.atan(0.5 / Gc0 - 0.5 * Gc0)
    
    fp_div_fz = np.linspace(1, 1000)
    phase_max = np.arctan(0.5 * np.sqrt(fp_div_fz) - 0.5 / np.sqrt(fp_div_fz)) 
    phase_max = phase_max * 180 / np.pi

    plt.grid(which="both", ls=':')
    plt.xlabel('fp/fz')
    plt.ylabel('Phase [deg]')
    plt.semilogx(fp_div_fz, phase_max)
    plt.show()
    
    fmax = math.sqrt(fz * fp)
    print('Maximum phase at', fmax, 'Hz')
    
    lead_comp = ctrl.tf([1/omega_z, 1], [1/omega_p, 1])
    print(lead_comp)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(lead_comp, 
                               ctrl.logspace(1, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    return lead_comp

"""
Lag (PI) Compensator

This type of compensator is used to increase the low-frequency loop gain, such 
that the output is better regulated at dc and at frequencies well below the 
loop crossover frequency. An inverted zero is added to the loop gain, at 
frequency fl.

Gc(s) = Gc∞ * (1 + wl/s)
      = Gc∞ * (s/wl + 1) / (s/wl)
      
If fl is sufficiently lower than the loop crossover frequency fc, then the 
phase margin is unchanged. This type of compensator is also called a 
proportional-plus-integral, or PI controller. At low frequencies, the inverted 
zero causes the compensator to integrate the error signal.

"""
def Lag_Compensator(fl):
    omega_l = 2 * np.pi * fl
    
    lag_comp = ctrl.tf([1/omega_l, 1], [1/omega_l, 0])
    print(lag_comp)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(lag_comp, 
                               ctrl.logspace(1, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    return lag_comp

"""
Combined (PID) Compensator

The advantages of the lead and lag compensators can be combined, to obtain both 
wide bandwidth and zero steady-state error.

               (1 + wl/s) * (1 + s/wz)
Gc(s) = Gcm * -------------------------
              (1 + s/wp1) * (1 + s/wp2)
         
The inverted zero at frequency fl functions in the same manner as the PI 
compensator. The zero at frequency fz adds phase lead in the vicinity of the 
crossover frequency, as in the PD compensator. The high-frequency poles at 
frequencies fp1 and fp2 must be present in practical compensators, to cause the 
gain to roll off at high frequencies and to prevent the switching ripple from 
disrupting the operation of the pulse-width modulator.

The loop gain crossover frequency fc is chosen to be greater than fl and fz, 
but less than fp1 and fp2

(fl, fz) < fc < (fp1, fp2)

"""
def PID_Compensator(fl, fz, fp1, fp2):
    omega_l = 2 * np.pi * fl
    omega_z = 2 * np.pi * fz
    omega_p1 = 2 * np.pi * fp1
    omega_p2 = 2 * np.pi * fp2
    
    # Create a variable 's' to allow algebra operations for SISO systems
    s = ctrl.tf('s')
    pid_comp = (1 + omega_l/s) * (1 + s/omega_z) / ((1 + s/omega_p1) * (1 + s/omega_p2))
    print('Gpid =', pid_comp)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(pid_comp, 
                               ctrl.logspace(1, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True,
                               color = 'k',
                               linestyle = '--')
    return pid_comp

"""
Parallel Impedance

Total impedance when two device connected in parallel.
"""
def Parallel_Impedance(z1, z2):
    return 1 / (1/z1 + 1/z2)

"""
The Relationship Between Phase Margin and Closed-Loop Damping Factor
Let us consider a loop gain T(s) which is well-approximated, in the vicinity 
of the crossover frequency, by the following function
              1                   1
T(s) = ----------------- = ----------------
       (s/w0)(1 + s/w2)    s/w0 + s^2/w0w2
       
  T(s)              1                     1
-------- = -------------------- = -------------------
1 + T(s)   1 + (s/w0)(1 + s/w2)   1 + s/w0 + s^2/w0w2

Looks carefully at above two equation, and the difference is denominator added 
"1 +" which will majorly affect at low frequency gain (= 1, 0dB).
      _____
wc = V w0w2 = 2*pi*fc
             ______
Q = w0/wc = V w0/w2

We can obtain such a relationship by finding the frequency at which the magnitude
of T is exactly equal to unity. We then evaluate the exact phase of T at this 
frequency, and compute the phase margin. This phase margin is a function of the 
ratio f0/ f2, or Q2.We can then solve to find Q as a function of the phase margin.
     ________
Q = V cos(ϕm) / sin(ϕm)

"""
def PhaseMargin2DampingFactor(f0, f2):
    '''
    omega_0 = 2 * np.pi * f0
    omega_2 = 2 * np.pi * f2
    
    Ts = ctrl.tf([1], [1/(omega_0*omega_2), 1/omega_0, 0])
    
    Gs = Ts / (1 + Ts)
    
    omega_c = math.sqrt(omega_0 * omega_2)
    Q = math.sqrt(omega_0 / omega_2)
    print('Cross frequency of Gs is around', (omega_c / 2 * np.pi),
          '\nQ-factory is', Q)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(Ts, 
                               ctrl.logspace(1, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True,
                               color = 'b')
    
    gain, phase, omega = ctrl.bode(Gs, 
                               ctrl.logspace(1, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True,
                               color = 'g')
    '''
    Phi = np.linspace(1, 90)          # φ in degree.
    Q = 20 * np.log10(np.sqrt(np.cos(Phi * np.pi / 180)) / np.sin(Phi * np.pi / 180))
    plt.grid(which="both", ls=':')
    plt.xlabel('Phase Margin φ (degree)')
    plt.ylabel('Q-factor (dB)')
    plt.ylim(-30, 30)
    plt.plot(Phi, Q)
    plt.show()
    
"""
Transient Response vs. Damping Factor
One can solve for the unit-step response of the T/(1 + T) transfer function, 
by multiplying below equation by 1/s and then taking the inverse Laplace transform.
       
  T(s)      1            1              1             1              1
-------- * -- = -------------------- * -- = --------------------- * --
1 + T(s)   s    1 + (s/w0)(1 + s/w2)   s    1 + s/Qwc + (s/wc)^2    s

"""
def Transient2DampingFactor(Q, omega_c):
    s, t = sympy.symbols('s, t')
    # Q = sympy.symbols('Q', real = True,  positive=True)
    # omega_c = sympy.symbols('omega_c', real = True,  positive=True)
    # expression = 1 / (s * (1 + s/(Q*omega_c) + (s/omega_c)**2))
    a = sympy.symbols('a', real = True,  positive=True)
    expression = 1 / (a + s)
    
    f = sympy.inverse_laplace_transform(expression, s, t)
    print(f)
    
"""
Load Step Response vs. Damping Factor
we also are interested in the response of the output voltage to a step change 
in load current. Let us consider the case where the closed-loop output impedance 
can be well approximated by a second-order function of the form
       
         1            sR/wc           1
Zo(s) * -- = --------------------- * --
        s    1 + s/Qwc + (s/wc)^2    s

"""
def Loadstep2DampingFactor(Q, omega_c, R):
    s, t = sympy.symbols('s, t')
    # Q = sympy.symbols('Q', real = True,  positive=True)
    # omega_c = sympy.symbols('omega_c', real = True,  positive=True)
    # expression = (s * R / omega_c) / (s * (1 + s/(Q*omega_c) + (s/omega_c)**2))
    a = sympy.symbols('a', real = True,  positive=True)
    expression = 1 / (a + s)**2
    
    f = sympy.inverse_laplace_transform(expression, s, t)
    print(f)

# Single_Pole(1000, 0.00001)
# Single_Zero(100)
# RHP_Zero(1000)
# Inverted_Pole(1000)
# Inverted_Zero(1000)
# Quadratic_Pole(10, 0.00016, 0.00016)
# Series_Impedances(10, 0.000001)
# Series_Resonant(10, 0.001, 0.0000001)
# Parallel_Resonant(1000, 0.001, 0.0000001)
# Lead_Compensator(100, 10000)
# Lag_Compensator(1000)
# PID_Compensator(500, 1000, 5000, 10000)
# PhaseMargin2DampingFactor(10, 1000)

