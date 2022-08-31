# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:58:30 2022

In this Module, Basic buck, boost, and buck–boost converters with current 
programmed control are summarized.

Control-to-output and line-to-output transfer functions for more accurate model.

For each case, the salient features are expressed based on the corresponding quantity 
of duty-cycle control, multiplied by a factor that accounts for current-programmed 
control.

Vern's OneNote Power Electronics --> Current-Programmed Control

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

import numpy as np
import math as math
import control.matlab as ctrl
import matplotlib.pyplot as plt
import plotBode as pB

# Global Variables definition.
Switching_Freq = 100e3          # Hz, Switching frequency fs of convertor.
Ts = 1 / Switching_Freq               # Second, Switching period.

"""
Duty to Conversion rate analysis of CCM and DCM

No matter in CCM or DCM, the output voltage with certain load R was depending on
output current DC I, as V = IR, while Δi dose not depend on load R, convertor will
enter DCM when Δi > I.

Once steady state achieved in CCM, I keep no change, duty cycle depends only on
input / output volage, for keeping Δi fo full duty = 0.

"""
def CCM_DCM_Conversion_Rate(R, L, type_c):
    
    D = np.linspace(0.01, 0.99)             # Duty cycle range from 1 ~ 99%.
    
    K = 2*L/(R*Ts)                          # Dimensionless parameter.
    
    # Klightload = 2*L/(10*R*Ts)
    # Koverload = 2*L/(0.1*R*Ts)
    
    if type_c == 'buck':
        # The critical value of K at the boundary between CCM and DCM.
        Kcrit = 1 - D
        # The critical value of K at the boundary between CCM and DCM.
        Rcrit = 2*L / ((1 - D)*Ts)
        # Equilibrium conversion ratio in CCM.
        Mccm = D
        # Equilibrium conversion ratio in DCM.
        Mdcm = 2 / (1 + np.sqrt(1 + 4*K/D**2))
        
    elif type_c == 'boost':
        '''
        For boost convertor, CCM can be achieved with small duty becuase M is samll,
        voltage applied on inductor (V -Vg) is also too small to decrease i to zero.
        '''
        # The critical value of K at the boundary between CCM and DCM.
        Kcrit = D * (1 - D)**2
        # The critical value of K at the boundary between CCM and DCM.
        Rcrit = 2*L / (D * Ts * (1 - D)**2)
        '''
        DCM or CCM can be achieved within full range of D (0 ~ 100%), but getting 
        different output voltage depends on M(D) or M(D, K).
        
        How to explain the curve across of DCM and CCM? which means same duty of
        DCM and CCM can achieve same conversion ratio?
        
        The key point is average cuurent is much different, during CCM, there is
        high DC value of I, and D is independ to R, as well as load condition.
        So, with load decrease (increase R), convertor enters DCM at cross point, 
        then as load continue decrease, slope of Mdcm increase, as well as Mdcm
        for certain D, to maintain same gain, duty must start decrease, on the 
        opposite side, with load increase (small R), Mdcm decrease, and duty keep 
        increase, till enter CCM.
        
        Noted that with (heavy load) small R, there is a maximum limit of M for 
        DCM boost, because with DCM, energy won't accumulated cross (trough DC I) 
        switching cycles (i drop to zero for each cycle), AS Ts is fixed, so 
        there is maximum limit of di of 1st interval, energy was limited.But when 
        increase R to make K smaller, the M of DCM could increase significantly 
        larger than M of CCM.
        '''
        # Equilibrium conversion ratio in CCM.
        Mccm = 1 / (1 - D)
        # Equilibrium conversion ratio in DCM.
        Mdcm = (1 + np.sqrt(1 + 4*D**2/K)) / 2
        
    elif type_c =='buck-boost':
        # The critical value of K at the boundary between CCM and DCM.
        Kcrit = (1 - D)**2
        # The critical value of R at the boundary between CCM and DCM.
        Rcrit = 2*L / (Ts * (1 - D)**2)
        # Equilibrium conversion ratio in CCM.
        Mccm = -D / (1 - D)
        # Equilibrium conversion ratio in DCM.
        Mdcm = -D / np.sqrt(K)
        
    else:
        print('Convertor type', type_c, 'does not support!')
        
    fig, ax = plt.subplots(figsize = (11, 8))
    ax.grid(which="both", ls=':')
    ax.set_xlabel('D, duty cycle')
    ax.set_ylabel('M, equilibrium conversion ratio')
    ax.set_ylim(0, 10)
    
    ax.plot(D, Mccm, label = 'M$_{ccm}$')
    ax.plot(D, Mdcm, label = 'M$_{dcm}$')
    
    ax.plot(D, Kcrit, label = 'K$_{crit}$')
    ax.plot(D, Rcrit, label = 'R$_{crit}$')
    ax.legend(loc=6, fontsize=16, shadow=True)
    # Draw plot.
    plt.show()

"""
Control duty according to output voltage usually called Voltage Mode Control.
"""
def Voltage_Mode_TF(omega_0, omega_z, Q, Gvg0, Gvd0):
    
    # Create a variable 's' to allow algebra operations for SISO systems
    s = ctrl.tf('s')
    
    # Denominator polynomial of transfer function.
    den = 1 + s / (Q * omega_0) + (s / omega_0)**2
    
    # Line to output transfer function.
    Gvg = Gvg0 / den
    
    # Control to output transfer function, notice the RHZ.
    if math.isinf(omega_z):
        Gvd = Gvd0 / den
    else:
        Gvd = Gvd0 * (1 - s/omega_z) / den
        
    return Gvg, Gvd, den

"""
Duty-to-inductor current and line-to-inductor current transfer function.
"""
def Inductor_Current_TF(Gig0, Gid0, omega_z_ig, omega_z_id, den):
    
    # Create a variable 's' to allow algebra operations for SISO systems
    s = ctrl.tf('s')
    
    # Current-programmed closed-loop control-to-output transfer function.
    Gig = Gig0 * (1 + s / omega_z_ig) / den
    
    # Current-programmed line-to-output transfer function
    Gid = Gid0 * (1 + s / omega_z_id) / den
    
    return Gig, Gid

"""
Instead control duty directly, current-programmed mode control inductor current
and duty cycle depends on inductor current and control target current, as well 
as capacitor voltage and input voltage.

An advantage is simpler dynamics, as CPM control move one of the pole to high
frequency, so that wide-bandwith output voltage control can be achieve even without
PD compensator networks

Another advantage is CPM handles naturally on transformer saturation problems
in isolated converters (full-bridge or push-pull).
"""
def Current_Mode_TF(omega_c, omega_z, omega_gz, Qc, Gg0, Gc0):
    
    # Create a variable 's' to allow algebra operations for SISO systems
    s = ctrl.tf('s')
    
    # Denominator polynomial of transfer function.
    den_c = 1 + s / (Qc * omega_c) + (s / omega_c)**2
    
    # Line to output transfer function.
    if math.isinf(omega_gz):
        Gvg_cpm = Gg0 / den_c
    else:
        Gvg_cpm = Gg0 * (1 + s/omega_gz)/ den_c
        
    # Control to output transfer function, notice the RHZ.
    if math.isinf(omega_z):
        Gvc = Gc0 / den_c
    else:
        Gvc = Gc0 * (1 - s/omega_z) / den_c
        
    return Gvg_cpm, Gvc, den_c

"""


"""
def Buck_Plant_TF(R, L, C, Vg, V):
    D = V / Vg                              # Duty cycle by voltage-second balance.
    Gvg0 = D                                # DC gain of line-to-output TF.
    Gvd0 = V / D                            # DC gain of duty-to-output TF.
    omega_0 = 1 / math.sqrt(L * C)          # Corner frequency in rad/sec.
    Q = R * math.sqrt(C / L)                # Q factory of voltage response.
    omega_z = math.inf                      # Zero of voltage response.
    
    # measure of the tendency to DCM.
    K = 2*L/(R*Ts)                          # Dimensionless parameter.
    Rcrit = 2*L / ((1 - D)*Ts)              # Critical R for enter DCM.
    print('Enters DCM around <', V/Rcrit, 'A \nK = ', K, '\nRcrit = ', Rcrit)
    
    # Voltage mode transfer functions.
    Gvg, Gvd, den = Voltage_Mode_TF(omega_0, omega_z, Q, Gvg0, Gvd0)
    
    # Current-Programmed Control.
    Gig0 = D / R                            # DC gain of line-to-current TF.
    Gid0 = V / (D*R)                        # DC gain of duty-to-current TF.
    omega_z_ig = 1 / (R * C)                # Zero of line-to-current TF.
    omega_z_id = omega_z_ig                 # Zero of duty-to-current TF.
    
    # Line and duty to inductor transfer functions.
    Gig, Gid = Inductor_Current_TF(Gig0, Gid0, omega_z_ig, omega_z_id, den)
    
    M1 = (Vg - V) / L                       # Slope during 1st interval.
    M2 = V / L                              # Slope during 2nd interval.
    
    # Slope of artificial ramp.
    if M1 > M2:
        Ma = 0                              # Duty can be limited to < 50%.
    else:
        Ma = M2/2                           # Minimum artificial ramp when duty > 50%.
        
    Fg = D * (1 - D) * Ts / (2 * L)         # Line correct factor to duty.
    # Fv = 0                                  # Output voltage correct factor to duty.
    Fm = 1 / ((Ma + (M1 - M2)/2) * Ts)      # CPM modulator gain Fm
    
    # Correct factors for current mode transfer function.
    iFactor = 1 + Fm * V / (D * R)
    omega_0_factor = math.sqrt(iFactor)
    omega_gz = math.inf
    Q_factor = omega_0_factor / (1 + R*C*Fm*V / (D*L))
    Gvc_cpm_factor = Fm / iFactor
    Gvg_cpm_factor = (1 - Fm*Fg*V / D**2) / iFactor
    
    # Apply correct factors to get current mode parameters.
    omega_c = omega_0 *omega_0_factor
    Qc = Q * Q_factor
    Gc0 = Gvd0 * Gvc_cpm_factor
    Gg0 = Gvg0 * Gvg_cpm_factor
    
    # Curren mode transfer functions.
    Gvg_cpm, Gvc, den_c = Current_Mode_TF(omega_c, omega_z, omega_gz, Qc, Gg0, Gc0)
    
    return Gvg, Gvd, Gig, Gid, Gvg_cpm, Gvc


def Boost_Plant_TF(R, L, C, Vg, V):
    D = (V - Vg) / V                        # Duty cycle by voltage-second balance.
    Gvg0 = 1 / (1 - D)                      # DC gain of line-to-output TF.
    Gvd0 = V / (1 - D)                      # DC gain of duty-to-output TF.
    omega_0 = (1 - D) / math.sqrt(L * C)    # Corner frequency in rad/sec.
    Q = (1 - D) * R * math.sqrt(C / L)      # Q factory of voltage response.
    omega_z = ((1 - D)**2) * R / L          # Zero of voltage response.
    
    # measure of the tendency to DCM.
    K = 2*L/(R*Ts)                          # Dimensionless parameter.
    Rcrit = 2*L / (D*Ts*(1 - D)**2)         # Critical R for enter DCM.
    print('\nEnters DCM around <', V/Rcrit, 'A \nD = ', D, '\nK = ', K, '\nRcrit = ', Rcrit)
    
    Gvg, Gvd, den = Voltage_Mode_TF(omega_0, omega_z, Q, Gvg0, Gvd0)
    
    # Current-Programmed Control
    Gig0 = 1 / (((1 - D)**2)*R)             # DC gain of line-to-current TF.
    Gid0 = 2 * V * Gig0                     # DC gain of duty-to-current TF.
    omega_z_ig = 1 / (R * C)                # Zero of line-to-current TF.
    omega_z_id = 2 * omega_z_ig             # Zero of duty-to-current TF.
    
    # Line and duty to inductor transfer functions.
    Gig, Gid = Inductor_Current_TF(Gig0, Gid0, omega_z_ig, omega_z_id, den)
    
    M1 = Vg / L                             # Slope during 1st interval.
    M2 = (V - Vg) / L                       # Slope during 2nd interval.
    
    # Slope of artificial ramp.
    if M1 > M2:
        Ma = 0                              # Duty can be limited to < 50%.
    else:
        Ma = M2/2                           # Minimum artificial ramp when duty > 50%.
        
    # Fg = 0                                  # Line correct factor to duty.
    Fv = D * (1 - D) * Ts / (2 * L)         # Output voltage correct factor to duty.
    Fm = 1 / ((Ma + (M1 - M2)/2) * Ts)      # CPM modulator gain Fm
    
    # Correct factors for current mode transfer function.
    iFactor = 1 + 2*Fm*V / (((1 - D)**2)*R) + Fm*Fv*V / (1 - D)
    omega_0_factor = math.sqrt(iFactor)
    omega_gz = math.inf
    Q_factor = omega_0_factor / (1 + R*C*Fm*V / L - Fm*Fv*V / (1 - D))
    Gvc_cpm_factor = Fm / iFactor
    Gvg_cpm_factor = (1 + Fm*V / (((1 - D)**2)*R)) / iFactor
    
    # Apply correct factors to get current mode parameters.
    omega_c = omega_0 *omega_0_factor
    Gc0 = Gvd0 * Gvc_cpm_factor
    Qc = Q * Q_factor
    Gg0 = Gvg0 * Gvg_cpm_factor
    
    # Curren mode transfer functions.
    Gvg_cpm, Gvc, den_c = Current_Mode_TF(omega_c, omega_z, omega_gz, Qc, Gg0, Gc0)
    
    return Gvg, Gvd, Gig, Gid, Gvg_cpm, Gvc


def Buck_Boost_Plant_TF(R, L, C, Vg, V):
    D = V / (V - Vg)                        # Duty cycle by voltage-second balance.
    Gvg0 = -D / (1 - D)                     # DC gain of line-to-output TF.
    Gvd0 = V / (D * (1 - D))                # DC gain of duty-to-output TF.
    omega_0 = (1 - D) / math.sqrt(L * C)    # Corner frequency in rad/sec.
    Q = (1 - D) * R * math.sqrt(C / L)      # Q factory of voltage response.
    omega_z =  ((1 - D)**2) * R / (D * L)   # Zero of voltage response.
    
    # measure of the tendency to DCM.
    K = 2*L/(R*Ts)                          # Dimensionless parameter.
    Rcrit = 2*L / (Ts*(1 - D)**2)           # Critical R for enter DCM.
    print('Enters DCM around <', V/Rcrit, 'A \nK = ', K, '\nRcrit = ', Rcrit)
    
    Gvg, Gvd, den = Voltage_Mode_TF(omega_0, omega_z, Q, Gvg0, Gvd0)
    
    # Current-Programmed Control
    Gig0 = D / (R*(1-D)**2)                 # DC gain of line-to-current TF.
    Gid0 = abs(V)*(1+D)/(D*R*(1-D)**2)      # DC gain of duty-to-current TF.
    omega_z_ig = 1 / (R * C)                # Zero of line-to-current TF.
    omega_z_id = (1 + D) * omega_z_ig       # Zero of duty-to-current TF.
    
    # Line and duty to inductor transfer functions.
    Gig, Gid = Inductor_Current_TF(Gig0, Gid0, omega_z_ig, omega_z_id, den)
    
    # Current-Programmed Control
    M1 = Vg / L                             # Slope during 1st interval.
    M2 = - V / L                            # Slope during 2nd interval.
    
    # Slope of artificial ramp.
    if M1 > M2:
        Ma = 0                              # Duty can be limited to < 50%.
    else:
        Ma = M2/2                           # Minimum artificial ramp when duty > 50%.
        
    Fg = D * (1 - D) * Ts / (2 * L)         # Line correct factor to duty.
    Fv = -D * (1 - D) * Ts / (2 * L)        # Output voltage correct factor to duty.
    Fm = 1 / ((Ma + (M1 - M2)/2) * Ts)      # CPM modulator gain Fm
    
    # Correct factors for current mode transfer function.
    iFactor = 1 + Fm*abs(V)*(1 + D) / (D*((1 - D)**2)*R) - Fm*Fv*abs(V) / (D*(1 - D))
    omega_0_factor = math.sqrt(iFactor)
    omega_gz_factor = 1 + Fm*abs(V) / R*(1 - D)**2 - Fm*Fg*abs(V) / D**2
    omega_gz = (D*R*(1 - D)**2 / (abs(V)*L*Fm*Fg)) * omega_gz_factor
    Q_factor = omega_0_factor / (1 + R*C*Fm*abs(V) / (D*L) + Fm*Fv*abs(V) / (1 - D))
    Gvc_cpm_factor = Fm / iFactor
    Gvg_cpm_factor = -(D/(1 - D)) * (omega_gz_factor / iFactor)
    
    # Apply correct factors to get current mode parameters.
    omega_c = omega_0 *omega_0_factor
    Gc0 = Gvd0 * Gvc_cpm_factor
    Qc = Q * Q_factor
    Gg0 = Gvg0 * Gvg_cpm_factor
    
    # Curren mode transfer functions.
    Gvg_cpm, Gvc, den_c = Current_Mode_TF(omega_c, omega_z, omega_gz, Qc, Gg0, Gc0)
    
    return Gvg, Gvd, Gig, Gid, Gvg_cpm, Gvc


"""
Take 1400W buck design as example:
    
    Vg = 48V
    V = 12.2V
    
    Full load impedance
    R = 12.2^2 / 1400 = 0.1 ohm
    
    initial inductance L0= 5.3uH, 
    full load inductance Lf = 4.2uH at 1400W / 12.2
    
    Output capacitance C = 90uF * 37


Gvg, Gvd, Gig, Gid, Gvg_cpm, Gvc = Buck_Plant_TF(0.1, 4.2e-6, 37*90e-6, 48, 12.2)

G = (Gvg, Gvd, Gig, Gid, Gvg_cpm, Gvc)
Gstring = ('Gvg', 'Gvd', 'Gig', 'Gid', 'Gvg_cpm', 'Gvc')

Tf2Bodeplot(G, Gstring)
"""

"""
Take 360W boost design as example:
    
    Vg = 12.5V
    V = 56V
    
    Full load impedance
    R = 56^2 / 360 = 8.7 ohm
    
    Inductance L = 6uH, 
    
    Output capacitance C = 150uF * 3

"""

CCM_DCM_Conversion_Rate(8.7, 6e-6, 'boost')

Gvg, Gvd, Gig, Gid, Gvg_cpm, Gvc = Boost_Plant_TF(8.7, 6e-6, 3*150e-6, 12.2, 56)

G = (Gvg, Gvd, Gig, Gid, Gvg_cpm, Gvc)
Gstring = ('G$_{vg}$', 'G$_{vd}$', 'G$_{ig}$', 'G$_{id}$', 'G$_{vgcpm}$', 'G$_{vc}$')

pB.Tf2Bodeplot(G, Gstring)


