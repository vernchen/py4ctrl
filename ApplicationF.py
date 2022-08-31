# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:19:16 2022

Transfer function for actual cicurt models.

General theorem can be found at Vern's Noteboot --> Power Electronics.

@author: Vern Chen

"""
import math as math
import numpy as np
import control.matlab as ctrl
import ResurrectionF as resf

"""
Op Amp PD Compensator:
    
    C    R2
  ——||———[]——
 |          |     R3
 |————[]————|————[]—————————————————
 |    R1    |           |     |    +
 |          | - _       |     |
(+)vin      ———| --_    |  RL[]   vout
(-)            |   _>———|     |
 |          ———|_--     |     |
 |          | +         |     |    -
 ———————————————————————————————————

let us analyze the op amp circuit with an ideal op amp, this lead-lag circuit 
exhibits a transfer function having a zero and pole, and is suitable as a PD 
compensator in feedback loops requiring improvement of phase margin.

The design was based on null double injection theorem.

Vern's OneNote Power Electronics --> The Feedback Theorem

"""
def OpAmp_Compensator():
    R1 = 1600                         # ohm
    R2 = 16                           # ohm
    R3 = 1600                         # ohm
    RL = 100                          # ohm
    C = 0.1 * 10**-6                  # F
    
    '''
    Model the op amp using the equivalent circuit.The positive and negative input 
    ports are modeled with infinite input impedance, and a Thevenin-equivalent 
    circuit models the output port.
    
    ——————                ————[]——————
    v-                   |    R0    +
                       +_|_ 
        Gop(s)(v+ - v-)|__|        vout
                      - |
    v+                  |          -
    ——————              —————————————
    
    '''
    R0 = 50                           # ohm, Op Amp output impedence.
    f1 = 10                           # Hz, Op Amp cut off frequency.
    Gop0 = 100000                     # DC gain of Op Amp.
    
    omega_1 = 2*np.pi*f1              # rad/sec, Radian frequency.
    # Create a variable 's' to allow algebra operations for SISO systems
    s = ctrl.tf('s')
    
    # The openloop control-to-output transfer function.
    Gop = Gop0 / (1 + s/omega_1)
    print('Gop =', Gop)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
    gain, phase, omega = ctrl.bode(Gop, 
                               ctrl.logspace(1, 7), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    '''
    To apply the feedback theorem, we first identify an ideal injection point. 
    The error signal of this op amp feedback circuit can be taken to be the 
    op amp differential input voltage, and hence we can employ voltage injection 
    immediately following op amp model.
    This will cause vy to be directly proportional to the error signal (v+ - v-).
        C    R2
      ——||———[]——
     |          |     R3                        ifb
     |————[]————|————[]————————————————————————<—————————————————
     |    R1    |                       vz      R0    |     |    +
     |          |            |————————(- +)—————[]——>—      |
 vin(+)      v- o          +_|_     -      -       iR0   RL []   vout
    (-)           -Gop(s)v-|__|     vy     vx               |
     |       v+ o         - |                               |
     |          |           |       +      +                |    -
     —————————————————————————————————————————————————————————————
     
     The ideal forward gain Ginf is found by null vy.the dependent voltage source
     −Gopv− is also nulled, which implies that v− is nulled. Hence, the current ifb
     can be related to the input voltage vin as follows:
             
                   Vin - v-              Vin
        ifb = - ——————————————— = - ———————————————
               R1||(R2 + 1/sC)     R1||(R2 + 1/sC)
                    
    The null condition also allows us to easily relate the output voltage vout to ifb
    
        Vout = v- + ifb*R3 = ifb*R3
        
    '''
    omega_2 = 1 / ((R1 + R2)*C)       # Zero frequency.
    omega_3 = 1 / (R2*C)              # Pole frequency.
    f2 = omega_2 / (2*np.pi)          # Zero frequency.
    f3 = omega_3 / (2*np.pi)          # Pole frequency.
    
    Ginf = - (R3 /R1) * (1 + s/omega_2)/(1 + s/omega_3)
    print('Ginf =', Ginf,
          '\nZero at', f2, 'Pole at', f3)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
    gain, phase, omega = ctrl.bode(Ginf, 
                               ctrl.logspace(1, 7), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    '''
    The direct forward transmission gain G0(s) is found as defined by null vx.
    In the presence of the input vin we adjust the injection source vz such that 
    vx is nulled. Under these conditions, the dependent voltage source −Gopv− 
    does not influence the output:
        
                           R0||RL
        G0(s) = —————————————————————————————
                R0||RL + R3 + R1||(R2 + 1/sC)
                
                   R0||RL                 1 + sC(R1 + R2)
              = ———————————————— * ——————————————————————————————
                R1 + R3 + R0||RL   1 + sC(R2 + R1||(R3 + R0||RL))
        
        C    R2
      ——||———[]——
     |          |     R3  
     |————[]————|————[]———————————————————————
     |    R1                 |        |      +
     |                       |        |
 vin(+)                   R0 []   RL []     vout
    (-)                      |        |
     |                       |        |
     |                       |        |      -
     —————————————————————————————————————————
     
    '''
    G00 = (R0*RL / (R0 + RL)) / (R1 + R3 + (R0*RL / (R0 + RL)))
    omega_3 = omega_2
    omega_4 = 1 / (C * (R2 + resf.Parallel_Impedance(R1, (R3 + resf.Parallel_Impedance(R0, RL)))))
    f3 = f2
    f4 = omega_4 / (2*np.pi)
    
    G0 = G00 * (1 + s/omega_3) / (1 + s/omega_4)
    print('G0 =', G0,
          '\nDC gain of G0 =', 20*np.log10(G00), 'dB\nPole at', f4)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
    gain, phase, omega = ctrl.bode(G0, 
                               ctrl.logspace(1, 7), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    '''
    Null loop gain Tn(s): In the presence of the input vin, the signal vz is 
    injected and is adjusted as necessary to null the output vout.
             vy(s)
    Tn(s) = ———————  with vout null to 0
             vx(s)
    The null condition implies that there is no voltage across the load resistor 
    RL and hence there is no current through the load resistor. The op amp 
    output current is
    
        iR0 = vx / R0
    
    Since the load current is zero, the current iR0 flows through R3. Since the
    load voltage is zero, we can express v− as:
        
        v- = -iR0 * R3
    
    The voltage vy is related to v− by the op amp gain Gop:
    
        vy = Gop(s)*v-
    
    Hence, we can express the null loop gain as:
        
                iR0    v-    vy    1
        Tn(s) = ——— * ——— * ——— = ——— * (-R3) * Gop(s)
                vx    iR0    v-   R0
    '''
    Tn = -R3 * Gop / R0
    print('Tn =', Tn)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
    gain, phase, omega = ctrl.bode(Tn, 
                               ctrl.logspace(1, 7), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    
    # According to reciprocity relationship, loop gain T can be found.
    T = G0 * Tn / Ginf
    print('T =', T)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
    gain, phase, omega = ctrl.bode(T, 
                               ctrl.logspace(1, 7), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    
    # Finally, the closed-loop transfer function G = vout/vin can be found.
    G = Ginf * T / (1 + T) + G0 / (1 + T)
    print('G =', G)
    
    # Get the Gain and Phase and angular frequency 10Hz ~ 10MHz.
    gain, phase, omega = ctrl.bode(G, 
                               ctrl.logspace(1, 7), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    
"""
Design Application CCM Buck

let us consider the design of a combined PID compensator for the dc–dc buck 
converter system

"""
def Buck_Converter():
    # Design requirement.
    fc = 5000                         # Hz, cross frequency.
    phase_l = 52                      # Degree, phase lead at cross frequency.
    Theta = phase_l * np.pi / 180     # Convert to radian for trigonometric cal.
    
    # Circuit parameters.
    Vg = 28                           # volt, DC input voltage.
    V = 15                            # volt, DC output voltage.
    I = 5                             # amp, DC output current.
    Vref = 5                          # volt, reference voltage.
    Vm = 4                            # volt, PWM modulator peak voltage.
    L = 0.00005                       # H, Buck inductor.
    C = 0.0005                        # F, Bulk capacitor.
    # fs = 100000                     # Hz, Switching frequency.
    
    R = V / I                         # ohm, DC output load.
    H = Vref / V                      # Sensor gain.
    D = V / Vg                        # The quiescent duty cycle.
    Vc = D * Vm                       # The quiescent control voltage.
    print('SensorGain =', H, '\nDutyCycle =', D, '\nControlVoltage =', Vc)
    
    Gd0 = V / D                       # DC gain.
    omega_0 = 1 / math.sqrt(L * C)    # Corner frequency in rad/sec.
    f0 = omega_0 / (2 * np.pi)        # Corner frequency in Hz.
    Q0 = R * math.sqrt(C/L)           # Q-factor.
    print('DC_gain =', Gd0, 
          '\nCornerFreq =', f0, 'Hz', 
          '\nQ-factor =', Q0, '\n')
    
    # Create a variable 's' to allow algebra operations for SISO systems
    s = ctrl.tf('s')
    
    # The openloop control-to-output transfer function.
    Gvd = Gd0 / (1 + s/(Q0 * omega_0) + (s/omega_0)**2)
    print('Gvd =', Gvd)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(Gvd, 
                               ctrl.logspace(0, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    '''
    # The open-loop line-to-output transfer function.
    Gvg = D / (1 + s/(Q0 * omega_0) + (s/omega_0)**2)
    print(Gvg)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(Gvg, 
                               ctrl.logspace(0, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    
    # The open-loop output impedance of the buck converter
    Zout = Parallel_Resonant(R, L, C)
    '''
    # The uncompensated DC loop gain
    Tu0 = H*V/(D*Vm)
    
    # At required cross frequency 5kHz, the uncompensated loop gain
    Tu5k = Tu0 * ((f0/fc)**2)
    
    fl = fc / 10                           # Hz, Inverted zero.
    fz = resf.Zero_Frequency(fc, Theta)    # Hz
    fp = resf.Pole_Frequency(fc, Theta)    # Hz
    fp2 = fc * 10                          # Hz, high frequency pole.
    
    # To meet the cross frequency requirment
    Gc0 = math.sqrt(fz/fp) * ((fc/f0)**2) / Tu0
    print('Zero at', fz, 'Hz',
          '\nPole at', fp, 'Hz',
          '\nUncompensated DC gain =', Tu0, 
          '\nLoop gain at 5kHz =', Tu5k,
          '\nCompanstor DC gain =', Gc0, '\n')
    
    Gc = Gc0 * resf.PID_Compensator(fl, fz, fp, fp2)
    # The loop gain of the system
    Ts = Gc * Gvd * H/Vm
    print('Ts =', Ts)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(Ts, 
                               ctrl.logspace(0, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True,
                               linestyle = '-.')
    
    # The close loop ref-to-output gain of the system
    gain, phase, omega = ctrl.bode(Ts / (H * (1 + Ts)), 
                               ctrl.logspace(0, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True,
                               linestyle = ':')
    
"""
Design Application CCM Current Mode Boost

let us consider the design of a combined PID compensator for the CCM current 
mode boost converter system.

Vern's OneNote Power Electronics --> Boost Control Model

ig   2V                  sRC/2 + 1
-- = -- * ---------------------------------------
d'   R    (LCs^2 + s(L/R + Rl*C) + Rl/R + (D')^2)
                  _____________
      Vg      1  / Vg^2     Rl
D' = --- +/- -- / ----- - 4---
     2V      2 V   V^2      R
"""
def Boost_iLoop_Converter():
    # Design requirement.
    fc = 10000                        # Hz, cross frequency.
    phase_l = 52                      # Degree, phase lead at cross frequency.
    Theta = phase_l * np.pi / 180     # Convert to radian for trigonometric cal.
    
    # Circuit parameters.
    Vg = 36                           # volt, DC input voltage.
    V = 120                           # volt, DC output voltage.
    Po = 750                          # watt, DC output power.
    Vref = 3.3                        # volt, reference voltage.
    L = 0.00006                       # H, Buck inductor.
    Rl = 0.05                         # ohm, Inductor series resistance.
    C = 0.0001                        # F, Bulk capacitor.
    fs = 100000                       # Hz, Switching frequency.
    
    R = V**2 / Po                     # ohm, DC output load.
    
    # Digitalize ADC gain
    ADC_bits = 15                     # Length of ADC buffer register.
    ADC_scale = 2**ADC_bits           # Maximum digital value of analog input.
    ADC_gain = ADC_scale / Vref       # ADC value per volt.
    # ADC_resolution = 1 / ADC_gain     # ADC resolution volt per bit.
    
    # ADC_Max_V = V * 1.2               # 120% of regulation voltage for overshot.
    # V_div = 3.3 / ADC_Max_V           # Resistor dividor rate for sampling.
    # Hv = V_div * ADC_gain             # Output voltage sensor gain.
    
    Rs = 0.01                         # ohm, Boost current sensing resistor.
    I_gain = Rs * 10                  # Sensing resistor multiply Op-amp gain.
    Hi = I_gain * ADC_gain            # Boost current sensor gain.
    
    PWM_freq = 100*(10**6)            # PWM base frequency 100MHz.
    PWM_period = PWM_freq / fs        # PWM counting number for one control period.
    
    # The quiescent duty cycle, rD = 1 - D.
    rD = Vg / (2 * V) + math.sqrt(Vg**2 / V**2 - 4*Rl/R) / 2
 
    print('SensorGain =', Hi, '\nDutyCycle =', 1 - rD, '\nPWM Count =', PWM_period)
    
    # Create a variable 's' to allow algebra operations for SISO systems
    s = ctrl.tf('s')
    
    Gid = (2*V/R) * (s*R*C/2 + 1) / (L*C*s**2 + s*(L/R + Rl*C) + Rl/R + rD**2)
    print('Gid =', Gid, '\nGid0 =', (2*V/R)/(Rl/R + rD**2))
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(Gid, 
                               ctrl.logspace(0, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True)
    
    # The uncompensated DC loop gain
    Tu0 = ((2*V/R)/(Rl/R + rD**2)) * (Hi/PWM_period)
    
    fl = fc / 10                           # Hz, Inverted zero.
    fz = resf.Zero_Frequency(fc, Theta)    # Hz
    fp = resf.Pole_Frequency(fc, Theta)    # Hz
    fp2 = fc * 10                          # Hz, high frequency pole.
    
    # f0 can be got from bode plot where gain drop to Tu0 with 20dB/dec.
    f0 = 2000
    
    # To meet the cross frequency requirment
    Gc0 = math.sqrt(fz/fp) * (fc/f0) / Tu0
    
    print('Zero at', fz, 'Hz',
          '\nPole at', fp, 'Hz',
          '\nUncompensated DC gain =', 20*np.log10(Tu0), 'dB',
          '\nCompanstor DC gain =', 20*np.log10(Gc0), 'dB\n')
    
    Gc = Gc0 * resf.PID_Compensator(fl, fz, fp, fp2)
    Ts = Gc * Gid * (Hi/PWM_period)
    print('Ts =', Ts)
    
    # Get the Gain and Phase and angular frequency 1Hz ~ 1MHz.
    gain, phase, omega = ctrl.bode(Ts, 
                               ctrl.logspace(0, 6), 
                               dB = True, 
                               Hz = True,
                               deg = True, 
                               plot = True,
                               linestyle = '-.')
    
# OpAmp_Compensator()
# Buck_Converter()
Boost_iLoop_Converter()


