# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 08:34:33 2022

Current-Programmed Control, including:
    > Peak current mode(PCM) control.
    > Average current mode (ACM) control

Vern's OneNote Power Electronics --> Current-Programmed Control

                 Outputs: duty cycle d
                               o
        _______________________|_______________________
       |                       d                      |        
       |                                              |   Parameters:
       |                      CPM                     |   Rf, fs, L, Va
       |                                              |
       | Control     Current         1            2   |
       ————————————————————————————————————————————————
           |            |            |            |
           o            o            o            o
Inputs: <vc(t)>Ts   Rf<iL(t)>Ts   <v1(t)>Ts   <v2(t)>Ts
Where:
    <vc(t)>Ts = Rf * <ic(t)>Ts, Average control input.
    Rf<iL(t)>Ts: Sensed average inductor current.
    <v1(t)>Ts:   Average voltage appplied across the inductor during 1st interval.
    <v2(t)>Ts:   Average voltage appplied across the inductor during 2nd interval.
    Rf:          Equivalnt current sense reistance.
    fs:          Switching frequency.
    L:           Inductance.
    Va = ma * Ts * Rf, Artificial ramp amplitude.

@author: Vern Chen

"""

