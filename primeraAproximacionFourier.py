#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 09:13:05 2023

@author: alberto
"""

from numpy import trapz, linspace, pi, sin, zeros, cos, round, fft


def dameLosCoeficientes(f, thetaArray, period, nCoefs):

    # Como quiero nCoefs defino a tal efecto el vector dnd los voy a guardar
    # Notese que lo he definido como vector complejo pq mis coeficietes lo son!
    coefVector = zeros(nCoefs, dtype = 'complex128') 
    for n in range(0, nCoefs):
        cosTerm = cos(n*thetaArray); sinTerm = sin(n*thetaArray)
        a_n = 2*trapz(f*cosTerm, thetaArray)/period; b_n = 2*trapz(f*sinTerm, thetaArray)/period;
        coefVector[n] = (a_n - 1j*b_n)/2
        print('El coeficiente del armonico %s es %s' % (n, round(coefVector[n], 3)))
    
    return coefVector


def main():

    nCoefs = 8; period = 2*pi
    thetaArray = linspace(0, period, 100)
    f = sin(thetaArray)**3
    #f = cos(thetaArray)**2
    coefVector = dameLosCoeficientes(f, thetaArray, period, nCoefs)
    
    
    

if __name__ == '__main__':
	main()


