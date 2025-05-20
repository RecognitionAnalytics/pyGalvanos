##################################################
# CVsim.m - Cyclic voltammetry simulation
# Peter Attia
# Based on Bard and Faulkner, Appendix B
# EC mechanism
# Updated September 20, 2020
##################################################

import matplotlib.pyplot as plt
import glob
import numpy as np
import os
from lmfit import create_params, minimize
import random


#Known variables
D      = .726E-5   # [=] cm^2/s, O & R diffusion coefficient. Default = 1E-5
etai   = +0.6   # [=] V, initial overpotential. Default = +0.4
etaf   = -0.3   # [=] V, final overpotential. Default = 0
v      = .1     # [=] V/s, sweep rate. Default = 1E-3
C      = 1e-3   # [=] mol/L, initial concentration of O. Default = 1.0
T      = 298.15 # [=] K, temperature. Default = 298.15
n      = 1.0    # [=] number of electrons transfered. Default = 1
alpha  = 0.51    # [=] dimensionless charge-transfer coefficient. Default = 0.5
k0     = np.power(10,-1.3)   # [=] cm/s, electrochemical rate constant. Default = 1E-2
kc     = 1E-3   # [=] 1/s, chemical rate constant. Default = 1E-3

## Fitting VARIABLES ##
area      =  .5*1 # [=] cm^2 surface area of electrode    
OCP       = .195  # [=] V vs Ag/AgCl open circuit potential (E0)
nCarbons  = 2     # [=] number of carbons in alkane film
pinhole   = 0     # [=] fraction of surface exposed by pinholes

## PHYSICAL CONSTANTS ##
F      = 96485   # [=] C/mol, Faraday's constant
R      = 8.3145  # [=] J/mol-K, ideal gas constant
f      = F/(R*T) # [=] 1/V, normalized Faraday's constant at room temperature

## SIMULATION VARIABLES ##
L      = 500    # [=] number of iterations per t_k (pg 790). Default = 500
DM     = 0.45   # [=] model diffusion coefficient (pg 788). Default = 0.45

etai   -= OCP   # [=] V, initial overpotential (relative to redox potential). Default = +0.2
etaf   -= OCP   # [=] V, final overpotential (relative to redox potential). Default = -0.2
C = C / 1000    # Convert C from mol/L to mol/cm3

def SimulateCurve(gap,area,OCP, pinholes =0):
  ## DERIVED CONSTANTS ##
  tk  = 2*(etai-etaf)/v             # [=] s, characteristic exp. time (pg 790). In this case, total time of fwd and rev scans
  Dt  = tk/L                        # [=] s, delta time (Eqn B.1.10, pg 790)
  Dx  = np.sqrt(D*Dt/DM)            # [=] cm, delta x (Eqn B.1.13, pg 791)
  j   = int(np.ceil(4.2*L**0.5)+5)  # number of boxes for diffusion (pg 792-793). If L~200, j=65

  ## REVERSIBILITY PARAMETERS ##
  ktk    = kc*tk               # dimensionless kinetic parameter (Eqn B.3.7, pg 797)
  km     = ktk/L               # normalized dimensionless kinetic parameter (see bottom of pg 797)
  #Lambda = k0/(D*f*v)**0.5     # dimensionless reversibility parameter (Eqn 6.4.4, pg. 236-239)
  
  k = np.linspace( 0,L,L)                # time index vector
  t = Dt * k             # time vector
  eta1 = etai - v*t      # overpotential vector, negative scan
  eta2 = etaf + v*t      # overpotential vector, positive scan
  eta = np.concatenate([eta1[eta1>etaf] ,eta2[eta2<=etai],eta1[eta1>etaf] ,eta2[eta2<=etai]]) # overpotential scan, both directions
  Enorm = eta*f          # normalized overpotential
  kf = k0*np.exp(  -alpha *n*Enorm)*gap # [=] cm/s, fwd rate constant (pg 799)
  kb = k0*np.exp((1-alpha)*n*Enorm)*gap# [=] cm/s, rev rate constant (pg 799)
  
  kfclean = k0*np.exp(  -alpha *n*Enorm) # [=] cm/s, fwd rate constant (pg 799)
  kbclean = k0*np.exp((1-alpha)*n*Enorm)# [=] cm/s, rev rate constant (pg 799)

  O = C*np.ones((2*L,j)) # [=] mol/cm^3, concentration of O
  R = np.zeros((2*L,j) ) # [=] mol/cm^3, concentration of R
  JO = np.zeros((2*L)) # [=] mol/cm^2-s, flux of O at the surface

  CO = np.zeros((2*L)) # [=] mol/cm^2-s, c of O at the surface
  CR = np.zeros((2*L)) # [=] mol/cm^2-s, c of O at the surface

  for i_time in range(0,2*L-1):
      # Update bulk concentrations of O and R
      for i_dist in range( 1,int(j)-1):
          O[i_time+1,i_dist] = O[i_time,i_dist] + DM*(O[i_time,i_dist+1]+O[i_time,i_dist-1]-2*O[i_time,i_dist])
          R[i_time+1,i_dist] = R[i_time,i_dist] + DM*(R[i_time,i_dist+1]+R[i_time,i_dist-1]-2*R[i_time,i_dist]) - km * R[i_time,i_dist]


      # flux assuming that the surface can be split into two areas, one with the film and one without
      # the total flux from the surface adjacent solution is then the sum of the fluxes from each area
      # this is unlikely to model a surface with a lot of different species on the surface, which should
      #give a general broading of the peak
      cleanO = O[i_time+1,1]*pinholes
      cleanR = R[i_time+1,1]*pinholes
      dirtyO= O[i_time+1,1]*(1-pinholes)
      dirtyR=R[i_time+1,1]*(1-pinholes)
      dirtyJ= ( kf[i_time+1]*dirtyO - kb[i_time+1]*dirtyR ) / (1 + Dx/D*(kf[i_time+1] + kb[i_time+1]) )
      cleanJ= ( kfclean[i_time+1]*cleanO - kbclean[i_time+1]*cleanR ) / (1 + Dx/D*(kfclean[i_time+1] + kbclean[i_time+1]) )
      JO[i_time+1] = dirtyJ+cleanJ
      
      #flux equation without the pinhole approximation
      #JO[i_time+1]   = ( kf[i_time+1]*O[i_time+1,1] - kb[i_time+1]*R[i_time+1,1] ) / (1 + Dx/D*(kf[i_time+1] + kb[i_time+1]) )

      # Update surface concentrations
      O[i_time+1,0] = O[i_time+1,1] - JO[i_time+1]*Dx/D
      CO[i_time+1]  = R[i_time+1,1]
      R[i_time+1,0] = R[i_time+1,1] + JO[i_time+1]*Dx/D - km*R[i_time+1,0]
      CR[i_time+1]  = O[i_time+1,1]


  # Calculate current density, Z, from flux of O
  Z = -n*F*JO * 1000 # [=] A/cm^2 -> mA/cm^2, current density

  ## PLOT RESULTS ##
  # Sometimes length(eta) = length(Z) + 1. If this is the case, truncate last value
  if len(eta) > len(Z):
      eta = eta[1:-1]
  return eta[int(len(Z)/2):]+OCP,Z[int(len(Z)/2):]*area


eta,Z=SimulateCurve(np.exp(-.9*nCarbons),area,OCP,pinholes=pinhole)

maxV=eta[ np.argmax(Z)]
minV=eta[ np.argmin(Z)]
print((maxV-minV)*1000,"mV")

plt.plot(eta  ,Z,label='Simulation' )

 
plt.xlabel('Overpotential (V)')
plt.ylabel('Current density (mA)')
plt.legend()
plt.show()