#! /usr/bin/env python
# -*- coding:utf-8 -*-

# ====================================================================
# This script reads Hybrid mean BL fields and extract some quantities
# It works only in python 2 -> python2 plot_extract_bl_profiles_all.py
# Author: P. S. Volpiani
# Date: 07/07/2022
# ====================================================================

import os, sys
import matplotlib.mathtext as mathtext
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import operator

from math import *
from pylab import *
from numpy import *
from matplotlib import *
from matplotlib import rc

c = 0

cases = [
         "../BL-M2/m228_twtr100_bl/bl_1536x384x128-1-242files",
         "../BL-M2/m228_twtr050_bl/bl1_1800x480x256-1-240files",
         "../BL-M2/m228_twtr190_bl/case2_1536x288x128-5-600files",
         "../BL-M2/m228_twtr100_bl_hRe/bl_1800x480x256-1-241files",
         "../BL-M5/TwTr08/bl_twtr08-2-152files",
         "../BL-M5/TwTr08/bl_twtr08-5-135files",
         "../BL-M5/TwTr08/bl_twtr08-6-102files",
         "../BL-M5/TwTr19/bl_twtr19-5-142files",
         ]

for case in cases:
    
    print case
    
    file = case+".stats"
    
    if    (case=="../BL-M5/TwTr19/bl_twtr19-5-142files"):  # 1.9 - Very high Re
      muref=0.0000120 ; muexp=0.75; TW=10.355; Tr=TW/1.9; xplot=40.; xsp=65; dz=5./140.;
    elif (case=="../BL-M5/TwTr08/bl_twtr08-2-152files"):    # 0.8 - High Re
      muref=0.0000180 ; muexp=0.75; TW=4.36; Tr=TW/0.8; xplot=40.; xsp=65; dz=5./200.;
    elif  (case=="../BL-M5/TwTr08/bl_twtr08-5-135files"):  # 0.8 - Very high Re test
      muref=0.0000090 ; muexp=0.75; TW=4.36; Tr=TW/0.8; xplot=40.; xsp=65; dz=5./340.;
    elif  (case=="../BL-M5/TwTr08/bl_twtr08-6-102files"):  # 0.8 - Between high and very high Re test
      muref=0.0000120 ; muexp=0.75; TW=4.36; Tr=TW/0.8; xplot=40.; xsp=65; dz=5./300.;
    elif  (case=="../BL-M5/TwTr19_Mu-025/bl_twtr19-3-107files"):  # 1.9 - Very high Re nu-0.25
      muref=0.000050 ; muexp=-0.25; TW=10.355; Tr=TW/1.9; xplot=45.; xsp=55; dz=5./200.;
    elif  (case=="../BL-M5/TwTr19_Mu+025/bl_twtr19-2-71files"):  # 1.9 - Very high Re nu+0.25
      muref=0.000030 ; muexp=0.25; TW=10.355; Tr=TW/1.9; xplot=50.; xsp=56; dz=5./160.;
    elif  (case=="../BL-M2/m228_twtr100_bl/bl_1536x384x128-1-242files"):    # Adi
      muref=0.000122 ; muexp=0.75; TW=1.920; Tr=TW/1.0; xplot=35.; xsp=60; dz=5./128;
    elif  (case=="../BL-M2/m228_twtr190_bl/case2_1536x288x128-5-600files"): # Hot
      muref=0.000120 ; muexp=0.75; TW=3.658; Tr=TW/1.9; xplot=35.; xsp=60; dz=5./128;
    elif(case=="../BL-M2/m228_twtr050_bl/bl1_1800x480x256-1-240files"):     # Cold
      muref=0.000129 ; muexp=0.75; TW=0.963; Tr=TW/0.5; xplot=35.; xsp=60; dz=5./256;
    elif  (case=="../BL-M2/m228_twtr100_bl_hRe/bl_1800x480x256-1-241files"):# Adi - High Reynolds
      muref=0.000060 ; muexp=0.75; TW=1.920; Tr=TW/1.0; xplot=35.; xsp=60; dz=5./256;

    with open(file, 'rb') as f:
    
        #  Read until 3 consequtive %%%
        lastThree = f.read(3)
        while ( lastThree != ['%','%','%'] ):
            next = f.read(1); lastThree = [ lastThree[1], lastThree[2], next[0] ]
        # Read new line
        next = f.read(1)
        # Read integers nx, ny, nv
        arg = np.fromfile(f, count=3, dtype='int32')
        nx = arg[0]; ny = arg[1]; nv = arg[2];
        # Reed coordinates
        X   = np.fromfile(f, count=nx, dtype='float64')
        Y   = np.fromfile(f, count=ny, dtype='float64')
        # Read quantities
        avg = np.fromfile(f, count=nx*ny*nv, dtype='float64')
        avg = avg.reshape([nx, ny, nv], order='F');


    # Load variables
    r   = np.zeros((nx, ny));   r[:,:] = avg[:,:,0];
    p   = np.zeros((nx, ny));   p[:,:] = avg[:,:,4];
    u   = np.zeros((nx, ny));   u[:,:] = avg[:,:,10];
    v   = np.zeros((nx, ny));   v[:,:] = avg[:,:,11];
    w   = np.zeros((nx, ny));   w[:,:] = avg[:,:,12];
    T   = np.zeros((nx, ny));   T[:,:] = avg[:,:,13];
    rr  = np.zeros((nx, ny));  rr[:,:] = avg[:,:,15];
    pp  = np.zeros((nx, ny));  pp[:,:] = avg[:,:,16];
    uu  = np.zeros((nx, ny));  uu[:,:] = avg[:,:,21];
    vv  = np.zeros((nx, ny));  vv[:,:] = avg[:,:,22];
    ww  = np.zeros((nx, ny));  ww[:,:] = avg[:,:,23];
    TT  = np.zeros((nx, ny));  TT[:,:] = avg[:,:,24];
    uv  = np.zeros((nx, ny));  uv[:,:] = avg[:,:,25];
    uT  = np.zeros((nx, ny));  uT[:,:] = avg[:,:,28];
    vT  = np.zeros((nx, ny));  vT[:,:] = avg[:,:,29];
    S12 = np.zeros((nx, ny)); S12[:,:] = avg[:,:,35];


    # Compute/assume global quantities
    gamma = 1.4
    Pr = 0.7
    gasR = mean(mean(p[:,:]/r[:,:]/T[:,:]))
    cp = gamma*gasR/(gamma-1.)


    # Compute wall quantities
    Tw = np.zeros(nx)
    if (TW>0.001): Tw[:] = TW*np.ones(nx)
    else: Tw[:] = T[:,0]

    # Dynamic viscosity at the wall
    muw  = np.zeros(nx);  muw[:] = muref * Tw[:]**muexp ;
    # Dynamic viscosity at first cell
    mu0  = np.zeros(nx);  mu0[:] = muref * T[:,0]**muexp ;
    # Mean dynamic viscosity
    muef = np.zeros(nx); muef[:] = ( muw[:] + mu0[:] )/2. ;
    # Wall shear stress
    tauw = np.zeros(nx); tauw[:] = muef[:] * u[:,0] / Y[0] ;
    # Heat transfer at the wall
    qw   = np.zeros(nx);   qw[:] = cp * muef[:] / Pr  * ( T[:,0] - Tw[:]) / Y[0] ;
    # Density at the wall
    rhow = np.zeros(nx); rhow[:] = p[:,0] / Tw[:] / gasR ;
    # Friction velocity
    utau = np.zeros(nx); utau[:] = sqrt( abs(tauw[:]) / rhow[:] ) ;
    # Viscous length scale
    lv   = np.zeros(nx);   lv[:] = muw[:] / rhow[:] / utau[:] ;


    # Find reference location
    min_index, min_value = min(enumerate(abs(X[:]-xplot)), key=operator.itemgetter(1)); iplot = min_index;


    # Compute integral quantities
    delta  = np.zeros(nx);
    dstar  = np.zeros(nx);
    theta  = np.zeros(nx);
    uvd    = np.zeros((nx, ny));
    ystar  = np.zeros((nx, ny));
    ustar  = np.zeros((nx, ny));
    mu     = np.zeros(ny);
    
    for I in range(0, nx):


        # Find edge of boundary layer, erring on the side of a bit too much outside
        f = np.zeros(ny); f[:]=r[I,:] * ww[I,:]/abs(tauw[I]);
        max_index, max_value = max(enumerate(f), key=operator.itemgetter(1)); J = max_index;
        while (J<ny) and (f[J]>0.02): J=J+1;
        U  = np.zeros(ny);  U[:] = u[I,:] / u[I,J] ;
        RU = np.zeros(ny); RU[:] = r[I,:] * u[I,:] / r[I,J] / u[I,J] ;


        # Find delta99(I)
        J=0;
        while ( U[J]<0.99 ): J=J+1;
        J=J-1;
        delta[I] = Y[J] + ( Y[J+1] - Y[J] ) / ( U[J+1] - U[J] ) * ( 0.99 - U[J] ) ;


        # Find delta*(I)
        f[:] = 1. - RU[:] ;  fw = 1. ;
        dstar[I] = (fw+f[0])/2. * Y[0] ;
        J=0;
        while (Y[J]<delta[I]):
            J=J+1;
            dstar[I] = dstar[I] + mean(f[J-1:J]) * (Y[J]-Y[J-1]) ;


        # Find deltaTheta
        f[:] = RU[:] * (1.-U[:]) ;  fw = 0 ;
        theta[I] = (fw+f[0])/2. * Y[0] ;
        J=0;
        while (Y[J]<delta[I]):
            J=J+1;
            theta[I] = theta[I] + mean(f[J-1:J]) * (Y[J]-Y[J-1]) ;


        # Find Van Driest transformed velocities
        f[:] = 1/utau[I] * sqrt(r[I,:]/rhow[I]) ;  fw = 1/utau[I] ;
        uvd[I,0] = ( fw + f[0] )/2. * ( u[I,0] - 0. );
        for J in range(1, ny):
            uvd[I,J] = uvd[I,J-1] + mean(f[J-1:J]) * (u[I,J]-u[I,J-1]) ;


        # Find Trettel & Larsson transformed velocities
        Te=1;
        mu[:] = muref * (T[I,:]/Te)**muexp ;
        ystar[I,:] = sqrt( r[I,:]*abs(tauw[I]) ) / mu[:] * Y[:] ;
        f = mu/abs(tauw[I]) ;  fw = muw[I]/abs(tauw[I]) ;
        ustar[I,0] = (fw+f[0])/2. * (ystar[I,0]-0.)/(Y[0]-0.) * u[I,0] ;
        for J in range(1, ny):
            dystardy = ( ystar[I,J] - ystar[I,J-1] ) / ( Y[J] - Y[J-1] ) ;
            ustar[I,J] = ustar[I,J-1] + mean(f[J-1:J]) * dystardy * (u[I,J]-u[I,J-1]) ;



    # Reynolds numbers
    Ue=1; rhoedge=1;
    Re_delta  = np.zeros(nx);  Re_delta[:]  = rhoedge * Ue * delta[:] / muref ;
    Re_theta  = np.zeros(nx);  Re_theta[:]  = rhoedge * Ue * theta[:] / muref ;
    Re_delta2 = np.zeros(nx);  Re_delta2[:] = rhoedge * Ue * theta[:] / muw[:] ;
    Re_tau    = np.zeros(nx);  Re_tau[:]    = delta[:] / lv[:] ;
    Re_star   = np.zeros(nx);  Re_star[:]   = delta[:] * sqrt(rhoedge * abs(tauw[:])) / muref ;
    

    # Compute dx+, dy+, dz+ at xplot
    dx  = ( X[iplot+1] - X[iplot] ) / lv[iplot]
    dy1 = ( Y[0] - 0. ) / lv[iplot]
    dyN = ( Y[-1]-Y[-2] ) / lv[iplot]
    dz  = dz / lv[iplot]

    # Print on screen
    print "nx        = ", arg[0]
    print "ny        = ", arg[1]
    print "delta99   = ", delta[iplot]
    print "delta*    = ", dstar[iplot]
    print "Re_delta  = ", Re_delta[iplot]
    print "Re_theta  = ", Re_theta[iplot]
    print "Re_delta2 = ", Re_delta2[iplot]
    print "Re_tau    = ", Re_tau[iplot]
    print "Re_star   = ", Re_star[iplot]
    print "dx+       = ", dx
    print "dy1+      = ", dy1
    print "dyN+      = ", dyN
    print "dz+       = ", dz


    # Write data file
    
    name = case.replace('/', "-")
    name = name.replace('..', "")
    print (name)
    # Write tecplot file
    data_file = open("./file"+name+".txt","w")

    if ('BL-M2' in case):
      M_inf = 2.28
      data_file.write('# Reference : "Volpiani, P. S., Bernardini, M., & Larsson, J. (2018). \n')
      data_file.write('#              Effects of a nonadiabatic wall on supersonic shock/boundary-layer interactions. \n')
      data_file.write('#              Physical Review Fluids, 3(8), 083401."\n')
    else:
      M_inf = 5.0
      data_file.write('# Reference : "Volpiani, P. S., Bernardini, M., & Larsson, J. (2020). \n')
      data_file.write('#              Effects of a nonadiabatic wall on hypersonic shock/boundary-layer interactions. \n')
      data_file.write('#              Physical Review Fluids, 5(1), 014602."\n')

    data_file.write('TITLE="Mean profiles for compressible boundary layer"\n\n')
    
    data_file.write('M_inf     =  {:>20}'.format( M_inf )              )
    data_file.write( '\n' )
    data_file.write('Tw/Tr     =  {:>20}'.format( TW/Tr )              )
    data_file.write( '\n' )
    data_file.write('Cf        =  {:>20}'.format( 2.*tauw[iplot] )     )
    data_file.write( '\n' )
    data_file.write('Re_delta  =  {:>20}'.format( Re_delta[iplot] )    )
    data_file.write( '\n' )
    data_file.write('Re_theta  =  {:>20}'.format( Re_theta[iplot] )    )
    data_file.write( '\n' )
    data_file.write('Re_delta2 =  {:>20}'.format( Re_delta2[iplot] )   )
    data_file.write( '\n' )
    data_file.write('Re_tau    =  {:>20}'.format( Re_tau[iplot] )      )
    data_file.write( '\n' )
    data_file.write('Re_star   =  {:>20}'.format( Re_star[iplot] )     )
    data_file.write( '\n' )
    data_file.write('dx+       =  {:>20}'.format( dx )                 )
    data_file.write( '\n' )
    data_file.write('dyw+      =  {:>20}'.format( dy1*2. )             )
    data_file.write( '\n' )
    data_file.write('dz+       =  {:>20}'.format( dz )                 )
    data_file.write( '\n\n' )
    data_file.write('VARIABLES = "Y/DELTA" , "Y+" , "U+" , "R11/TAU_W" , "R22/TAU_W" , "R33/TAU_W" , "R12/TAU_W" , "RHO/RHO_W" , "T/T_W", "M/M_W" , "Prms/P_INF" \n')
    for j in xrange(len(Y)):
        data_file.write('{:>25}'.format( str(Y[j]/delta[iplot])      ) )
        data_file.write('{:>25}'.format( str(Y[j]/lv[iplot])         ) )
        data_file.write('{:>25}'.format( str(u[iplot,j]/utau[iplot]) ) )
        data_file.write('{:>25}'.format( str(r[iplot,j]*uu[iplot,j]/tauw[iplot]) ) )
        data_file.write('{:>25}'.format( str(r[iplot,j]*vv[iplot,j]/tauw[iplot]) ) )
        data_file.write('{:>25}'.format( str(r[iplot,j]*ww[iplot,j]/tauw[iplot]) ) )
        data_file.write('{:>25}'.format( str(r[iplot,j]*uv[iplot,j]/tauw[iplot]) ) )
        data_file.write('{:>25}'.format( str(r[iplot,j]/rhow[iplot]) ) )
        data_file.write('{:>25}'.format( str(T[iplot,j]/Tw[iplot])   ) )
        data_file.write('{:>25}'.format( str(muref*(T[iplot,j]/Te)**muexp/muw[iplot]) ) )
        #data_file.write('{:>25}'.format( str(sqrt(pp[iplot,j])/gasR) ) )
        data_file.write( '\n' )
    data_file.close()


    c = c + 1
