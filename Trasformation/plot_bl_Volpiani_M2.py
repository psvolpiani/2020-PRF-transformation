#! /usr/bin/env python
# -*- coding:utf-8 -*-

# ====================================================================
# This script reads Hybrid stats fields and plot BL quantities
# It works only in python 2 -> python2 plot_bl_Volpiani_M2.py
# Author: P. S. Volpiani
# Date: 07/07/2022
# ====================================================================

import os, sys
import matplotlib as mpl
mpl.use('tkagg') # for Mac
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

rc('text', usetex=True)
rc('font', family='serif')
rc('lines', linewidth=2)
rc("font", size=16)
plt.rc('legend',**{'fontsize':14})
matplotlib.rcParams.update({'axes.labelsize': 20})

colormap  = 'jet' #'RdBu'
linestyle = ["-","-.","--","-.", "--"]
cor       = ["k", "C3", "C0", "g", "m"]
print_figure = "yes"
c = 0

# First file
cases = [
         "../BL-M2/m228_twtr100_bl/bl_1536x384x128-1-242files",
#         "../BL-M2/m228_twtr190_bl/case2_1536x288x128-5-600files",
#         "../BL-M2/m228_twtr050_bl/bl1_1800x480x256-1-240files",
        ]

fileXP = "../EXP/Schlatter2010_Retheta1410.txt"
y2     = numpy.loadtxt(fileXP,skiprows=14, usecols=(0,))
yp2    = numpy.loadtxt(fileXP,skiprows=14, usecols=(1,))
up2    = numpy.loadtxt(fileXP,skiprows=14, usecols=(2,))
urmsp2 = numpy.loadtxt(fileXP,skiprows=14, usecols=(3,))
vrmsp2 = numpy.loadtxt(fileXP,skiprows=14, usecols=(4,))
wrmsp2 = numpy.loadtxt(fileXP,skiprows=14, usecols=(5,))
uvp2   = numpy.loadtxt(fileXP,skiprows=14, usecols=(6,))
markers2 = np.arange(1, len(yp2), 5)

for case in cases:
    
    print case
    
    file = case+".stats"
    
    if (case=="../BL-M5/TwTr08/bl_twtr08-2-152files"):    # 0.8 - High Re
      muref=0.0000180 ; muexp=0.75; TW=4.36; Tr=TW/0.8; xplot=40.; xsp=65; dz=5./200.;
    elif  (case=="../BL-M5/TwTr08/bl_twtr08-5-135files"):  # 0.8 - Very high Re test
      muref=0.0000090 ; muexp=0.75; TW=4.36; Tr=TW/0.8; xplot=40.; xsp=65; dz=5./340.;
    elif  (case=="../BL-M5/TwTr08/bl_twtr08-6-102files"):  # 0.8 - Between high and very high Re test
      muref=0.0000120 ; muexp=0.75; TW=4.36; Tr=TW/0.8; xplot=40.; xsp=65; dz=5./300.;
    elif  (case=="../BL-M5/TwTr19/bl_twtr19-5-142files"):  # 1.9 - Very high Re
      muref=0.0000120 ; muexp=0.75; TW=10.355; Tr=TW/1.9; xplot=40.; xsp=65; dz=5./140.;
    elif  (case=="../BL-M2/m228_twtr100_bl/bl_1536x384x128-1-242files"):    # Adi
      muref=0.000122 ; muexp=0.75; TW=1.920; Tr=TW/1.0; xplot=48.; xsp=60; dz=5./128;
    elif  (case=="../BL-M2/m228_twtr190_bl/case2_1536x288x128-5-600files"): # Hot
      muref=0.000120 ; muexp=0.75; TW=3.658; Tr=TW/1.9; xplot=53.; xsp=60; dz=5./128;
    elif(case=="../BL-M2/m228_twtr050_bl/bl1_1800x480x256-1-240files"):     # Cold
      muref=0.000129 ; muexp=0.75; TW=0.963; Tr=TW/0.5; xplot=55.; xsp=60; dz=5./256;
    elif  (case=="../BL-M2/m228_twtr100_bl_hRe/bl_1800x480x256-1-241files"):# Adi - High Reynolds
      muref=0.000060 ; muexp=0.75; TW=1.920; xplot=35.; xsp=60; dz=5./256;


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
    lplus  = np.zeros((nx, ny));
    ushe   = np.zeros((nx, ny));
    uc     = np.zeros((nx, ny));
    yc     = np.zeros((nx, ny));
    uj     = np.zeros((nx, ny));
    yj     = np.zeros((nx, ny));
    
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

        # Apply She et al. (2017) formula
        l0=1.03; yplus_sub=9.7; kappa=0.45; yplus_buf=41.;
        yplus   = Y[:]/lv[I]
        lplus[I,:] = l0*(yplus/yplus_sub)**(3./2.) * ( 1 + (yplus/yplus_sub)**(4.) )**(1./8.) * ( 1 + (yplus/yplus_buf)**(4.) )**(-1./4.)
        f[:] = (-1.+sqrt(1.+4.*lplus[I,:]**2.) )/ (2.*lplus[I,:]**2.) ;  fw = 0. ;
        ushe[I,0] = 0.
        for J in range(1, ny):
            ushe[I,J] = ushe[I,J-1] + mean(f[J-1:J]) * ( yplus[J] - yplus[J-1] )
              
        # Find Brun et al. transformed velocities
        mu[:] = muref * (T[I,:]/Te)**muexp ;
        R = np.zeros(ny); R = r[I,:] /rhow[I]
        N = np.zeros(ny); N = mu / muw[I] / R
        f = 1./(R*N) ;  fw = 1.;
        yc[I,0] = ( fw + f[0] )/2. * ( Y[0] - 0. );
        for J in range(1, ny):
          yc[I,J] = yc[I,J-1] + mean(f[J-1:J]) * (Y[J]-Y[J-1]) ;
        f = 1./( N * sqrt(R)) * Y[:]/yc[I,:] ;  fw = 0. ;
        uc[I,0] = ( fw + f[0] )/2. * ( u[I,0] - 0. )/utau[I] ;
        for J in range(1, ny):
            uc[I,J] = uc[I,J-1] + mean(f[J-1:J]) * (u[I,J]-u[I,J-1])/utau[I] ;
              
        # Find NEW transformation
        alpha=1.5
        beta =0.5
        M = np.zeros(ny); M = muref * (T[I,:]/Te)**muexp / muw[I] ;
        R = np.zeros(ny); R = r[I,:] /rhow[I]
        N = np.zeros(ny); N = M / R

        f = (R/M)**beta/M**(alpha-beta) ;  fw = 1.;
        yj[I,0] = ( fw + f[0] )/2. * ( Y[0] - 0. );
        for J in range(1, ny):
          yj[I,J] = yj[I,J-1] + mean(f[J-1:J]) * (Y[J]-Y[J-1]) ;

        f = (R/M)**beta/M**(alpha-beta) * M ;  fw = 1.;
        uj[I,0] = ( fw + f[0] )/2. * ( u[I,0] - 0. )/utau[I] ;
        for J in range(1, ny):
            uj[I,J] = uj[I,J-1] + mean(f[J-1:J]) * (u[I,J]-u[I,J-1])/utau[I] ;
              
    # Reynolds numbers
    Ue=1; rhoedge=1;
    Re_delta  = np.zeros(nx);  Re_delta[:]  = rhoedge * Ue * delta[:] / muref ;
    Re_theta  = np.zeros(nx);  Re_theta[:]  = rhoedge * Ue * theta[:] / muref ;
    Re_delta2 = np.zeros(nx);  Re_delta2[:] = rhoedge * Ue * theta[:] / muw[:] ;
    Re_tau    = np.zeros(nx);  Re_tau[:]    = delta[:] / lv[:] ;
    Re_star   = np.zeros(nx);  Re_star[:]   = delta[:] * sqrt(rhoedge * abs(tauw[:])) / muref ;
    
    Re_b   = np.zeros(nx);
    Re_j   = np.zeros(nx);
    for I in range(1, nx):
      min_index, min_value = min(enumerate(abs(Y[:]-delta[I])), key=operator.itemgetter(1));
      Re_b[I]      = yc[I,min_index] / lv[I] ;
      Re_j[I]      = yj[I,min_index] / lv[I] ;

    # Compute dx+, dy+, dz+ at xplot
    dx  = ( X[iplot+1] - X[iplot] ) / lv[iplot]
    dy1 = ( Y[0] - 0. ) / lv[iplot]
    dyN = ( Y[-1]-Y[-2] ) / lv[iplot]
    dz  = dz / lv[iplot]

    # Print on screen
    print "iplot     = ", iplot
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
    print "Re_b      = ", Re_b[iplot]
    print "Re_j      = ", Re_j[iplot]
    

    # ======= #
    # Figures #
    # ======= #

    if (c==0): fig2, axes2 = plt.subplots(figsize=(5.5, 4))
    y_min=0.; y_max=25;
    t = np.arange(1, 10000, 1)
    tlin = np.arange(0.5, 11, 0.1)
    tlog = np.zeros(size(t));
    tlog[:] = 5.2+(1./0.41)*log(t[:]);

    axes2.semilogx(yj[iplot,:]/lv[iplot],uj[iplot,:], linestyle = linestyle[c], color = cor[c]);
    axes2.semilogx(tlin, tlin, 'k:')
    axes2.semilogx(t[10:-1], tlog[10:-1], 'k:')
    axes2.scatter(yp2[markers2], up2[markers2], color='gray', marker='^', facecolors='white', linewidth='1')

    axes2.set_xlabel( r"$y^+_{V}$");
    axes2.set_ylabel(r"$u^+_{V}$");
    axes2.set_xlim([0.5, 1000])
    axes2.set_ylim([y_min, y_max])
    fig2.subplots_adjust(left=0.2)
    fig2.subplots_adjust(bottom=0.18)
    plt.minorticks_on()
    #fig2.savefig("fig_uv_M2_NEW.pdf", dpi=300)

    if (c==0): fig3, axes3 = plt.subplots(figsize=(5.5, 4))
    x_min=0.; x_max=1.2;

    axes3.scatter(yp2[markers2], urmsp2[markers2]**2, color='gray', marker='^', facecolors='white', linewidth='1')
    axes3.scatter(yp2[markers2], vrmsp2[markers2]**2, color='gray', marker='^', facecolors='white', linewidth='1')
    axes3.scatter(yp2[markers2], wrmsp2[markers2]**2, color='gray', marker='^', facecolors='white', linewidth='1')
    axes3.scatter(yp2[markers2], - uvp2[markers2]**2, color='gray', marker='^', facecolors='white', linewidth='1')

    axes3.semilogx(yj[iplot,:]/lv[iplot],r[iplot,:]*uu[iplot,:]/tauw[iplot], linestyle = '-' , color = cor[c]);
    axes3.semilogx(yj[iplot,:]/lv[iplot],r[iplot,:]*vv[iplot,:]/tauw[iplot], linestyle = '--', color = cor[c]);
    axes3.semilogx(yj[iplot,:]/lv[iplot],r[iplot,:]*ww[iplot,:]/tauw[iplot], linestyle = '-.', color = cor[c]);
    axes3.semilogx(yj[iplot,:]/lv[iplot],r[iplot,:]*uv[iplot,:]/tauw[iplot], linestyle = ':' , color = cor[c]);

    axes3.set_xlabel( r"$y^+_{V}$");
    axes3.set_ylabel(r"$\bar{\rho}\widetilde{u''_i u''_j} / \tau_w$");
    axes3.set_xlim([0.5, 1000])
    fig3.subplots_adjust(left=0.2)
    fig3.subplots_adjust(bottom=0.18)
    plt.minorticks_on()
    #fig3.savefig("fig_Rij_M2_NEW.pdf", dpi=300)


    c = c + 1


plt.show()

