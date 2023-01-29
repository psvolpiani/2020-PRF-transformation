#! /usr/bin/env python
# -*- coding:utf-8 -*-

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

rc('text', usetex=True)
rc('font', family='serif')
rc('lines', linewidth=1.5)
rc("font", size=16)
plt.rc('legend',**{'fontsize':14})
matplotlib.rcParams.update({'axes.labelsize': 20})

colormap  = 'jet' #'RdBu'
linestyle = ["--","-","--","-.", "--"]
cor       = ["r", "k", "b", "g", "m"]

# Numerical data from Duan et al.
fileXP  = "../Duan/M14Tw018_Stat.dat.txt"

U_inf   = 1882.2
rho_inf = 0.017
T_inf   = 47.4
T_sut   = 273.15
T_w     = 300.0
mu_inf  = 1.458e-6*(T_inf)**(3./2.)/(T_inf+110.4)
mu_w    = 1.458e-6*(T_w)**(3./2.)/(T_w+110.4)
u_tau   = 67.6
lv      = 0.1024 /1000.

Y   = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(0,), delimiter=" ")
U   = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(4,), delimiter=" ")
P   = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(5,), delimiter=" ")
T   = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(6,), delimiter=" ")
MU  = 1.458e-6*(T*T_inf)**(3./2.)/((T*T_inf)+110.4)
R   = P/T
tauw = u_tau**2*R[0]*rho_inf

# Plus quantities
yp  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(2,), delimiter=" ")
up  = U*U_inf/u_tau
R_Rw  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(7,), delimiter=" ")
T_Tw  = T/T[0]
M_Mw  = MU/MU[0]
N_Nw  = M_Mw/R_Rw

# Van Driest quantities
yd  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(2,), delimiter=" ")
ud  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(24,), delimiter=" ")

# Trettel & Larsson quantities
yt  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(3,), delimiter=" ")
ut  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(25,), delimiter=" ")

# Reynolds stresses
R11  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(8,), delimiter=" ")
R22  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(9,), delimiter=" ")
R33  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(10,), delimiter=" ")
R12  = numpy.loadtxt(fileXP,comments='%',skiprows=0, usecols=(15,), delimiter=" ")

# Find delta99
J=0;
while ( U[J]<0.99 ): J=J+1;
J=J-1;
delta = Y[J] + ( Y[J+1] - Y[J] ) / ( U[J+1] - U[J] ) * ( 0.99 - U[J] ) ;
print (delta)

# Find delta*
f = np.zeros(len(Y));
RU = R*U
f[:] = 1. - RU[:] ;
dstar = 0. ;
J=0;
while (Y[J]<delta):
  J=J+1;
  dstar = dstar + mean(f[J-1:J]) * (Y[J]-Y[J-1]) ;
print (dstar)

# Find deltaTheta
f[:] = RU[:] * (1.-U[:]) ;  fw = 0 ;
theta = (fw+f[0])/2. * Y[0] ;
J=0;
while (Y[J]<delta):
  J=J+1;
  theta = theta + mean(f[J-1:J]) * (Y[J]-Y[J-1]) ;
print (theta)


# Transformation kernels
ny = len(Y)
f  = np.zeros(ny);
g  = np.zeros(ny);

# Find New Volpiani et al. transformed velocities
# ===============================================
yn = np.zeros(ny);
un = np.zeros(ny);

alpha=1.5
beta =0.5

f = (R_Rw/M_Mw)**beta/M_Mw**(alpha-beta)
yn[0] = 0.
for J in range(1, ny):
  yn[J] = yn[J-1] + mean(f[J-1:J]) * (yp[J]-yp[J-1]) ;

g = (R_Rw/M_Mw)**beta/M_Mw**(alpha-beta) * M_Mw
un[0] = 0.
for J in range(1, ny):
  un[J] = un[J-1] + mean(g[J-1:J]) * (up[J]-up[J-1]);

min_index, min_value = min(enumerate(abs(Y[:]-delta)), key=operator.itemgetter(1));

# Reynolds numbers
# ================
Re_delta  = rho_inf * U_inf * delta / mu_inf ;
Re_theta  = rho_inf * U_inf * theta / mu_inf ;
Re_delta2 = rho_inf * U_inf * theta / mu_w ;
Re_tau    = delta / lv ;
Re_star   = delta * sqrt(rho_inf * abs(tauw)) / mu_inf ;
Re_n      = yn[min_index] ;

# Print on screen
print ("Re_delta  = ", Re_delta)
print ("Re_theta  = ", Re_theta)
print ("Re_delta2 = ", Re_delta2)
print ("Re_tau    = ", Re_tau)
print ("Re_star   = ", Re_star)
print ("Re_n      = ", Re_n)

# ======== #
# Figure u #
# ======== #

fig, axes = plt.subplots(figsize=(6, 5))
y_min=0.; y_max=30;
t = np.arange(1, 10000, 1)
tlin = np.arange(0.5, 11, 0.1)
tlog = np.zeros(size(t));
tlog[:] = 5.2+2.44*log(t[:]);

#axes.scatter(yp4[markers4], up4[markers4], color='gray', marker='o', facecolors='none', linewidth='1')
axes.semilogx(yd,ud, linestyle = '-', color = 'k');
axes.semilogx(yt,ut, linestyle = '--', color = 'g');
axes.semilogx(yn,un, linestyle = '-.', color = 'r');

axes.semilogx(tlin, tlin, linestyle = ':', color = 'k')
axes.semilogx(t[10:-1], tlog[10:-1], linestyle = ':', color = 'k')

axes.set_xlabel( r"$y^+$");
axes.set_ylabel(r"$u^+$");
axes.set_xlim([0.5, 10000])
axes.set_ylim([y_min, y_max])
fig.subplots_adjust(left=0.2)
fig.subplots_adjust(bottom=0.18)
grid(b=True, which='major', color='gray', linestyle=':')
plt.minorticks_on()
fig_name = "./bl_u_M14_Duan_NEW.pdf"
#fig.savefig(fig_name, dpi=300)
plt.show()

