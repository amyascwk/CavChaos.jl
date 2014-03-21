"""
List of refractive index distributions
"""

import numpy as np
import matplotlib.pyplot as plt

#Generic functions for speed =================================================
def sq(x): return x*x
def cb(x): return x*x*x


#Refractive index functions ==================================================
#List of indices (returns [n, dn/dr, dn/dtheta])
def n_uniform(r,theta,n0=1,D=True):
    """Refractive index of uniform cavity at radial position r and angle theta"""
    n = n0*np.ones(np.shape(r))
    if D:
        dn = np.zeros(np.shape(r))
        return n, dn, dn
    else:
        return n

def n_gaussianxy(r,theta,n0=2.,A=-0.2,sigx=0.5,sigy=0.3,x0=0.,y0=0.,D=True):
    """Refractive index of gaussian variation at radial position r and angle theta
        Parameters: n0=2, A=-0.2, sigx=0.5,sigy=0.3,x0=0.,y0=0.
    """
    x,y = r*np.cos(theta),r*np.sin(theta)
    n = n0 + A*np.exp(-0.5*(sq((x-x0)/sigx) +sq((y-y0)/sigy)))
    if D:
        #dn_r = \frac{\partial n}{\partial r},
        dn_r = (n-n0)/(np.sqrt(sq(x)+sq(y)))*(x*(x0-x)/sq(sigx)+y*(y0-y)/sq(sigy))
        #dn_theta = \frac{\partial n}{\partial \theta}}
        dn_theta = (n-n0)*(y*(x-x0)/sq(sigx)-x*(y-y0)/sq(sigy))
        return n, dn_r, dn_theta
    else:
        return n

def n_core(r,theta,n0=1.,A=0.4,sigma=0.5,D=True):
    """Refractive index of symmetric gaussian variation at radial position r and angle theta
        Parameters: n0=1., A=0.4, sigma=0.5
    """
    n = n0 + A*np.exp(-0.5*sq(r/sigma))
    if D:
        dn_r = -(n-n0)*(r/sq(sigma))
        dn_theta = np.zeros(np.shape(n))
        return n, dn_r, dn_theta
    else:
        return n

def n_test(r,theta,D=True):
    n = -0.8*r+1
    if D:
        dn_r = -0.8
        dn_theta = 0
        return n, dn_r, dn_theta
    else:
        return n

