"""
List of cavities
"""

import numpy as np
import matplotlib.pyplot as plt

#Generic functions for speed =================================================
def sq(x): return x*x
def cb(x): return x*x*x

#Boundary functions ==========================================================
#List of boundaries (returns [r, phi_normal])
def circbnd(theta,par=None,shownormal=False):
    '''Radial position of circular boundary at angle theta'''
    assert type(theta) != list
    r = np.ones(np.shape(theta))
    if shownormal: return [r, theta]
    else: return r

def ellipbnd(theta,e=0.8,shownormal=False):
    '''Radial position of elliptical boundary of eccentricity e at angle theta'''
    assert type(theta) != list
    r = 1./np.sqrt(sq(np.cos(theta)) + sq(np.sin(theta))/(1-sq(e)))
    if shownormal:
        dr_theta = -0.5*np.sin(2*theta)*sq(e)*cb(r)/(1-sq(e))
        return [r, theta-np.arctan2(dr_theta,r)]
    else: return r

def polygbnd(theta,n=3,rounded=0,shownormal=False):
    '''Radial position of polygonal boundary of n sides at angle theta'''
    assert type(theta) != list
    angdiv = 2*np.pi/n
    thetap = np.mod(theta,angdiv)-0.5*angdiv
    r = (1-rounded)/np.cos(thetap) + rounded/np.cos(0.5*angdiv)
    if shownormal:
        if type(theta) != np.ndarray:
            theta = np.array([theta])
            thetap = np.array([thetap])
            r = np.array([r])
        cnrmiss = np.mod(theta,angdiv).nonzero()[0]
        normang = theta
        dr_theta = (1-rounded)*np.sin(thetap[cnrmiss])/sq(np.cos(thetap[cnrmiss]))
        normang[cnrmiss] -= np.arctan2(dr_theta,r[cnrmiss])
        return [r, normang]
    else: return r

def limaconbnd(theta,eps=0.2,shownormal=False):
    '''Radial position of limacon boundary with quadrupolar deformations of amplitude eps at angle theta'''
    assert type(theta) != list
    r = 1. + eps*np.cos(2*theta)
    if shownormal:
        dr_theta = -2*eps*np.sin(2*theta)
        return [r, theta-np.arctan2(dr_theta,r)]
    else: return r

def genbnd(theta,epsA=[0,0.05],epsB=[0],shownormal=False):
    '''Radial position of generalized boundary with deformations of amplitudes epsilon_i for order i,
    where r = 1 + \sum_{i=1}^m{\epsilon_{A,i}\cos(i\theta)} + \sum_{i=1}^n{\epsilon_{B,i}\sin(i\theta)}'''
    assert type(theta) != list
    orderA = np.arange(1,len(epsA)+1)[:,np.newaxis]
    orderB = np.arange(1,len(epsB)+1)[:,np.newaxis]
    epsA = np.array(epsA)[np.newaxis,:]
    epsB = np.array(epsB)[np.newaxis,:]
    qA = np.dot(orderA,theta[np.newaxis,:])
    qB = np.dot(orderB,theta[np.newaxis,:])
    r = 1. + np.dot(epsA,np.cos(qA)) + np.dot(epsB,np.sin(qB))
    if shownormal:
        dr_theta = -np.dot(epsA*orderA.T,np.sin(qA)) + np.dot(epsB*orderB.T,np.cos(qB))
        return [r, theta-np.arctan2(dr_theta,r)]
    else: return r

def arcbnd(theta,R=0.89,eps=0.28,eps1=0.06,d=0.06,shownormal=False):
    '''Radial position of Asymmetric Resonant Cavity boundary in PRL 105,103902'''
    assert type(theta) != list
    r = R*(1+eps*np.cos(theta))*(1-eps1*np.cos(2*theta)) + d
    if shownormal:
        dr_theta = R*((0.5*eps1-1)*eps*np.sin(theta)+2*eps1*np.sin(2*theta)+1.5*eps*eps1*np.sin(3*theta))
        return [r, theta-np.arctan2(dr_theta,r)]
    else: return r

