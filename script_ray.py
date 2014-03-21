#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Ray tracing for dielectic cavity
"""

__date__ = "Date: 2013/05/28"
__author__= "Author: Michael Pasek, Yidong Chong"

import numpy as npy
from matplotlib import rc
import matplotlib.pyplot as plt
import scipy.integrate
import scipy.optimize
from mpl_toolkits.axes_grid.inset_locator import inset_axes, zoomed_inset_axes, mark_inset
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from scipy.integrate import ode
import sys

rc('text', usetex=True)

def circular_cavity_boundary(theta):
    '''Return radial distance of a circle at polar angle theta.'''
    return 1.0

def limacon_cavity_boundary(theta):
    '''Return radial distance of a limacon at polar angle theta.'''
    return 1.0 + 0.2 * npy.cos(theta)

def ray_integrate(r0, tet0, pphi0, t1, dt, index_fun, boundary_fun):
    '''Integrate the ray by solving the Eikonal equation.
    r0, tet0, pphi0:  initial values of r, theta, and k_phi.
    t1:               time to integrate to (starting from t=0).
    dt:               time step for the integrator.
    index_fun:        Function which accepts polar coordinates (r,p)
                      and returns an array [n, drn, dtn].
    boundary_fun:     Function which accepts one argument (theta)
                      and returns the boundary R(theta).'''

    # Set up the initial coordinate vector y, which consists of
    # [r, \phi, k_r, k_\phi].  The initial value of k_r is obtained
    # from eikonal equation, (grad S)^2 = n^2 = p_r^2 + p_\phi^2 / r^2
    # (Hakan p.41)
    n0 = index_fun(r0, tet0)
    prsq = n0[0]*n0[0] - pphi0 * pphi0 / (r0 * r0)
    assert prsq > 0
    y0 = [r0, tet0, npy.sqrt(prsq), pphi0]

    # For details on the Runge-Kutta integrator, see
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
    # or E. Hairer et al.,
    #    Solving Ordinary Differential Equations i. Nonstiff Problems.
    #    2nd edition. Springer Series in Computational Mathematics,
    #    Springer-Verlag (1993)

    # Derivative function for the initial value problem,
    #     \frac{d \vec{y}}{d s} = f(s, \vec{y})
    # See (2.15) and (2.16) in Hakan's thesis.
    def f(s, y):
        n = index_fun(y[0], y[1])
        return [y[2]/n[0], y[3]/(y[0]*y[0]*n[0]), n[1]+(1./n[0])*y[3]*y[3]/y[0]**3, n[2]]

    r = ode(f, jac=None).set_integrator('dopri5', atol=1e-12, rtol=1e-12)
    r.set_initial_value(y0, 0.0)

    pss_vec = npy.empty(0)
    rvec, tvec = [], []
    r0, p0 = y0[0], y0[1]
    R0 = boundary_fun(p0)
    while r.successful() and r.t < t1:
        r.integrate(r.t + dt)
        R1 = boundary_fun(r.y[1])
        if r.y[0] > R1:
            ## The ray has passed the boundary of the cavity.
            assert r.y[2] > 0.0
            ynew, pss = ray_bounce(r0, p0, r.y, R0, R1, index_fun)
            r.y[:] = ynew
            pss_vec = npy.append(pss_vec, pss)
        R0 = R1
        r0, p0 = r.y[0], r.y[1]
        rvec.append(r0)
        tvec.append(p0)

        # Print to verify conservation of Hamiltonian function
        H = r.y[2]*r.y[2]+r.y[3]*r.y[3]/(r.y[0]*r.y[0]) - index_fun(r.y[0],r.y[1])[0]**2
        sys.stdout.write("\r%9.3f , H = %10.6f" % (r.t, H))
        sys.stdout.flush()
    sys.stdout.write('\n')
    return pss_vec, npy.array(rvec), npy.array(tvec)

def ray_bounce(r0, p0, pos, R0, R1, index_fun):
    '''Re-compute a bounced ray coordinate vector.
    r0, p0:    Position coordinates at previous time step
    pos:       Ray coordinate vector, which has overshot the
               boundary at theta=pos[1].
    R0, R1:    r-values of boundary at p0 and pos[1].
    index_fun: Function returning optical index array.'''

    r1, p1, kr, kp = pos[0], pos[1], pos[2], pos[3]

    c0 = npy.cos(p0)
    x0, X0 = r0 * c0, R0 * c0
    s0 = npy.sin(p0)
    y0, Y0 = r0 * s0, R0 * s0
    c1 = npy.cos(p1)
    x1, X1 = r1 * c1, R1 * c1
    s1 = npy.sin(p1)
    y1, Y1 = r1 * s1, R1 * s1

    ## Intersection of trajectory with cavity boundary
    u, v = (y1-y0)/(x1-x0), (Y1-Y0)/(X1-X0)
    x = (y0 - Y0 + X0 * v - x0 * u) / (v - u)
    y = y0 + (x-x0) * u
    r, p = npy.sqrt(x*x + y*y), npy.arctan2(y, x)
    phi = npy.mod(p, 2.0 * npy.pi) # For PSS

    ## u is parallel to surface, v is normal to surface.
    u = npy.array([X1-X0, Y1 - Y0])
    u = u / npy.sqrt(npy.dot(u, u))
    v = npy.array([u[1], -u[0]])
    if v[0] * npy.cos(p) + v[1] * npy.sin(p) < 0:
        v = -v

    d = npy.array([x1 - x, y1 - y])
    d = d / npy.sqrt(npy.dot(d, d))
    sinchi, coschi = -(d[0]*v[1] - d[1]*v[0]), d[0]*v[0] + d[1]*v[1] # For PSS

    ## kphi is conserved.  Re-compute kr from the Hamiltonian.
    n  = index_fun(r, p)[0]
    #kr = - npy.sqrt(n*n - kp*kp/ (r*r))
    vr = v[0]*npy.cos(p)+v[1]*npy.sin(p)
    vp = v[1]*npy.cos(p)-v[0]*npy.sin(p)
    dr = d[0]*npy.cos(p)+d[1]*npy.sin(p)
    dp = d[1]*npy.cos(p)-d[0]*npy.sin(p)
    #From v x (Eq.(18)), Chap.3.2 of Born&Wolf - Principles of Optics (6th ed.), setting n1=n2
    #and using sqrt(kr**2+(kp/r)**2)=n from eikonal equation
    kr=n*(dr-2.*coschi*vr)
    kp=n*r*(dp-2.*coschi*vp)

    return npy.array([r, p, kr, kp]), npy.array([phi,sinchi])

def test_index(y1,y2):
    ''' Return optical index n inside the cavity, together with its gradient components.'''
    n0 = 2.
    x0=0.0
    y0=0.0
    sigx=0.50
    sigy=0.30
    A = -0.2
    x = y1*npy.cos(y2)
    y = y1*npy.sin(y2)
    n = n0 + A*npy.exp(-((x-x0)**2/(2.*sigx**2)+(y-y0)**2/(2.*sigy**2)))
    # drn is \frac{\partial n}{\partial r}, and dtn is \frac{\partial n}{\partial phi}}
    drn = (n-n0)/(npy.sqrt(x**2+y**2))*(x*(x0-x)/sigx**2+y*(y0-y)/sigy**2)
    dtn = (n-n0)*(y*(x-x0)/sigx**2-x*(y-y0)/sigy**2)
    return [n, drn, dtn]

def plot_ray(rvec, tvec, boundary_fun):
    '''Plot the ray and the cavity boundary.'''
    rays_plot = plt.figure(1, [8,8]).add_subplot(111, aspect='equal')
    rays_plot.plot(rvec * npy.cos(tvec), rvec * npy.sin(tvec), 'b-')

    nt = 1000
    tbound, rbound = npy.linspace(0.,2.*npy.pi,nt), npy.zeros(nt)
    for ii in range(nt):
        rbound[ii] = boundary_fun(tbound[ii])
    rays_plot.plot(rbound * npy.cos(tbound), rbound * npy.sin(tbound),'k-')

def plot_poincare(pss_vec, fontsize=18):
    '''Plot the Poincare surface of section.'''
    poinc_plot = plt.figure(2, [8, 7]).add_subplot(111)
    pss_vec = npy.reshape(pss_vec,(-1,2))
    phip, sin_chip = pss_vec[:,0], pss_vec[:,1]
    poinc_plot.plot(phip,sin_chip,'o')

    plt.xlabel(r"$\displaystyle \phi$", fontsize=fontsize)
    plt.ylabel(r"$\displaystyle \sin \chi$", fontsize=fontsize)
    plt.xticks([0, npy.pi/2, npy.pi, 3*npy.pi/2, 2*npy.pi],[r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

def plot_index(boundary_fun,index_fun):
    '''Generate a colored contour plot of the cavity's refractive index.'''
    index_plot = plt.figure(3, [8, 8]).add_subplot(111, aspect='equal')
    nt= 100
    nr = 80
    tbound, rbound = npy.linspace(0.,2.*npy.pi,nt), npy.zeros(nt)
    for ii in range(nt):
        rbound[ii] = boundary_fun(tbound[ii])
    index_plot.plot(rbound * npy.cos(tbound), rbound * npy.sin(tbound),'k-')
    tn, rn = npy.linspace(0.,2.*npy.pi,nt), npy.zeros((nt,nr))
    for ii in range(nt):
        rn[ii,:] = npy.linspace(0., boundary_fun(tn[ii]), nr)
    tn2 = npy.array([tn,]*nr)
    XN = rn*npy.cos(npy.transpose(tn2))
    YN = rn*npy.sin(npy.transpose(tn2))
    index_plot.contourf(XN, YN, index_fun(rn,npy.transpose(tn2))[0], 25, cmap=cm.Blues)

def prog():
    r0 = 0.82
    tet0 = 0.2
    pphi0 = 0.5
    t1 = 200
    dt = 1e-3
    boundary = limacon_cavity_boundary
    index = test_index

    pss_vec, rvec, tvec = ray_integrate(r0, tet0, pphi0, t1, dt, index, boundary)
    plot_ray(rvec, tvec, boundary)
    plot_poincare(pss_vec)
    plot_index(boundary,index)

    plt.show()

if __name__ == '__main__':
    prog()
    print ("Hasta la vista, baby.")
