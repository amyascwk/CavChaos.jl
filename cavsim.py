#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cavity simulator and mode finder, version 0.0.2
"""

__date__ = "Date: 2013/07/24"
__author__= "Author: Michael Pasek, Yidong Chong, Amyas Chew"

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import ode
import sys, os, time, re

from cavitybnds import *
from cavityidxs import *

#Generic functions for speed =================================================
def sq(x): return x*x
def cb(x): return x*x*x


#Boundary functions ==========================================================
#Add boundary
def addbnd(bndname):
    accept = False
    while not accept:
        newbnd = ''
        print('\nDefine the boundary function below.')
        line = 'def '+bndname+'(theta'
        line += raw_input('\n>>>'+line)
        while line != '':
            newbnd += '\n'+line
            line = raw_input('...')
        try:
            exec(newbnd)
            bndplot = plt.figure('Test boundary',figsize=(5,5))\
                      .gca(aspect='equal',title=bndname)
            plotbnd(bndplot,bndfunc=eval(bndname))
            plt.show()
            accept = (raw_input('\nAccept boundary? (y/n): ')+'n')[0].lower()=='y'
        except(SyntaxError):
            print('\nFunction is not valid, boundary not added. Try again.\n')

    #accepted
    cavs = open('cavitybnds.py','a')
    cavs.write(newbnd+'\n\n')
    cavs.close()

def listbnds():
    cavs = open('cavitybnds.py','r')

def showbnd(bndname):
    cavs = open('cavitybnds.py','r')

#Plot boundaries
def plotbnd(ploth,bndfunc=circbnd):
    """Plots the cavity boundary of function bndfunc"""
    theta = np.linspace(0.,2*np.pi,500)
    r = bndfunc(theta)
    ploth.plot(r*np.cos(theta),r*np.sin(theta),'g')
    ploth.axis([-r[0],r[0],-r[0],r[0]])


#Refractive index functions ==================================================
def addidx(nname):
    accept = False
    while not accept:
        newidx = ''
        print('\nDefine the refractive index distribution function below.')
        line = 'def '+nname+'(r,theta'
        line += raw_input('\n>>>'+line)
        while line != '':
            newidx += '\n'+line
            line = raw_input('...')
        try:
            exec(newidx)
            idxplot = plt.figure('Test refractive index',figsize=(5,5))\
                      .gca(aspect='equal',title=nname)
            x,y = np.meshgrid(np.linspace(-1,1,100),np.linspace(-1,1,100))
            z = eval(nname+'(np.sqrt(sq(x)+sq(y)),np.arctan2(y,x))')
            idxplot.contourf(x,y,z,100,cmap='bone_r')
            contobj = idxplot.contour(x,y,z,100,cmap='bone_r',alpha=0.2)
            levels = np.array(contobj.levels)[0:-1:10]
            plt.clabel(contobj,levels,colors='k')
            
            plt.show()
            accept = (raw_input('\nAccept refractive index? (y/n): ')+'n')[0].lower()=='y'
        except(SyntaxError):
            print('\nFunction is not valid, refractive index not added. Try again.\n')

    #accepted
    cavs = open('cavityidxs.py','a')
    cavs.write(newidx+'\n\n')
    cavs.close()


#Plot index
def plotidx(ploth,bndfunc,nfunc):
    '''Plots a contour plot of the cavity interior's refractive index of function nfunc.'''
    thetares,rres = 100,100
    theta = np.linspace(0.,2*np.pi,thetares)
    rbnd = bndfunc(theta)[np.newaxis,:]
    interv = np.linspace(0.,1.,rres)[:,np.newaxis]
    r = np.dot(interv,rbnd)
    theta2 = np.array([theta]*rres)
    x,y,z = r*np.cos(theta2),r*np.sin(theta2),nfunc(r,theta2,D=False)
    
    ploth.contourf(x,y,z,100,cmap='bone_r')
    contobj = ploth.contour(x,y,z,100,cmap='bone_r',alpha=0.2)
    levels = np.array(contobj.levels)[0:-1:10]
    plt.clabel(contobj,levels,colors='k')


#Ray tracing =================================================================
#Compute p coordinates from direction
def phi2p(r,theta,phi,nfunc):
    #pr,ptheta obtained from eikonal eqn (See Hakan p.41)
    #(\del S)^2 = n^2 = p_r^2 + \frac{p_\theta^2}{r^2}
    n = nfunc(r,theta,D=False)
    pr,ptheta = n*np.cos(phi-theta),n*r*np.sin(phi-theta)
    return pr,ptheta

#Ray reflection
def raybounce(r0,theta0,S_y,bndfunc,nfunc):
    '''Re-compute a bounced ray coordinate vector.
    r0, theta0:     Pre-collision position
    S_y:            Post-collision coordinate vector
    bndfunc:        Cavity boundary function
    nfunc:          Refractive index function'''

    #Setup coordinate values
    r,theta,pr,ptheta = S_y
    R0 = bndfunc(theta0)
    R = bndfunc(theta)
    
    c0,s0 = np.cos(theta0),np.sin(theta0)
    x0,y0 = r0*c0, r0*s0
    X0,Y0 = R0*c0, R0*s0
    c,s = np.cos(theta),np.sin(theta)
    x,y = r*c, r*s
    X,Y = R*c, R*s

    #Find intersection of trajectory with cavity boundary
    l = ((Y-Y0)*(x0-X0)-(X-X0)*(y0-Y0))/((X-X0)*(y-y0)-(Y-Y0)*(x-x0))
    xc = x0 + l*(x-x0)
    yc = y0 + l*(y-y0)
    rc, thetac = np.sqrt(sq(xc)+sq(yc)), np.mod(np.arctan2(yc,xc),2*np.pi)

    #Compute reflected ray
    Rc,normal = bndfunc(thetac,shownormal=True)
    if np.abs(1-rc/Rc) > 2e-7:
        print('Warning: Error in radial position is '+str(1-rc/Rc)+'\n')
    chi = np.arctan2(y-y0,x-x0) - normal
    phi = np.pi - chi + normal
    prc,pthetac = phi2p(rc,thetac,phi,nfunc)
    
    return np.array([rc,thetac,prc,pthetac]), np.array([[thetac,np.sin(chi)]])

#Ray traversal
def raytrav(r0=0,theta0=0,phi0=0,bndfunc=circbnd,nfunc=n_uniform,dt=0.001,tmax=10):
    """
    Numerically integrate the ray by solving the Eikonal equation
    r0, theta0, ptheta0:    Initial values of r, theta, and ptheta.
    tmax:                   Total time for integration.
    dt:                     Time step for the integration.
    nfunc:                  Function which accepts polar coordinates (r,theta)
                              and returns an array [n, drn, dtn].
    bndfunc:                Function which accepts one argument (theta)
                              and returns the boundary R(theta).
    """
    
    #Initiation
    pr0,ptheta0 = phi2p(r0,theta0,phi0,nfunc)
    y0 = [r0, theta0, pr0, ptheta0]
    
    def odefunc(t,y):
        """
        The derivative function f for the traversal of the ray, such that
            y' = f(y)
        where y = [r,theta,pr,ptheta] is 
        expresses the parameters of the ray at a given point.
        (See (2.15) and (2.16) in Hakan's thesis)
            r' = \frac{p_r}{n}, \theta' = \frac{p_\theta}{r^2 n},
            p_r' = \frac{\partial n}{\partial r} + \frac{p_\theta^2}{n r^3},
            p_\theta' = \frac{\partial n}{\partial \theta}
        """
        n = nfunc(y[0],y[1])
        return [y[2]/n[0], y[3]/(n[0]*sq(y[0])), n[1]+sq(y[3])/(n[0]*cb(y[0])), n[2]]
    
    
    #Integration
    S = ode(odefunc).set_integrator('dopri5', atol=1e-12, rtol=1e-12)
    S.set_initial_value(y0)
    Nmax = np.floor(tmax/dt)+2
    rval, thetaval = np.zeros(Nmax), np.zeros(Nmax)
    rval[0], thetaval[0] = r0, theta0
    pssval = np.empty((0,2))
    idx = 0
    while S.successful() and idx < Nmax-1:
        S.integrate(S.t + dt)
        R = bndfunc(S.y[1])
        if S.y[0] > R:
            #Collision with boundary
            S.y[:], pss = raybounce(rval[idx], thetaval[idx],S.y,bndfunc,nfunc)
            pssval = np.append(pssval,pss,0)
        #Record values
        idx += 1
        rval[idx] = S.y[0]
        thetaval[idx] = S.y[1]

        #Verify conservation of Hamiltonian function
        if np.mod(idx,1000) == 0:
            H = S.y[2]*S.y[2]+S.y[3]*S.y[3]/(S.y[0]*S.y[0]) - sq(nfunc(S.y[0],S.y[1])[0])
            if np.abs(H) > 1e-9:
                print('Warning: Hamiltonian is %0.8f, might not be conserved\n'%(H))
            sys.stdout.write('\rt = %0.3f'%S.t)
            sys.stdout.flush()
    return pssval, rval, thetaval

#Plot rays
def plotrays(ploth,rval,thetaval):
    '''Plot the rays'''
    ploth.plot(rval*np.cos(thetaval),rval*np.sin(thetaval),'y-',aa=False)


#Poincare surface of section =================================================
def plotpss(ploth,pssvalval,bndfunc=None,nfunc=None):
    '''Plot the Poincare surface of section.'''
    mpl.rc('text', usetex=True)

    #Emphasize first set of points
    pssval = np.array(pssvalval[0])
    ploth.plot(pssval[:,0],pssval[:,1],'.',ms=8,mec='k',mfc='r',mew=1.5,picker=True,zorder=len(pssvalval))

    #Scatter remaining points
    for i in range(1,len(pssvalval)):
        pssval = np.array(pssvalval[i])
        color = colorfunc(clusteredness(pssval))
        ploth.plot(pssval[:,0],pssval[:,1],'.',color=color,markersize=2,picker=True)
    ploth.axis([0.,2*np.pi,-1,1])
    plt.xlabel(r'$\displaystyle\theta$', fontsize=18)
    plt.ylabel(r'$\displaystyle\sin\chi$', fontsize=18)

    #Show critical angle
    if bndfunc != None and nfunc != None:
        theta = np.linspace(0.,2*np.pi,200)
        sin_chi0 = 1./nfunc(bndfunc(theta),theta,D=False)
        ploth.plot(theta,sin_chi0,'m',theta,-sin_chi0,'m')


#Configure parameters ========================================================
def writecfg(r0val,theta0val,phi0val,bndname,bndpar,nname,npar,dt,tmax):
    '''Writes the configuration file based on given parameters'''
    cfg = open('config_default.txt','w')
    cfg.write('##Cavity parameters\n'\
                 +'def bndfunc(theta,**kw): return '+bndname+'(theta,'+bndpar+',**kw)\n'\
                 +'def nfunc(r,theta,**kw): return '+nname+'(r,theta,'+npar+',**kw)\n'\
                 +'\n##Simulation parameters\n'\
                 +'dt = '+str(dt)+'\ntmax = '+str(tmax)+'\n'\
                 +'\n##Initial conditions\n'\
                 +'theta0val = '+theta0val+'\nphi0val = '+phi0val+'\nr0val = '+r0val)
    cfg.close()

def stdincfg():
    '''Specify parameters as arguments'''
    #Defaults
    bndname,bndpar = 'circbnd','*[]'
    nname,npar = 'n_uniform','*[]'
    dt,tmax = 1e-3,20

    #Stdin defined
    theta0val = sys.argv[2]
    sinchi0val = sys.argv[3]
    if len(sys.argv) > 4: bndname = sys.argv[4]
    if len(sys.argv) > 5: bndpar = sys.argv[5]
    if len(sys.argv) > 6: nname = sys.argv[6]
    if len(sys.argv) > 7: npar = sys.argv[7]
    if len(sys.argv) > 8: dt = float(sys.argv[8])
    if len(sys.argv) > 9: tmax = float(sys.argv[9])
    def bndfunc(theta,**kw):
        return eval(bndname+'(theta,*'+bndpar+',**kw)')
    r0val,normangval = bndfunc(np.array(eval(theta0val)),shownormal=True)
    r0val = str(r0val.tolist())
    phi0val = str((np.arcsin(np.array(eval(sinchi0val))) + normangval).tolist())

    #Write config file
    writecfg(r0val,theta0val,phi0val,bndname,bndpar,nname,npar,dt,tmax)
    
def manualcfg():
    '''Specify parameters interactively'''
    #Get input
    print('Cavity Parameters '+40*'='+'\n')
    bndname = raw_input('Boundary function (default is circbnd): ')
    if bndname == '': bndname = 'circbnd'
    bndpar = raw_input('Boundary parameters (default is *[]): ')
    if bndpar == '': bndpar = '*[]'

    nname = raw_input('Refractive index function(default is n_uniform): ')
    if nname == '': nname = 'n_uniform'
    npar = raw_input('Refractive index parameters(default is *[]): ')
    if npar == '': npar = '*[]'

    print('\nSimulation Parameters '+40*'='+'\n')    
    dt = raw_input('Time step (default is 1e-3): ')
    if dt == '': dt = 1e-3
    else: dt = float(dt)

    tmax = raw_input('Time of simulation (default is 20): ')
    if tmax == '': tmax = 20
    else: tmax = float(tmax)

    print('\nInitial Conditions '+40*'='+'\n')
    theta0val = raw_input('Initial theta values: ')
    phi0val = raw_input('Initial phi values: ')
    r0val = raw_input('Initial r values (default is bndfunc(np.array(theta0val))): ')
    if r0val == '': r0val = 'bndfunc(np.array(theta0val))'

    #Write config file
    writecfg(r0val,theta0val,phi0val,bndname,bndpar,nname,npar,dt,tmax)


#Runs ========================================================================
#Read config
def cfgrun(grprunname=None,cfgname='default',parallel=[1,0]):
    '''Collates data from the simulation on a set of initial conditions specified in a config file'''
    #Read config
    if not os.path.exists('config_'+cfgname+'.txt'):
        print('\nConfig file not found. Run \'python cavsim.py c\' to generate config, then edit file where necessary.')
        exit()
    cfg = open('config_'+cfgname+'.txt','r')
    lines = list(cfg)
    clbl1 = lines.index('##Cavity parameters\n')
    clbl2 = lines.index('##Initial conditions\n')
    label = ''.join(lines[clbl1:clbl2])
    cfg.seek(0)
    exec(cfg)
    cfg.close()

    #Cleanup initial conditions
    if type(r0val) == list: r0val = np.array(r0val)
    if type(theta0val) == list: theta0val = np.array(theta0val)
    if type(phi0val) == list: phi0val = np.array(phi0val)
    #Convert single entry to list
    #(lengths must match or must be scalar)
    Nruns = max(np.shape(r0val),np.shape(theta0val),np.shape(phi0val))
    if Nruns == (): Nruns = (1L)
    ones = np.ones(Nruns)
    r0val = r0val*ones
    theta0val = theta0val*ones
    phi0val = phi0val*ones

    #Select parallelised job
    r0val = r0val[parallel[1]::parallel[0]]
    theta0val = theta0val[parallel[1]::parallel[0]]
    phi0val = phi0val[parallel[1]::parallel[0]]

    #Prepare results record
    if grprunname == None:
        grprunname = str(int(time.time()))
    grprunpath = os.path.join('results','grprun_'+grprunname)
    if not os.path.exists(grprunpath): os.makedirs(grprunpath)
    resultpath = os.path.join(grprunpath,'result_'+grprunname)
    if parallel[0] > 1: resultpath += '_'+str(parallel[1])
    result = open(resultpath+'.txt','a')
    result.write(55*'#'+'\nlabel = """'+label+'"""\n'+55*'#'+'\n\n')

    #Plot cavity
    rayplot = plt.figure(figsize=(5,5)).gca(aspect='equal')
    plotbnd(rayplot,bndfunc)
    plotidx(rayplot,bndfunc,nfunc)
    plt.gcf().tight_layout()
    cavaxis = plt.axis()
    if not os.path.exists(os.path.join(grprunpath,'cavity.jpg')):
        plt.savefig(os.path.join(grprunpath,'cavity.jpg'))
    plt.close()

    #Run!!!
    for i in range(len(theta0val)):
        print('Running '+str(i+1)+' of '+str(len(theta0val))+' runs:\n')
        pssval, rval, thetaval = raytrav(r0val[i],theta0val[i],phi0val[i],bndfunc,nfunc,dt,tmax)
        rayplot = plt.figure(figsize=(5,5)).add_subplot(111,aspect='equal')
        plotrays(rayplot,rval,thetaval)
        plt.axis(cavaxis)
        plt.axis('off')
        plt.gcf().tight_layout()
        
        if not np.isnan(pssval).any():
            #Write results
            runname = str(int(np.mod(time.time(),1e5)))
            while True:
                figpath = os.path.join(grprunpath,runname+'.png')
                if not os.path.exists(figpath): break
                else: runname = str(int(runname)+1)
            open(figpath,'a').close()
            plt.savefig(figpath,transparent=True)
            result.write('##Run '+runname+' '+40*'='\
                         +'\nimgid{0}{1})\nr0val{0}{2})\ntheta0val{0}{3})\nphi0val{0}{4})\n\npssvalval{0} {5} )\n\n\n'\
                         .format('.append(',runname,r0val[i],theta0val[i],phi0val[i],pssval.tolist()))
        
        plt.close()
    result.close()
    return grprunname


#Plot results ================================================================
def onpick(ev):
    '''Interactive behaviour of results plot: emphasis, highlighting, recomputation'''
    axs = plt.gcf().get_axes()
    if len(axs) >= 3:
        #Figure is not closed
        mev = ev.mouseevent
        art = ev.artist
        if mev.button == 1:
            #Turn off previous emphasis
            for a in axs[2].findobj(mpl.lines.Line2D):
                if a.get_ms() == 8:
                    idx = int(re.subn(r'\D*(\d+)',r'\1',a.get_label())[0])
                    pssval = pssvalval[idx]
                    color = colorfunc(clusteredness(pssval))
                    a.set(mec=color,mfc=color,mew=None,ms=2,zorder=0)

            #Emphasize current artist
            art.set(mec='k',mfc='r',mew=1.5,ms=8,zorder=len(imgid))

            #Display rayplot corresponding to current artist
            emph = int(re.subn(r'\D*(\d+)',r'\1',art.get_label())[0])
            print('Run '+str(imgid[emph])+' selected.')
            rayimg = plt.imread(os.path.join('results','grprun_'+grprunname,str(imgid[emph])+'.png'))
            axs[1].get_children()[3].remove()
            init = r'($r_0$, $\theta_0$, $\phi_0$) = (%0.4f,  %1.4f,  %2.4f)'%(r0val[emph],theta0val[emph],phi0val[emph])\
                   +'\n'+r'($\theta_0$, $\sin\chi_0$) = (%0.4f,  %1.4f)'%(pssvalval[emph][0][0],pssvalval[emph][0][1])
            axs[1].set_xlabel(init)
            axs[1].imshow(rayimg)
            print('Clusteredness score: '+str(clusteredness(pssvalval[emph])))
            
        elif mev.button == 3:
            if art.get_ms() == 10:
                #remove from recompute list
                theta = art.get_xdata()[0]
                sinchi = art.get_ydata()[0]
                art.remove()
                recomp.remove([theta,sinchi])
                print('Removed ('+str(theta)+','+str(sinchi)+') from recompute list')
            elif art.get_ms() != 4:
                #Highlight current artist
                art.set(mec='k',mfc='m',ms=4,mew=1,zorder=len(imgid))
            else:
                #Unhighlight current artist
                tag = int(re.subn(r'\D*(\d+)',r'\1',art.get_label())[0])
                color = colorfunc(clusteredness(pssvalval[tag]))
                art.set(mec=color,mfc=color,ms=2,mew=None,zorder=0)

        plt.draw()
        
    else: plt.close()

def onpress(mev):
    '''Interactive behaviour of results plot: recomputation'''
    axs = plt.gcf().get_axes()
    if mev.dblclick == True and mev.button == 1:
        ax = plt.gca()
        theta = mev.xdata
        sinchi = mev.ydata
        
        #Add to recompute list and mark
        recomp.append([theta,sinchi])
        ax.plot([theta],[sinchi],'+r',ms=10,picker=2)
        plt.draw()
        print('Added ('+str(theta)+','+str(sinchi)+') to recompute list')

def plotresults(grprunname):
    '''Plot the results of simulations on a cavity'''
    #Get results
    imgid, pssvalval = [], []
    r0val,theta0val,phi0val = [],[],[]
    grprunpath = os.path.join('results','grprun_'+grprunname)
    result = open(os.path.join(grprunpath,'result_'+grprunname+'.txt'),'r')
    exec(result)
    exec(label)
    result.close()

    #Prepare Cavity Label text
    label0 = label
    label = re.subn(r'##(.+)',r'\\textbf{\1}:',label)[0] #bold section comments
    label = re.subn(r'\ntmax',', tmax',label)[0]
    label = re.subn(r'\n+','\n',label)[0] #remove empty lines
    label = re.subn(r'\n$','',label)[0]
    label = re.subn('    ',r'........',label)[0] #identify indents
    label = re.subn('([^ ]+)('+16*r'\.'+')','\\1\n\\2',label)[0] #assume multiple indents means next line
    label = re.subn(r'_','\_',label)[0] #escape underscores
    fontsize = np.round(70./(len(label.split('\n'))-1))

    #Generate Plots
    fig = plt.gcf()
    fig.subplots_adjust(left=0.05,bottom=0.1,right=0.95,top=0.9,wspace=0.05,hspace=0.3)
    cavlabel = plt.subplot2grid((5,2),(0,0),title='\\Large{\\textbf{Group run: '+grprunname+'}}')
    rayplot = plt.subplot2grid((5,2),(1,0),rowspan=4,aspect='equal')
    pssplot = plt.subplot2grid((5,2),(0,1),rowspan=5,title='Poincare Surface of Section')
    fig.canvas.mpl_connect('pick_event',onpick)
    fig.canvas.mpl_connect('button_press_event',onpress)

    #Ray plot
    init = r'($r_0$, $\theta_0$, $\phi_0$) = (%0.4f,  %1.4f,  %2.4f)'%(r0val[0],theta0val[0],phi0val[0])\
           +'\n'+r'($\theta_0$, $\sin\chi_0$) = (%0.4f,  %1.4f)'%(pssvalval[0][0][0],pssvalval[0][0][1])
    rayplot.set(title='Ray Plot',xlabel=init,xticks=[],yticks=[])
    #Cavity background
    cavimg = plt.imread(os.path.join(grprunpath,'cavity.jpg'))
    rayplot.imshow(cavimg)
    #Rays foreground
    rayimg = plt.imread(os.path.join(grprunpath,str(imgid[0])+'.png'))
    rayplot.imshow(rayimg)
    

    #PSS plot and label
    plotpss(pssplot,pssvalval,bndfunc,nfunc)
    mpl.rc('text',usetex=True)
    cavlabel.text(0,1,label,ha='left',va='top',size=str(min(fontsize,12)))
    cavlabel.axis('off')
    
    return imgid, label0, pssvalval, r0val, theta0val, phi0val


#Analyse results ==============================================================
def colorfunc(clst):
    cnorm = mpl.colors.Normalize(-0.2,1.2)
    scalarmap = mpl.cm.ScalarMappable(norm=cnorm,cmap='spectral_r')
    return scalarmap.to_rgba(clst)[0:3]
    
def clusteredness(pssval,scaling=1e4):
    '''Heuristic function that characterizes the degree of clustering in a set of points,
    from 0 (localized on discrete points) to 1 (all over the place)'''
    pssval = np.array(pssval)
    bounces = len(pssval)
    clst = 2
    m = 2
    while m < 20:
        iclst = np.zeros(m);
        for i in range(m):
            pseudomean = np.angle(np.mean(np.exp(1j*pssval[i::m,0])))
            theta1 = (np.mod(pssval[i::m,0]-pseudomean+np.pi,2*np.pi)-np.pi)/np.pi
            x = 0.5*(np.var(theta1)+np.var(pssval[i::m,1]))*m/bounces
            iclst[i] = (1 - 1./(1+scaling*x))
            #iclst[i] = (1-np.exp(-10*x)) #Alternative function
        mclst = np.mean(iclst)
        if mclst < clst: clst = mclst
        m += 1
    return clst
    

#Execute ======================================================================
if __name__ == '__main__':
    #Display usage -----------------------------------------------------------
    if len(sys.argv) < 2 or sys.argv[1][0]=='h':
        print('A raytracing and plotting tool for studying the chaotic properties of light in a cavity.'\
              +'\n\nUsage:\ncavsim.py <options> [<args> ...]'\
              +'\n\nOptions:\n'+(6*'{}: {}\n\n')\
              .format('c [<theta0> <sinchi0> [<pars> ...]]','With no additional arguments, interactively configures simulation parameters. '\
                          +'Additional arguments will be taken as inputs for the initial condition <theta0> and <sinchi0>, '\
                          +'and parameters <pars> in the order: <bndname>, <bndpar>, <nname>, <npar>, <dt>, <tmax>.',\
                      'r [<grn> [<cfg>]]','Runs the simulation based on the configured parameters in the file <cfg>, '\
                          +'and records results with name <grn>.',
                      'd<nproc> [<grn> [<cfg>]]','Same as r, but distributes the run into <nproc> parallel processes (default 3).',\
                      'p [<grn>]','Presents an interactive plot of the results with name <grn>. With no argument, lists results available.',\
                      'h','Prints this page.',\
                      'v','Displays the version number.'))
    else:
        #Version -------------------------------------------------------------
        if sys.argv[1][0]=='v':
            print('0.0.2')


        #Configure parameters ------------------------------------------------
        runnow = 'n'
        if sys.argv[1][0]=='c':
            if len(sys.argv)<3:
                manualcfg()
            else:
                stdincfg()
            if (raw_input('\nRun now? ')+'n')[0].lower() == 'y':
                grprunname = raw_input('Enter a label for this cavity type: ')
                if grprunname == '': grprunname = None
                sys.argv[1:] = ['r',grprunname]

        #Config-based run ----------------------------------------------------
        if sys.argv[1][0]=='r':
            np.set_printoptions(threshold=np.NaN)
            grprunname = cfgrun(*sys.argv[2:])


        #Distributed run -----------------------------------------------------
        if sys.argv[1][0]=='d':
            if sys.argv[1][1:] == '':
                print('Defaulting to 3 processes...')
                parallel = 3
            else: parallel = eval(sys.argv[1][1:])

            #Set defaults
            if len(sys.argv)<3: sys.argv.append(str(int(time.time())))
            if len(sys.argv)<4: sys.argv.append('default')
            grprunpath = os.path.join('results','grprun_'+sys.argv[2])
            
            if type(parallel)==int:
                #Become Delegator
                #Start branches
                for proc in range(parallel):
                    os.system('start python cavsim.py d[{0},{1}] {2} {3}'\
                              .format(parallel,proc,sys.argv[2],sys.argv[3]))
                    print('Started branch {} of {}.'.format(proc+1,parallel))
                    time.sleep(5)

                #Check branches
                for proc in range(parallel):
                    while not os.path.exists(os.path.join(grprunpath,'result_'+sys.argv[2]+'_'+str(proc)+'_done.txt')):
                        time.sleep(5)
                    print('Branch {} of {} completed.'.format(proc+1,parallel))

                #Merge branches
                resultpath = os.path.join(grprunpath,'result_'+sys.argv[2])
                result = open(resultpath+'.txt','a')
                for proc in range(parallel):
                    respartpath = resultpath+'_'+str(proc)+'_done.txt'
                    respart = open(respartpath,'r')
                    result.write(''.join(list(respart)))
                    respart.close()
                    os.remove(respartpath)
                result.close()
                    
                
            elif type(parallel)==list and len(parallel)==2 and type(parallel[0])==type(parallel[1])==int:

                #Become Receiver
                np.set_printoptions(threshold=np.NaN)
                grprunname = cfgrun(*sys.argv[2:],parallel=parallel)
                respartpath = os.path.join(grprunpath,'result_'+sys.argv[2]+'_'+str(parallel[1]))
                os.rename(respartpath+'.txt',respartpath+'_done.txt')


        #Plot the results ----------------------------------------------------
        if sys.argv[1][0]=='p':
            if len(sys.argv)<3:
                #List available runs to plot
                print('Results available for plotting: \n')
                print(', '.join([grn[7:] for grn in os.listdir('results') if grn[:7]=='grprun_']))
                
            else:
                #Plot specified
                #Generate Plot
                grprunname = sys.argv[2]
                os.system('title {}'.format('Group run: '+grprunname))
                plt.figure('Group run - '+grprunname,figsize=(14,8))
                imgid, label, pssvalval, r0val, theta0val, phi0val = plotresults(grprunname)
                recomp = []
                plt.show()

                if len(recomp) != 0:
                    #Recompute new points
                    print('\nRecomputing new points...\n')
                    recomp = np.array(recomp)
                    cfg = open('config_default.txt','w')
                    cfg.write(label+'##Initial conditions'\
                              +'\ntheta0val = np.array('+str(recomp[:,0].tolist())+')'\
                              +'\nr0val,normangval = bndfunc(theta0val,shownormal=True)'
                              +'\nphi0val = normangval + np.array('+str(np.arcsin(recomp[:,1]).tolist())+')')
                    cfg.close()
                    os.system('python cavsim.py r '+grprunname+' default')
                    replotnow = raw_input('\nPlot results of '+grprunname+' again? ')+'n'
                    if replotnow[0] == 'y':
                        os.system('python cavsim.py p '+grprunname)


        #Debugging purposes --------------------------------------------------
        if sys.argv[1][0]=='t':
            addidx(sys.argv[2])
