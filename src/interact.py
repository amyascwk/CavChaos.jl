import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys, re


#==============================================================================
#Matplotlib connect functions
#Size code: 2 is normal, 8 is in focus


#PSS plot picker (focus results from selected run)
def onpick(ev):
    
    axs = plt.figure("PSS Explorer").get_axes()
    mev = ev.mouseevent
    art = ev.artist
    cavdir = "../results/"+re.subn(r'\\',r'',axs[0].get_title()[23:-2])[0]
    
    #Read new ray data
    emph = eval(art.get_label())[0]
    f = open(cavdir+"/"+emph+".txt","r")
    lines = f.readlines()
    f.close()
    
    if mev.button == 1: #left click
        #Find old artist
        for a in axs[2].findobj(mpl.lines.Line2D):
            if a.get_ms() == 8: #emphasized previously
                #De-emphasize old artist
                unemph,uclst = eval(a.get_label())
                color = colorfunc(uclst)
                a.set(mec=color,mfc=color,mew=None,ms=2,zorder=0)
                break
        
        #Emphasize new artist
        art.set(mec='k',mfc='r',mew=1.5,ms=8,zorder=1000)
        
        #Display rayplot corresponding to current artist
        axs[1].get_children()[3].remove() #Remove old ray image
        init = readarray1d(lines[1][7:-1]) #(removed starting 'Float64', trailing '\n')
        bounces = readarray2d(lines[7][7:-1])
        initlabel = "Run {!s}: ($r_0$, $\\theta_0$, $\\phi_0$) = ({:.8f}, {:.8f}, {:.8f})\n($\\theta_0$, $\\sin\\chi_0$) = ({:.8f}, {:.8f})".format(emph,init[0],init[1],init[2],bounces[0,0], np.sin(bounces[0,1]))
        axs[1].set_xlabel(initlabel) #Set new ray plot label
        rayimg = plt.imread(cavdir+"/"+emph+".png") #Read new ray image
        axs[1].imshow(rayimg) #Show new ray image
        
        #Edit mode label
        modelabel = 'Selected path does not appear to\nbe near to a periodic mode.'
        if len(lines) >= 12:
            uclst, order = eval(lines[9][:-1]);
            dwdm, omega0_s, omega0_p = eval(re.subn(r'im',r'j',lines[11][:-1])[0])
            modelabel = 'Periodic mode of order {:d} ({:.4f}\%)'.format(order,50+500*(0.1-uclst))
            modelabel = modelabel+'\n$\\omega_{{m,s}}$ = {:.8f}*m + ({:.8f})'.format(dwdm,omega0_s)
            modelabel = modelabel+'\n$\\omega_{{m,p}}$ = {:.8f}*m + ({:.8f})'.format(dwdm,omega0_p)
            #modelabel = modelabel+'\nFarfield directionality: '
        
        axs[0].get_children()[3].set_text(modelabel)
        
        plt.draw()
        
        #Announce on stdout
        uclst, order = eval(lines[9][:-1]);
        print('\n{!s}: uclst is {:.8f} for order {:d}'.format(emph,uclst,order))
        
        #Farfield emission plot
        #Warning: Still testing
        if len(lines) >= 18:
            phi_t = readarray1d(lines[13][7:-1])
            output_s = readarray1d(lines[15][7:-1])
            output_p = readarray1d(lines[17][7:-1])
            
            fig2 = plt.figure('Far Field Emissions (still testing)')
            ax2 = plt.gca()
            ax2.cla()
            for i in range(len(phi_t)):
                ax2.plot([0,output_s[i]*np.cos(phi_t[i])],[0,output_s[i]*np.sin(phi_t[i])],'r')
                ax2.plot([0,output_p[i]*np.cos(phi_t[i])],[0,output_p[i]*np.sin(phi_t[i])],'b')
            ax2.axis('equal')
            
            plt.draw()
        
    elif mev.button == 3: #right click
        #Open for other features (more plots, etc)
        sys.stdout.write("\n")


#PSS plot clicker (get coordinates and tag modes)
def onpress(mev):
    axs = plt.figure("PSS Explorer").get_axes()
    cavdir = "../results/"+re.subn(r'\\',r'',axs[0].get_title()[23:-2])[0]
    
    if mev.button == 1:
        print('\nMouseclick at (x,y) = ('+str(mev.xdata)+', '+str(mev.ydata)+')')
    elif mev.button == 2:
        for a in axs[2].findobj(mpl.lines.Line2D):
            if a.get_ms() == 2 or a.get_ms() == 7:
                runnum,uclst = eval(a.get_label())
                if a.get_ms() == 7:
                    color = colorfunc(uclst)
                    a.set(mec=color,mfc=color,mew=None,ms=2,zorder=0)
                elif ((a.get_ms() == 2) and (uclst < 1e-4)):
                    a.set(mec='k',mfc='m',mew=1,ms=7,zorder=1000)
        plt.draw()


#PSS plot scroller (Chaotic sea visibility control)
def onscroll(mev):
    axs = plt.figure("PSS Explorer").get_axes()
    cavdir = "../results/"+re.subn(r'\\',r'',axs[0].get_title()[23:-2])[0]
    patch = axs[2].findobj(mpl.patches.Rectangle)[0]
    x = int(patch.get_label())
    newx = min(max(x + mev.step,0),10)
    
    #Update
    maxuclst = newx*newx*newx*0.001 + 1e-4
    for a in axs[2].findobj(mpl.lines.Line2D):
        if a.get_ms() == 2 or a.get_ms() == 8:
            runnum,uclst = eval(a.get_label())
            if uclst <= maxuclst:
                a.set_visible(True)
            else:
                a.set_visible(False)
    plt.draw()
    patch.set_label(str(newx))
    

#==============================================================================
#Miscellaneous tools


#Conversion of results data in julia syntax to numpy ndarray
def readarray1d(s0):
    s = s0
    s = re.subn(r'NaN',r'np.nan',s)[0]
    return np.array(eval(s))

def readarray2d(s0):
    s = s0
    s = re.subn(r'\[',r'[[',s)[0]
    s = re.subn(r'\]',r']]',s)[0]
    s = re.subn(r';',r'],[',s)[0]
    s = re.subn(r'(\d)\s+([\.\d-])',r'\1,\2',s)[0]
    s = re.subn(r'NaN',r'np.nan',s)[0]
    return np.array(eval(s))


#PSS plot colormap function
def colorfunc(uclst):
    cnorm = mpl.colors.Normalize(-0.2,1.2)
    scalarmap = mpl.cm.ScalarMappable(norm=cnorm,cmap='spectral_r')
    return scalarmap.to_rgba(uclst)[0:3]
