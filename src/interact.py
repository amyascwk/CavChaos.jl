import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys, re


#==============================================================================
#Matplotlib connect functions


#PSS plot picker (focus results from selected run)
def onpick(ev):
    axs = plt.figure("PSS explorer").get_axes()
    mev = ev.mouseevent
    art = ev.artist
    cavdir = "../results/"+re.subn(r'\\',r'',axs[0].get_title()[23:-2])[0]
    
    #Read new ray data
    emph = art.get_label()[5:]
    f = open(cavdir+"/"+emph.zfill(4)+".txt","r")
    lines = f.readlines()
    f.close()
    
    if mev.button == 1: #left click
        #Find old artist
        for a in axs[2].findobj(mpl.lines.Line2D):
            if a.get_ms() == 8: #emphasized previously
                unemph = a.get_label()[5:]
                #Read old ray data
                f0 = open(cavdir+"/"+unemph.zfill(4)+".txt","r") 
                lines0 = f0.readlines()
                f0.close()
                
                #Get color for old artist
                bounces = readarray2d(lines0[7][7:-1])
                uclst,order = unclusteredness(bounces)
                color = colorfunc(uclst)
                print('{!s}: uclst is {:.8f} for order {:d}\n'.format(unemph.zfill(4),uclst,order))
                #De-emphasize old artist
                a.set(mec=color,mfc=color,mew=None,ms=2,zorder=0)
        
        #Emphasize new artist
        art.set(mec='k',mfc='r',mew=1.5,ms=8,zorder=1000)
        
        #Display rayplot corresponding to current artist
        axs[1].get_children()[3].remove() #Remove old ray image
        init = readarray1d(lines[1][7:-1])
        bounces = readarray2d(lines[7][7:-1])
        initlabel = "Run {!s}: $(r_0,\\theta_0,\\phi_0)$ = ({:.8f}, {:.8f}, {:.8f})\n$(\\theta_0,\\sin\\chi_0)$ = ({:.8f}, {:.8f})".format(emph.zfill(4),init[0],init[1],init[2],bounces[0,0], np.sin(bounces[0,1]))
        axs[1].set_xlabel(initlabel) #Set new ray plot label
        rayimg = plt.imread(cavdir+"/"+emph.zfill(4)+".png") #Read new ray image
        axs[1].imshow(rayimg) #Show new ray image
        
        #Edit mode label
        modelabel = 'Selected path does not appear to\nbe near to a periodic mode.'
        if len(lines) >= 12:
            uclst, order = eval(lines[9][:-1]);
            dwdm, omega0_s, omega0_p = eval(re.subn(r'im',r'j',lines[11][:-1])[0])
            modelabel = 'Periodic mode of order {:d} ({:.4f}\%)'.format(order,50+500*(0.1-uclst))
            modelabel = modelabel+'\n$\\omega_{{m,s}}$ = {:.8f}*m + ({:.8f})'.format(dwdm,omega0_s)
            modelabel = modelabel+'\n$\\omega_{{m,p}}$ = {:.8f}*m + ({:.8f})'.format(dwdm,omega0_p)
            modelabel = modelabel+'\nFarfield directionality: Not computed yet'
        
        axs[0].get_children()[3].set_text(modelabel)
        
        plt.draw()
        
        #Farfield emission plot
        #Warning: Still testing, not properly normalized, not even between both     
        #directions
        if len(lines) >= 18:
            phi_t = readarray1d(lines[13][7:-1])
            output_s = readarray1d(lines[15][7:-1])
            output_p = readarray1d(lines[17][7:-1])
            
            fig2 = plt.figure('farfield')
            ax2 = plt.gca()
            ax2.cla()
            for i in range(len(phi_t)):
                ax2.plot([0,output_s[i]*np.cos(phi_t[i])],[0,output_s[i]*np.sin(phi_t[i])],'r')
                ax2.plot([0,output_p[i]*np.cos(phi_t[i])],[0,output_p[i]*np.sin(phi_t[i])],'b')
            ax2.axis('equal')
            
            plt.draw()
        
    elif mev.button == 3: #right click
        #Open to other features
        sys.stdout.write("\n")


#==============================================================================
#Miscellaneous tools


#Conversion of results data in julia syntax to ndarray
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


#Unclusteredness heuristic measurement of degree of unclusteredness
def unclusteredness(bounces,scaling=1e4):
    #Heuristic function that characterizes the degree of unclustering in a set of points,
    #from 0 (localized on discrete points) to 1 (all over the place)
    
    if len(bounces) < 40:
        return 1.0, -1
    
    uclst = 2.0 #initiate unclusteredness score
    order = 20
    for m in range (2,21):
        iuclst = np.zeros(m);
        for i in range(m):
            #Try to invent a notion of distance in modular arithmetic
            pseudomean = np.angle(np.mean(np.exp(1j*bounces[i::m,0])));
            u = (np.mod(bounces[i::m,0]-pseudomean+np.pi,2*np.pi)-np.pi)/np.pi;
            x = 0.5*(np.var(u)+np.var(np.sin(bounces[i::m,1])))*m/len(bounces);
            iuclst[i] = (1.0 - 1.0/(1+scaling*x));
            #iuclst[i] = (1-exp(-10*x)) #Alternative function
        muclst = np.mean(iuclst);
        if muclst < uclst:
            #Record smallest
            uclst = muclst;
            order = m;
    return uclst, order
