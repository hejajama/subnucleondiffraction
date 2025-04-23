import math
import numpy as np
import scipy.ndimage
import scipy.interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

zmax = 8.5
zmin = -zmax
rmax = 8.5
rmin = 0.0
scutoff = 0.0009
tcutoff = 0.0009
zpoints = 400
rpoints = 400
spinlw = 2.0

cmap='YlGn'

vmax = 0.1

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 18,
        }

labelsize = 15
boxprop = dict(boxstyle='round',facecolor='white')



def plotden(fname,txt):
    rc1,zc1,rhn1  = np.loadtxt('den_neutron_208Pb_'+fname+'.dat',dtype='float',usecols=[0,1,2], comments='#').T
    rc1,zc1,rhp1  = np.loadtxt('den_proton_208Pb_'+fname+'.dat',dtype='float',usecols=[0,1,2], comments='#').T

    # dublicate to -z values.
    rcb1 = np.append(rc1,rc1) ;    zcb1 = np.append(zc1,-zc1) ;
    rhn1 = np.append(rhn1,rhn1) ;  rhp1 = np.append(rhp1,rhp1);

    # dublicate to -r_rho values.
    rcb1 = np.append(rcb1,-rcb1) ; zcb1 = np.append(zcb1,zcb1) ;
    rhn1 = np.append(rhn1,rhn1)  ; rhp1 = np.append(rhp1,rhp1);

    # Set up a regular grid of interpolation points
    ri, zi = np.linspace(rmin, rmax, rpoints), np.linspace(zmin, zmax, zpoints)
    extent=(rmin, rmax, zmin, zmax)
    ri, zi = np.meshgrid(ri, zi)

    # use griddata
    rhon = scipy.interpolate.griddata((rcb1,zcb1),rhn1,(ri,zi),method='cubic')
    rhop = scipy.interpolate.griddata((rcb1,zcb1),rhp1,(ri,zi),method='cubic')

    # interpolation along z=0 line
    rhonln = scipy.interpolate.Rbf(rcb1,zcb1,rhn1)
    rhopln = scipy.interpolate.Rbf(rcb1,zcb1,rhp1)

    rvals = np.linspace(0.0,rmax,num=rpoints)
    rhonz0 = np.linspace(0.0,rmax,num=rpoints)
    rhopz0 = np.linspace(0.0,rmax,num=rpoints)
    for i in range(0,rpoints):
        rhonz0[i] = rhonln(rvals[i],0.0)
        rhopz0[i] = rhopln(rvals[i],0.0)


    # with colorbar
    fig = plt.figure(figsize=(9.2,8.5))
    gs = fig.add_gridspec(nrows=4,ncols=2,left=0.11,right=0.8,hspace=0.54,wspace=0.25,bottom=0.08,top=0.982)
    ax00 = fig.add_subplot(gs[:-1,:-1])
    ax01 = fig.add_subplot(gs[:-1,-1])
    ax10 = fig.add_subplot(gs[-1,:-1])
    ax11 = fig.add_subplot(gs[-1,-1])

    ax00.set_aspect('equal')
    ax00.set_xlim([rmin,rmax]) ; ax00.set_ylim([zmin,zmax]) ;
    ax00.pcolormesh(ri,zi,rhon,cmap=cmap,vmin=0.0,vmax=vmax,shading='auto')
    ax00.contour(ri, zi, rhon,levels=[0.02,0.04,0.06,0.08],colors='white',linewidths=2.5,linestyles='dashed')
    ax00.set_ylabel("z (fm)",fontdict=font)
    ax00.set_xlabel("y (fm)",fontdict=font)

    ax00.text(8.0,7.3,txt+r', ${\rho_{\rm n}}$',fontdict=font,bbox=boxprop,ha='right')

    ax01.set_aspect('equal')
    ax01.set_xlim([rmin,rmax]) ; ax01.set_ylim([zmin,zmax]) ;
    ax3 = ax01.pcolormesh(ri,zi,rhon,cmap=cmap,vmin=0.0,vmax=vmax,shading='auto')
    ax01.contour(ri, zi, rhop,levels=[0.02,0.04,0.06,0.08],colors='white',linewidths=2.5,linestyles='dashed')
    #ax01.set_ylabel("z (fm)",fontdict=font)
    ax01.set_xlabel("y (fm)",fontdict=font)
    ax01.text(8.0,7.3,txt+r', ${\rho_{\rm p}}$',fontdict=font,bbox=boxprop,ha='right')


    ax10.set_ylim([0,0.1]) ;
    ax10.set_xlim([0.0,rmax]) ;
    ax10.set_ylabel(r'$\rho(r_\rho,z=0)$',fontdict=font)
    ax10.set_xlabel(r'$r_\rho$ (fm)',fontdict=font)
    ax10.plot(rvals,rhonz0,color='tab:green',linewidth=2.0)
    ax10.text(8.1,0.08,r'${\rho_{\rm n}}$',fontdict=font,ha='right')


    ax11.set_ylim([0,0.1]) ;
    ax11.set_xlim([0.0,rmax]) ;
    ax11.set_xlabel(r'$r_\rho$ (fm)',fontdict=font)
    ax11.plot(rvals,rhopz0,color='tab:orange',linewidth=2.0)
    ax11.text(8.1,0.08,r'${\rho_{\rm p}}$',fontdict=font,ha='right')

    cbar_ax = fig.add_axes([0.82,0.33,0.025,0.649])
    cbar = fig.colorbar(ax3,cax=cbar_ax,shrink=1.0,ticks=[0,0.05,0.1])
    cbar.ax.tick_params(labelsize=labelsize)

    plt.savefig('densities_208Pb_'+fname+'.png',format='png',dpi=150)
    plt.show()

plotden('UNEDF0',r'$^{208}$Pb, UNEDF0')

