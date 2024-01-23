# This file contains the plot function for plotting the inversion results
# pass a list of time steps "timestep":
# if the length of timesteps is >1 then this will plot a png at each timestep (in movie directory)
# if the length of timesteps is == 1 then this will just plot a single png snapshot labeled "snap"
import sys,os
sys.path.insert(0, '../source')
from post_process import calc_dV, calc_dQ
import numpy as np
import matplotlib.pyplot as plt
from params import data_dir,H,t0,Nt,results_dir,lake_name
import matplotlib as mpl
from matplotlib.transforms import Bbox
from shapely.geometry import Point
import geopandas as gpd

def plot(outlines,t_ref,timesteps):
    if os.path.isdir(results_dir+'/movie')==False:
        os.mkdir(results_dir+'/movie') 

    w_inv = np.load(results_dir+'/w_inv.npy')
    h_obs = np.load(data_dir+'/h_obs.npy')
    x_d = np.load(data_dir+'/x_d.npy')
    y_d = np.load(data_dir+'/y_d.npy')
    xp,yp = np.meshgrid(x_d,y_d)


    h_lim = np.floor(np.max(np.abs(h_obs)))
    w_lim = np.floor(np.max(np.abs(w_inv)))   

    # define anomaly relative to reference time
    i0 = np.argmin(np.abs(t0-t_ref))
    h_ref = h_obs[i0,:,:] + 0*h_obs
    h_obs -= h_ref

    num_lakes = outlines.size

    dV = np.zeros((num_lakes,Nt))
    dQ = np.zeros((num_lakes,Nt))

    # define boundary for volume integration
    for l in  range(num_lakes):
        outline = outlines[l]
        bdry = 0*xp
        for i in range(x_d.size):
            for j in range(y_d.size):
                point = Point(x_d[i],y_d[j])
                bdry[i,j] = outline.contains(point)

        # calculate volume change time series:
        dV[l,:] = calc_dV(w_inv,bdry)           # volume change from inversion
    
        dQ[l,:] = calc_dQ(w_inv,bdry)


    colors =  mpl.cm.Dark2(np.linspace(0, 1, num_lakes))
    # plot everything
    for i in timesteps:

        fig = plt.figure(figsize=(8,14))
        plt.suptitle(r'$t=$'+format(t0[i],'.1f')+' yr',fontsize=26,y=1.02,bbox=dict(boxstyle='round', facecolor='w', alpha=1))
       
        plt.subplot(312)

        for l in range(num_lakes):
            dV_l = dV[l,:]*(H**2)/1e9
            plt.plot(t0[0:j],dV_l,color=colors[l],linewidth=3,label=r'$\Delta V_\mathrm{inv}$')
        plt.axhline(0,color='k',linestyle='--')
        plt.ylabel(r'$\Delta V$ (km$^3$)',fontsize=20)
        plt.xlim(0,t0.max())
        plt.gca().xaxis.set_ticklabels([])
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        plt.subplot(313)
  
        for l in range(num_lakes):
            dQ_l = dQ[l,:]*(H**2)/3.154e7
            plt.plot(t0,dQ_l,color=colors[l],linewidth=3)
        plt.axhline(0,color='k',linestyle='--')
        plt.ylabel(r'$\Delta Q$ (m$^3$ / s)',fontsize=20)
        plt.xlim(0,t0.max())
        plt.xlabel(r'$t$ (yr)',fontsize=20)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)

        plt.subplot(321)
        p=plt.contourf(xp,yp,h_obs[i,:,:],cmap='bwr_r',extend='both',levels=h_lim*np.linspace(-1,1,11))
        for l in range(num_lakes):
            outline = gpd.GeoSeries(outlines[l])
            outline.plot(edgecolor=colors[l],facecolor='none',ax=plt.gca(),linewidth=4)
        plt.ylabel(r'$y$ (km)',fontsize=20)
        plt.xlabel(r'$x$ (km)',fontsize=20)
        plt.xlim(x_d.min(),x_d.max())
        plt.ylim(y_d.min(),y_d.max())
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.gca().set_aspect('equal', 'box')
        cbar_ax = fig.add_axes([0.13, 0.9, 0.35, 0.02])
        cbar = fig.colorbar(p,orientation='horizontal',cax=cbar_ax)
        cbar.set_label(r'$\Delta h_\mathrm{obs}$ (m)',fontsize=24,labelpad=15)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')

        plt.subplot(322)
        p=plt.contourf(xp,yp,w_inv[i,:,:],cmap='bwr_r',levels=w_lim*np.linspace(-1,1,11),extend='both')
        for l in range(num_lakes):
            outline = gpd.GeoSeries(outlines[l])
            outline.plot(edgecolor=colors[l],facecolor='none',ax=plt.gca(),linewidth=4)
        plt.xlabel(r'$x$ (km)',fontsize=20)
        plt.ylabel(r'$y$ (km)',fontsize=20)
        plt.xlim(x_d.min(),x_d.max())
        plt.ylim(y_d.min(),y_d.max())
        plt.gca().yaxis.tick_right()
        plt.gca().yaxis.set_label_position("right")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.gca().set_aspect('equal', 'box')
        cbar_ax = fig.add_axes([0.55, 0.9, 0.35, 0.02])
        cbar = fig.colorbar(p,orientation='horizontal',cax=cbar_ax)
        cbar.set_label(r'$w_\mathrm{b}$ (m/yr)',fontsize=24,labelpad=15)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')
        plt.show()  
        plt.close()