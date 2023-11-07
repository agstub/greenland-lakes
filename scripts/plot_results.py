# This file contains the plot function for plotting the inversion results
# pass a list of time steps "timestep":
# if the length of timesteps is >1 then this will plot a png at each timestep (in movie directory)
# if the length of timesteps is == 1 then this will just plot a single png snapshot labeled "snap"
import sys,os
sys.path.insert(0, '../source')
from post_process import calc_dV_w,calc_dV_h
import numpy as np
import matplotlib.pyplot as plt
from params import data_dir,H,t0,x0,y0,Nt,Nx,Ny,x,y,results_dir,t,lake_name,u_d,v_d,t_r
from matplotlib.transforms import Bbox
from shapely.geometry import Point


def plot(t_ref,timesteps=range(Nt),h_lim=5,w_lim=5):
    if os.path.isdir(results_dir+'/movie')==False:
        os.mkdir(results_dir+'/movie')

    w_inv = np.load(results_dir+'/w_inv.npy')

    h_obs = np.load(data_dir+'/h_obs.npy')
    off_lake = np.load(data_dir+'/off_lake.npy')
    x_d = np.load(data_dir+'/x_d.npy')
    y_d = np.load(data_dir+'/y_d.npy')
    xp,yp = np.meshgrid(x_d,y_d)

    # define anomaly relative to reference time
    i0 = np.argmin(np.abs(t0-t_ref))
    h_ref = h_obs[i0,:,:] + 0*h_obs
    h_obs -= h_ref

    # define boundary for volume integration
    xc,yc = xp.mean(),yp.mean()
    alt_bdry = 0*xp+1
    rad = 3
    alt_bdry[np.sqrt((xp-xc)**2+(yp-yc)**2)>rad] = 0

    # calculate volume change time series:
    dV_inv = calc_dV_w(w_inv,alt_bdry)           # volume change from inversion
    dV_h = calc_dV_h(h_obs,alt_bdry)             # surface-based volume change
   
    # plot everything

    dV_max = np.max(np.array([np.max(dV_inv),np.max(dV_h)]))*(H**2)/1e9+np.sqrt(np.var(dV_inv))*(H**2)/1e9
    dV_min = np.min(np.array([np.min(dV_inv),np.min(dV_h)]))*(H**2)/1e9-np.sqrt(np.var(dV_inv))*(H**2)/1e9

    for i in timesteps:
        print('Saving image '+str(timesteps.index(i)+1)+' out of '+str(np.size(timesteps))+' \r',end='')

        fig = plt.figure(figsize=(12,12))
        plt.suptitle(r'data and inversion at $t=$'+format(t0[i],'.1f')+' yr',y=1.04,x=0.4,fontsize=26,bbox=dict(boxstyle='round', facecolor='w', alpha=1),clip_on=False)
        plt.subplot(212)
       
        if len(timesteps)>1:
            j = i
        else:
            j = -1    

        dV_invs = dV_inv[0:j]*(H**2)/1e9
        dV_alt = dV_h[0:j]*(H**2)/1e9
        label = r'ICESat-2 ($\Delta V_\mathrm{alt}$)'
        plt.plot(t0[0:j],dV_alt,color='indianred',linewidth=5,label=label)
        plt.plot(t0[0:j],dV_invs,color='royalblue',linewidth=4,label=r'inversion ($\Delta V_\mathrm{inv}$)')
        plt.ylabel(r'$\Delta V$ (km$^3$)',fontsize=20)
        plt.xlim(0,t0.max())
        plt.ylim(dV_min-0.1,0.9*dV_max+0.1)
        plt.xlabel(r'$t$ (yr)',fontsize=20)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=20,ncol=2,bbox_to_anchor=(0.875,-0.15))
        plt.title(r'volume change time series', fontsize=26,bbox=dict(boxstyle='round', facecolor='w', alpha=1),y=0.91,x=0.3075)


        plt.subplot(221)
        p=plt.contourf(xp,yp,h_obs[i,:,:],cmap='PuOr_r',extend='both',levels=h_lim*np.linspace(-1,1,11))
        circle = plt.Circle((xc, yc), rad, color='k',fill=False,linewidth=4)
        plt.gca().add_patch(circle)
        plt.ylabel(r'$y$ (km)',fontsize=20)
        plt.xlabel(r'$x$ (km)',fontsize=20)
        plt.xlim(x_d.min(),x_d.max())
        plt.ylim(y_d.min(),y_d.max())
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.gca().set_aspect('equal', 'box')
        cbar_ax = fig.add_axes([0.13, 0.9, 0.35, 0.02])
        cbar = fig.colorbar(p,orientation='horizontal',cax=cbar_ax)
        cbar.set_label(r'$\Delta h_\mathrm{a}$ (m)',fontsize=24,labelpad=15)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')

        plt.subplot(222)
        p=plt.contourf(xp,yp,w_inv[i,:,:],cmap='PuOr_r',levels=w_lim*np.linspace(-1,1,11),extend='both')
        circle = plt.Circle((xc, yc), rad, color='k',fill=False,linewidth=4)
        plt.gca().add_patch(circle)

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
        if len(timesteps)>1:
            plt.savefig(results_dir+'/movie/'+str(i),bbox_inches=Bbox([[0.15,-0.25],[12.5,14]]))
        elif len(timesteps)==1:    
            plt.savefig(results_dir+'/'+lake_name+'_snap',bbox_inches=Bbox([[0.15,-0.25],[12.5,14]]))
            # plt.tight_layout()
            plt.show()
        plt.close()