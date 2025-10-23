import matplotlib.animation as animation
from matplotlib.patches import ConnectionPatch
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import warnings
import ffmpeg
import glob
import re
warnings.filterwarnings('ignore')

def load_file(filename, Nx, Nz):
    file = h5py.File(filename, "r")
    utot = file['timeseries']['u']
    keys = utot.keys()
    key_list = [key for key in keys]
    key_list = key_list[:-1]
    def keyfunc(s):
        return int(s)
    sorted_keys = sorted(key_list, key=keyfunc)
    x = file['grid']['underlying_grid']['xᶜᵃᵃ']
    z = file['grid']['underlying_grid']['z']['cᵃᵃᶜ']
    x_array = np.array(x)
    z_array = np.array(z)
    bh = file['grid']['immersed_boundary']['bottom_height']
    bh_array = np.array(bh).reshape(Nx+8,)
    X, Z = np.meshgrid(x_array, z_array)
    Z_masked = Z
    for indx in range(len(x_array)):
        zz = Z[:, indx]
        zz[zz<bh_array[indx]] = np.nan
        Z_masked[:, indx] = zz


    time_array = np.zeros((len(sorted_keys),))
    for ind in range(len(sorted_keys)):
        time = file['timeseries']['t'][sorted_keys[ind]]
        time_array[ind] = np.array(time)/60
    
    return file, sorted_keys, x_array, z_array, bh_array, Z_masked, time_array

def extract_variables(file, sorted_keys, x_array, z_array, Z_masked, Nx, Nz, U0, key_ind=-1):
    #u = file['timeseries']['u′'][sorted_keys[key_ind]]
    utot = file['timeseries']['u'][sorted_keys[key_ind]]
    #u_array = np.array(u[:]).reshape(Nz+8,Nx+8)
    utot_array = np.array(utot[:]).reshape(Nz+8,Nx+8)

    w = file['timeseries']['w'][sorted_keys[key_ind]]
    w_array = np.array(w[:]).reshape(Nz+9,Nx+8)[:-1, :]

    N2 = file['timeseries']['N2'][sorted_keys[key_ind]]
    N2_array = np.array(N2[:]).reshape(Nz+9,Nx+8)[:-1, :]

    #u_array[np.isnan(Z_masked)] = np.nan
    utot_array[np.isnan(Z_masked)] = np.nan
    w_array[np.isnan(Z_masked)] = np.nan
    N2_array[np.isnan(Z_masked)] = np.nan


    return utot_array, w_array, N2_array

def ke_omega(file, sorted_keys, x_array, z_array, Z_masked, Nx, Nz, U0, key_ind=-1):
    #u = file['timeseries']['u′'][sorted_keys[key_ind]]
    utot = file['timeseries']['u'][sorted_keys[key_ind]]
    #u_array = np.array(u[:]).reshape(Nz+8,Nx+8)
    utot_array = np.array(utot[:]).reshape(Nz+8,Nx+8)

    w = file['timeseries']['w'][sorted_keys[key_ind]]
    w_array = np.array(w[:]).reshape(Nz+9,Nx+8)[:-1, :]

    #u_array[np.isnan(Z_masked)] = np.nan
    utot_array[np.isnan(Z_masked)] = np.nan
    w_array[np.isnan(Z_masked)] = np.nan

    ke = 0.5*((utot_array-U0)**2 + w_array**2)
    ke[np.isnan(Z_masked)] = np.nan

    wz, wx = np.gradient(w_array, z_array, x_array)
    uz, ux = np.gradient(utot_array, z_array, x_array)
    omega = (uz-wx)**2
    omega[np.isnan(Z_masked)] = np.nan


    return ke, omega

def save_fig_uwN(file, sorted_keys, x_array, z_array, bh_array, Z_masked, z0, Nx, Nz, time_array, U0, i, foldername):
    fig, (ax0, ax1, ax2) = plt.subplots(
    ncols=3,
    sharey=True,
    figsize=(18, 4)
    )
    utot_array, w_array, N2_array = extract_variables(file, sorted_keys, x_array, z_array,  Z_masked, Nx, Nz, U0,  key_ind=i)
    H = z_array.min()
    
    im = ax0.pcolor(x_array, H-z_array, utot_array-U0, vmin=-0.2, vmax=0.2, cmap='bwr')
    ax0.plot(x_array, H-bh_array, '-k', lw=2.0)
    ax0.plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax0)
    ax0.set(xlabel='$x$ (m)', ylabel='$z$ (m)', title='$u-U$ (m/s): '+ str(int(time_array[i]))+' mins')
    #ax0.set_ylim([-200, -150])
    
    
    im = ax1.pcolor(x_array, H-z_array, w_array, vmin=-0.1, vmax=0.1, cmap='bwr')
    ax1.plot(x_array, H-bh_array, '-k', lw=2.0)
    ax1.plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax1)
    ax1.set(xlabel='$x$ (m)', title='$w$ (m/s): '+ str(int(time_array[i]))+' mins')
    #ax1.set_ylim([-200, -150])
    
    
    im = ax2.pcolor(x_array, H-z_array, N2_array, vmin=-0.003, vmax=0.003, cmap='bwr')
    ax2.plot(x_array, H-bh_array, '-k', lw=2.0)
    ax2.plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax2)
    ax2.set(xlabel='$x$ (m)', title='$N^2\, (1/\mathrm{s}^2)$: '+ str(int(time_array[i]))+' mins')


    filename = foldername+'uwN_'+str(int(time_array[i]))+'.png'
    fig.savefig(filename)
    plt.close('all')

def save_fig_uwN_zoomin(file, sorted_keys, x_array, z_array, bh_array, Z_masked, z0, h0, Nx, Nz, time_array, U0, i, foldername):
    zlim = max(-60, 10*h0) 
    
    fig, ax = plt.subplots(2, 3, sharex=True, figsize=(18, 8))
    utot_array, w_array, N2_array = extract_variables(file, sorted_keys, x_array, z_array,  Z_masked, Nx, Nz, U0, key_ind=i)
    H = z_array.min()
    
    im = ax[0,0].pcolor(x_array, H-z_array, utot_array-U0, vmin=-0.05, vmax=0.05, cmap='bwr')
    ax[0,0].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[0,0].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[0,0])
    ax[0,0].set( ylabel='$z$ (m)', title='$u-U$ (m/s): '+ str(int(time_array[i]))+' mins')
    #ax[0,0].set_aspect('equal', 'box')
    #ax0.set_ylim([-200, -150])
    
    im = ax[0,1].pcolor(x_array, H-z_array, w_array, vmin=-0.01, vmax=0.01, cmap='bwr')
    ax[0,1].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[0,1].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[0,1])
    ax[0,1].yaxis.set_ticklabels([])
    ax[0,1].set(title='$w$ (m/s): '+ str(int(time_array[i]))+' mins')
    #ax[0,1].set_aspect('equal', 'box')
    #ax1.set_ylim([-200, -150])
    
    im = ax[0,2].pcolor(x_array, H-z_array, N2_array, vmin=-0.001, vmax=0.001, cmap='bwr')
    ax[0,2].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[0,2].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[0,2])
    ax[0,2].yaxis.set_ticklabels([])
    ax[0,2].set(title='$N^2\, (1/\mathrm{s}^2)$: '+ str(int(time_array[i]))+' mins')
    #ax[0,2].set_aspect('equal', 'box')

    im = ax[1,0].pcolor(x_array, H-z_array, utot_array-U0, vmin=-0.05, vmax=0.05, cmap='bwr')
    ax[1,0].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[1,0].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[1,0])
    ax[1,0].set(xlabel='$x$ (m)', ylabel='$z$ (m)')
    ax[1,0].set_ylim([zlim, 0])
    #ax[1,0].set_aspect('equal', 'box')
    
    im = ax[1,1].pcolor(x_array, H-z_array, w_array, vmin=-0.01, vmax=0.01, cmap='bwr')
    ax[1,1].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[1,1].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[1,1])
    ax[1,1].yaxis.set_ticklabels([])
    ax[1,1].set(xlabel='$x$ (m)')
    ax[1,1].set_ylim([zlim, 0])
    #ax[1,1].set_aspect('equal', 'box')
    
    im = ax[1,2].pcolor(x_array, H-z_array, N2_array, vmin=-0.001, vmax=0.001, cmap='bwr')
    ax[1,2].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[1,2].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[1,2])
    ax[1,2].yaxis.set_ticklabels([])
    ax[1,2].set(xlabel='$x$ (m)')
    ax[1,2].set_ylim([zlim, 0])
    #ax[1,2].set_aspect('equal', 'box')

    filename = foldername+'uwN_zoomin_'+str(int(time_array[i]))+'.png'
    fig.savefig(filename)
    plt.close('all')

def calc_cell_vol(x_array, z_array, Nx, Nz, Z_masked):
    dx = ((x_array.max()-x_array.min())/len(x_array))*np.ones((Nz+8, Nx+8))
    dz = ((z_array.max()-z_array.min())/len(z_array))*np.ones((Nz+8, Nx+8))
    dz[np.isnan(Z_masked)] = np.nan
    cell_vol = dx*dz
    tot_vol = np.nansum(cell_vol)
    return cell_vol, tot_vol

def calc_histogram(x_array, z_array, Nx, Nz, Z_masked, var1, var2, Nbins1, Nbins2, cell_vol, tot_vol, var1_lim1=-1, var1_lim2=-1, var2_lim1=-1, var2_lim2=-1):
    #cell_vol, tot_vol = calc_cell_vol(x_array, z_array, Nx, Nz, Z_masked)
    if var1_lim1 == -1:
        var1_lim1 = np.nanmax(var1)
    if var1_lim2 == -1:
        var1_lim2 = np.nanmin(var1)
    if var2_lim1 == -1:
        var2_lim1 = np.nanmax(var2)
    if var2_lim2 == -1:
        var2_lim2 = np.nanmin(var2)
    
    dv1 = (var1_lim1-var1_lim2)/Nbins1
    dv2 = (var2_lim1-var2_lim2)/Nbins2
    
    v1_i = var1_lim2 + 0.5*(2*np.arange(1, Nbins1+1, 1)-1)*dv1
    v2_j = var2_lim2 + 0.5*(2*np.arange(1, Nbins2+1, 1)-1)*dv2
    I_ij = np.zeros((Nbins1, Nbins2))

    for ind1 in range(Nbins1):
        for ind2 in range(Nbins2):
            valid_ind = np.where((var1>=v1_i[ind1]-0.5*dv1) & (var1<v1_i[ind1]+0.5*dv1) \
                                 & (var2>=v2_j[ind2]-0.5*dv2) & (var2<v2_j[ind2]+0.5*dv2)) 
            I_ij[ind1, ind2] = np.nansum(cell_vol[valid_ind])

    return v1_i, v2_j, I_ij

def save_fig_USP(file, sorted_keys, x_array, z_array, bh_array, Z_masked, z0, Nx, Nz, cell_vol, tot_vol, U0, time_array, i, foldername):
    fig, ax = plt.subplots(2, 2, figsize=(15, 10))
    _, _, N2_array, ke, omega = extract_variables(file, sorted_keys, x_array, z_array,  Z_masked, Nx, Nz, U0,  key_ind=i)
    N2_i, KE_j, I_N_KE = calc_histogram(x_array, z_array, Nx, Nz, Z_masked, N2_array, ke, 100, 100, cell_vol, tot_vol, var1_lim1=0.004, var1_lim2=-0.002, var2_lim1=0.02,  var2_lim2=0)
    N2_i, omega_j, I_N_omega = calc_histogram(x_array, z_array, Nx, Nz, Z_masked, N2_array, omega, 100, 150, cell_vol, tot_vol, var1_lim1=0.004, var1_lim2=-0.002, var2_lim1=0.02, var2_lim2=0)
    H = z_array.min()

    im = ax[0,0].pcolor(x_array, H-z_array, np.log(ke), vmin=-14, vmax = -2, cmap='hot_r')
    ax[0,0].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[0,0].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[0,0])
    ax[0,0].set(xlabel='$x$ (m)', ylabel='$z$ (m)', title='KE $(\mathrm{m}^2/\mathrm{s}^2)$: '+ str(int(time_array[i]))+' mins')
    
    im = ax[0,1].pcolor(x_array, H-z_array, np.log(omega), vmin=-14, vmax = -2, cmap='hot_r')
    ax[0,1].plot(x_array, H-bh_array, '-k', lw=2.0)
    ax[0,1].plot(x_array, z0+0*x_array, ':k', lw=2.0)
    plt.colorbar(im, ax=ax[0,1])
    ax[0,1].yaxis.set_ticklabels([])
    ax[0,1].set(xlabel='$x$ (m)', title='$\Omega\, (1/\mathrm{s}^2)$ : '+ str(int(time_array[i]))+' mins')

    im = ax[1,0].pcolor(N2_i, KE_j, np.log(I_N_KE/tot_vol).T, vmin=-14, vmax = -2, cmap='hot_r')
    plt.colorbar(im, ax=ax[1,0])
    ax[1,0].set(xlabel='$N^2\, (1/\mathrm{s}^2)$', ylabel='KE $(\mathrm{m}^2/\mathrm{s}^2)$')
    ax[1,0].set_ylim([0, 0.01])
    
    im = ax[1,1].pcolor(N2_i, omega_j, np.log(I_N_omega/tot_vol).T, vmin=-14, vmax = -2, cmap='hot_r')
    plt.colorbar(im, ax=ax[1,1])
    ax[1,1].set(xlabel='$N^2\, (1/\mathrm{s}^2)$', ylabel='$\Omega\, (1/\mathrm{s}^2)$')
    ax[1,1].set_ylim([0, 0.01])

    filename = foldername+'KE_omega_USP_'+str(int(time_array[i]))+'.png'
    fig.savefig(filename)
    plt.close('all')

    return np.nansum(cell_vol*ke)

def key_listname(s):
    nums = re.findall(r'\d+', s)
    return int(nums[-1])

def make_movie(foldername, filekey, txtname, moviename, ind_text=28):
    filenames = foldername + filekey + '*.png'
    imfiles = sorted(glob.glob(filenames))

    txtfile = foldername + txtname + '.txt'
    sorted_list = sorted(imfiles, key=key_listname)
    with open(txtfile, 'w') as f:
        for ind in range(len(sorted_list)):
            f.write("file " + sorted_list[ind][ind_text:]+"\n")

    moviefile = foldername + moviename
    (
    ffmpeg
    .input(txtfile, r='5', f='concat', safe='0')
    .output(moviefile, pix_fmt='yuv420p', vcodec='libx264')
    .run()
    )
