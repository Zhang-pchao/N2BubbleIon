#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.offsetbox import AnchoredText
import os
import MDAnalysis as mda

####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/BubbleIon/datafile"
trj_path    = "../"
save_path   = "./"
data_geo    = "D30L55N2_20H_atomic.data"
trj_name    = "bubble_adjust_center.lammpstrj"
delta_r     = 0.2 # bin slice (Angstrom)
trj_step    = 1
trj_skip    = 0
each_step   = 500
####################change above####################
c_mass      = 1.6605390666e-27 #  1/12C kg
o_mass      = 15.999 * c_mass * 1e3 # g
h_mass      = 1.008  * c_mass * 1e3 # g
n_mass      = 14.007 * c_mass * 1e3 # g

data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')

def min_box_length(u):
    length = []
    for i in range(len(u.trajectory)):
        length.append(u.trajectory[i].dimensions[0]) # iso: dimensions[0] = dimensions[1] = dimensions[2]
    return min(length)

def atom_counts(bin_centers,choose_frames):
    counts = np.full(bin_centers.size, fill_value=0.0)
    counts_list = counts
    for i in range(len(choose_frames)-1):
        counts_list = np.vstack((counts_list, counts))
    return counts_list

def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)

def get_se_up_low(avrg,se):
    up  = []
    low = []
    for i in range(len(avrg)):
        up.append(avrg[i]+se[i])
        low.append(avrg[i]-se[i])
    return up, low

def density_fig(bin_centers,
                bubble_density,water_density,nitrogen_density,
                se_bubble_density,se_water_density,se_nitrogen_density,
                trj_step,delta_r,save_path,ii):
    fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    
    ax.plot(bin_centers, bubble_density,   alpha=0.7, lw=3, label='Bubble',   c='mediumaquamarine')
    ax.plot(bin_centers, water_density,    alpha=0.7, lw=3, label='Water',    c='orange')
    ax.plot(bin_centers, nitrogen_density, alpha=0.7, lw=3, label='Nitrogen', c='royalblue')
    bubble_se1,bubble_se2  = get_se_up_low(bubble_density,se_bubble_density)
    ax.fill_between(bin_centers,bubble_se1,bubble_se2,
                    facecolor = 'mediumaquamarine', alpha = 0.2)
    water_se1,water_se2 = get_se_up_low(water_density,se_water_density)
    ax.fill_between(bin_centers,water_se1,water_se2,
                    facecolor = 'orange', alpha = 0.2)
    nitrogen_se1,nitrogen_se2 = get_se_up_low(nitrogen_density,se_nitrogen_density)
    ax.fill_between(bin_centers,nitrogen_se1,nitrogen_se2,
                    facecolor = 'royalblue', alpha = 0.2)
    ax.legend(fontsize = 12)
    
    ax.set_xlabel("Radial distance from bubble center "+r"$\ \rm (\AA)$",fontsize = 12)
    ax.set_ylabel("Density (g/" + '$\mathregular{cm^3)}$',fontsize = 12)
    fig.savefig(os.path.join(save_path, "Density_step{0}_dr{1}_{2}.png".format(trj_step,delta_r,ii)), dpi=600, bbox_inches='tight')    

box_length = min_box_length(u) # get minimum NPT box length
box_center  = [0,0,0] # iso
bin_edges   = np.linspace(0, int(box_length/2), int(int(box_length/2)/delta_r)+1)
bin_centers = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2

v_list = []
for i in range(len(bin_centers)):
    v = 4/3*np.pi*(np.power(bin_edges[i+1],3)-np.power(bin_edges[i],3))*1e-24 # cm^3
    v_list.append(v)
    
choose_frames_all = range(trj_skip,len(u.trajectory),trj_step)

for ii in range(int(len(u.trajectory)/each_step)):
    choose_frames = choose_frames_all[ii*each_step+1:(ii+1)*each_step+1]

    o_counts_list = atom_counts(bin_centers,choose_frames)
    h_counts_list = atom_counts(bin_centers,choose_frames)
    n_counts_list = atom_counts(bin_centers,choose_frames)
    
    for i in range(len(choose_frames)):
        u.trajectory[int(i*trj_step+ii*each_step+1)]        
        for k in range(len(u.atoms)):
            d_cent_atom = get_distance(u.atoms[k].position, box_center)
            hist, xxx = np.histogram(d_cent_atom, bins=bin_edges)
            if  u.atoms[k].type  == '1':            
                h_counts_list[i] += hist
            elif u.atoms[k].type == '3':            
                n_counts_list[i] += hist    
            else: # type =='2'
                o_counts_list[i] += hist
    
    h_density_avg = np.average(h_counts_list, axis=0)
    o_density_avg = np.average(o_counts_list, axis=0)
    n_density_avg = np.average(n_counts_list, axis=0)
    
    h_density_se = np.std(h_counts_list, axis=0, ddof=1)/np.sqrt(len(choose_frames))
    o_density_se = np.std(o_counts_list, axis=0, ddof=1)/np.sqrt(len(choose_frames))
    n_density_se = np.std(n_counts_list, axis=0, ddof=1)/np.sqrt(len(choose_frames))
    
    bubble_density   = []
    water_density    = []
    nitrogen_density = []
    se_bubble_density   = []
    se_water_density    = []
    se_nitrogen_density = []
    for i in range(len(h_density_avg)):
        d_all = (h_density_avg[i]*h_mass+o_density_avg[i]*o_mass+n_density_avg[i]*n_mass)/v_list[i]
        d_h2o = (h_density_avg[i]*h_mass+o_density_avg[i]*o_mass)/v_list[i]
        d_n2  = (n_density_avg[i]*n_mass)/v_list[i]
        se_all = (h_density_se[i]*h_mass+o_density_se[i]*o_mass+n_density_se[i]*n_mass)/v_list[i]
        se_h2o = (h_density_se[i]*h_mass+o_density_se[i]*o_mass)/v_list[i]
        se_n2  = (n_density_se[i]*n_mass)/v_list[i]
        bubble_density.append(d_all)
        water_density.append(d_h2o)
        nitrogen_density.append(d_n2)
        se_bubble_density.append(se_all)
        se_water_density.append(se_h2o)
        se_nitrogen_density.append(se_n2)
    
    Dfile = open(os.path.join(save_path, "Density_step{0}_dr{1}_{2}.txt".format(trj_step,delta_r,ii)), 'w+')
    Dfile.write('#%16s%16s%16s%16s%16s%16s%16s\n' %('Distance(A)','Bubble(g/cm^3)','StandErr_Bub','H2O(g/cm^3)','StandErr_H2O','N2(g/cm^3)','StandErr_N2'))
    for i in range(len(bin_centers)):
        Dfile.write('%16.4f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n' %(bin_centers[i],
                                                                    bubble_density[i],  se_bubble_density[i],
                                                                    water_density[i],   se_water_density[i],
                                                                    nitrogen_density[i],se_nitrogen_density[i]))
    Dfile.close()
    density_fig(bin_centers,
                    bubble_density,water_density,nitrogen_density,
                    se_bubble_density,se_water_density,se_nitrogen_density,
                    trj_step,delta_r,save_path,ii)