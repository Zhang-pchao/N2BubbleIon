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
trj_xyz     = "bubble_adjust_center.lammpstrj"
trj_vil     = "bubble_virial2xyz.lammpstrj"
trj_kec     = "bubble_ke2xyz.lammpstrj"
delta_r     = 1 # bin slice (Angstrom)
trj_step    = 1
trj_skip    = 0
each_step   = 500
####################change above###################

data_geo    = os.path.join(geo_path,data_geo)
trj_xyz     = os.path.join(trj_path,trj_xyz)
trj_vil     = os.path.join(trj_path,trj_vil)
trj_kec     = os.path.join(trj_path,trj_kec)

u_x = mda.Universe(data_geo,trj_xyz, atom_style='id type x y z',format='LAMMPSDUMP')
u_v = mda.Universe(data_geo,trj_vil, atom_style='id type x y z',format='LAMMPSDUMP')
u_k = mda.Universe(data_geo,trj_kec, atom_style='id type x y z',format='LAMMPSDUMP')

def min_box_length(u):
    length = []
    for i in range(len(u.trajectory)):
        length.append(u.trajectory[i].dimensions[0]) # iso: dimensions[0] = dimensions[1] = dimensions[2]
    return min(length)

def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)

def press_fig(save_path,stress_dict,virial_dict,kineti_dict,trj_step,delta_r,ii):
    fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    
    ax.plot(stress_dict.keys(), stress_dict.values(), alpha=0.7, lw=3, label='Press')
    ax.plot(virial_dict.keys(), virial_dict.values(), alpha=0.7, lw=3, label='Virial contribution')
    ax.plot(kineti_dict.keys(), kineti_dict.values(), alpha=0.7, lw=3, label='Kinetic contribution')
    ax.legend(fontsize = 12)
        
    ax.set_xlabel("Radial distance from bubble center "+r"$\ \rm (\AA)$",fontsize = 12)
    ax.set_ylabel("Press (MPa)",fontsize = 12)
    fig.savefig(os.path.join(save_path, "Press_step{0}_dr{1}_{2}.png".format(trj_step,delta_r,ii)), dpi=600, bbox_inches='tight')    

choose_frames_all = range(trj_skip,len(u_x.trajectory),trj_step)

box_length  = int(min_box_length(u_x))+4 # get minimum NPT box length
box_center  = [0,0,0] # iso
bin_edges   = np.linspace(0, box_length/2, int(box_length/2/delta_r)+1)
bin_centers = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
    
v_list = []
for i in range(len(bin_centers)):
    v = 4/3*np.pi*(np.power(bin_edges[i+1],3)-np.power(bin_edges[i],3)) # A^3
    v_list.append(v)

for ii in range(int(len(u_x.trajectory)/each_step)):
    choose_frames = choose_frames_all[ii*each_step+1:(ii+1)*each_step+1]
    
    stress_dict = {}
    virial_dict = {}
    kineti_dict = {}
    for i in bin_centers:
        stress_dict[i] = 0.0
        virial_dict[i] = 0.0
        kineti_dict[i] = 0.0
    
    for i in choose_frames:
        u_x.trajectory[i]   
        u_v.trajectory[i]
        u_k.trajectory[i]
        for k in range(len(u_x.atoms)):
            d_cent_atom = get_distance(u_x.atoms[k].position, box_center)
            for j in range(len(bin_centers)):
                if d_cent_atom >= bin_edges[j] and d_cent_atom < bin_edges[j+1]:
                    stress_dict[bin_centers[j]] = stress_dict[bin_centers[j]]+u_v.atoms[k].position[0]+u_v.atoms[k].position[1]+u_v.atoms[k].position[2]-u_k.atoms[k].position[0]-u_k.atoms[k].position[1]-u_k.atoms[k].position[2]
                    virial_dict[bin_centers[j]] = virial_dict[bin_centers[j]]+u_v.atoms[k].position[0]+u_v.atoms[k].position[1]+u_v.atoms[k].position[2]
                    kineti_dict[bin_centers[j]] = kineti_dict[bin_centers[j]]-u_k.atoms[k].position[0]-u_k.atoms[k].position[1]-u_k.atoms[k].position[2]
                    #break
    
    for i in range(len(bin_centers)):
        stress_dict[bin_centers[i]] = stress_dict[bin_centers[i]]/(3*v_list[i])/len(choose_frames)/10 # bar to MPa
        virial_dict[bin_centers[i]] = virial_dict[bin_centers[i]]/(3*v_list[i])/len(choose_frames)/10 # bar to MPa
        kineti_dict[bin_centers[i]] = kineti_dict[bin_centers[i]]/(3*v_list[i])/len(choose_frames)/10 # bar to MPa
    
    Pfile = open(os.path.join(save_path, "Press_step{0}_dr{1}_{2}.txt".format(trj_step,delta_r,ii)), 'w+')
    Pfile.write('# %16s%32s%32s%32s\n' %('Distance(A)','Press(MPa)','Virial(MPa)','Kinetic(MPa)'))
    for i in range(len(bin_centers)):
        Pfile.write('%16.4f%32.8f%32.8f%32.8f\n' %(bin_centers[i],stress_dict[bin_centers[i]],virial_dict[bin_centers[i]],kineti_dict[bin_centers[i]]))
    Pfile.close()
    
    press_fig(save_path,stress_dict,virial_dict,kineti_dict,trj_step,delta_r,ii)