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
trj_data    = "D30L55N2_20H_atomic.data"
n2__data    = "D30L55N2_20H_cluster.data"
trj_path    = "../"
trj_lmps    = "bubble_adjust_center.lammpstrj"
n2__lmps    = "bubble_cluster2xyz.lammpstrj"
#Ion_file    = "oh_o_xyz_dump1_cut1.275_bubble_adjust_center.txt"
Ion_file    = "h3o_o_xyz_dump1_cut1.25_bubble_adjust_center.txt"
trj_skip    = 0
trj_step    = 1
bond_cutoff = 10
box_center  = [0.0,0.0,0.0]
####################change above####################

trj_data    = os.path.join(geo_path,trj_data)
n2__data    = os.path.join(geo_path,n2__data)
trj_lmps    = os.path.join(trj_path,trj_lmps)
n2__lmps    = os.path.join(trj_path,n2__lmps) # 'id type cluster cluster cluster'
Ion_file    = os.path.join(trj_path,'FindIon',Ion_file)
Ion_txt     = np.loadtxt(Ion_file)

u = mda.Universe(trj_data,trj_lmps, atom_style='id type x y z',format='LAMMPSDUMP')
n = mda.Universe(n2__data,n2__lmps, atom_style='id type x y z',format='LAMMPSDUMP')
choose_frames = range(trj_skip,len(u.trajectory),trj_step)

def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)

def N_OO_Dict(u,n,choose_frames,Ion_txt,box_center,bond_cutoff):
    n_o_index  = {} # k = distance between N and center of bubble, v = O index list
    n_o_ratio  = {} # k = distance between N and center of bubble, v = [frame index, ion O num, water O num]
    for frame_idx in choose_frames:
        u.trajectory[frame_idx]
        n.trajectory[frame_idx]
        
        ion_o_idx = [] # ion O index in a frame
        for k in range(int(Ion_txt[frame_idx][0])):
            ion_o_idx.append(Ion_txt[frame_idx][2+4*k])
        
        for i in range(len(n.atoms)):
            if int(n.atoms[i].position[0]) != 1: # find N out of bubble
                len_nc = get_distance(u.atoms[i].position,box_center)
                n_o_index[len_nc] = []
                for j in range(len(u.atoms)):
                    # find O near N within cutoff
                    if u.atoms[j].type == '2':
                        bond_no = u.select_atoms('index {0}'.format(i),'index {0}'.format(j)).bond.value()
                        if bond_no <= bond_cutoff:
                            n_o_index[len_nc].append(j+1) # +1 for match index in findion_file
                n_o_ratio[len_nc] = [frame_idx,0,len(n_o_index[len_nc])]          
                for m in ion_o_idx:
                    if m in n_o_index[len_nc]:
                        n_o_ratio[len_nc][1] += 1
    return n_o_ratio,n_o_index

def N_O_Dict(u,n,choose_frames,Ion_txt,box_center,bond_cutoff):
    n_in__o_index  = {} # k = distance between N and center of bubble, v = O index list
    n_out_o_index  = {} # k = distance between N and center of bubble, v = O index list
    for frame_idx in choose_frames:
        u.trajectory[frame_idx]
        n.trajectory[frame_idx]
        
        ion_o_idx = [] # ion O index in a frame
        for k in range(int(Ion_txt[frame_idx][0])):
            ion_o_idx.append(Ion_txt[frame_idx][2+4*k])
        
        for i in range(len(n.atoms)):
            if int(n.atoms[i].position[0]) != 1: # find N out of bubble
                len_nc = get_distance(u.atoms[i].position,box_center)
                n_out_o_index[len_nc] = []
                for j in ion_o_idx:
                    # find O near N within cutoff
                    # j-1 for match index
                    bond_no = u.select_atoms('index {0}'.format(i),'index {0}'.format(int(j-1))).bond.value()
                    if bond_no <= bond_cutoff:
                        n_out_o_index[len_nc].append(int(j-1))
            else: # find N in bubble == 1
                len_nc = get_distance(u.atoms[i].position,box_center)
                n_in__o_index[len_nc] = []
                for j in ion_o_idx:
                    # find O near N within cutoff
                    # j-1 for match index
                    bond_no = u.select_atoms('index {0}'.format(i),'index {0}'.format(int(j-1))).bond.value()
                    if bond_no <= bond_cutoff:
                        n_in__o_index[len_nc].append(int(j-1))                                     
    return n_in__o_index,n_out_o_index

#n_o_ratio,n_o_index = N_OO_Dict(u,n,choose_frames,Ion_txt,box_center,bond_cutoff)
n_in__o_index,n_out_o_index = N_O_Dict(u,n,choose_frames,Ion_txt,box_center,bond_cutoff)

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

def no_num_scatter(Dict,bond_cutoff,scheme,savepath='./'):
    fig  = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax   = fig.add_subplot(1, 1, 1)
    x = []
    y = []
    for i in Dict.keys():
        x.append(i)
        #y.append(Dict[i][1])
        y.append(len(Dict[i]))
    ax.scatter(x,y,s=10)
    #ax.legend(fontsize = 12)
    ax.set_xlabel("Distance to the bubble center "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Number of Ions near the Nitrogen", fontsize = 12)
    fig.savefig(os.path.join(savepath, "no_num_scatter_cutoff{0}_{1}.png".format(bond_cutoff,scheme)), dpi=600, bbox_inches='tight')

def no_num_hist2d(Dict,bond_cutoff,scheme,binss=[150,15],savepath='./'):
    fig  = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax   = fig.add_subplot(1, 1, 1)
    x = []
    y = []
    for i in Dict.keys():
        x.append(i)
        #y.append(Dict[i][1])
        y.append(len(Dict[i]))
    h = ax.hist2d(x, y, bins=binss,range=[[-0.5,50.5], [-0.25,7.25]],cmap='Blues') 
    fig.colorbar(h[3],shrink=1.0) #orientation='horizontal',
    ax.set_xlabel("Distance to the bubble center "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Number of Ions near the Nitrogen", fontsize = 12)
    fig.savefig(os.path.join(savepath, "no_num_hist2d_cutoff{0}_{1}.png".format(bond_cutoff,scheme)), dpi=600, bbox_inches='tight')

no_num_scatter(n_in__o_index,bond_cutoff,"in")
no_num_hist2d( n_in__o_index,bond_cutoff,"in")
no_num_scatter(n_out_o_index,bond_cutoff,"out")
no_num_hist2d( n_out_o_index,bond_cutoff,"out")