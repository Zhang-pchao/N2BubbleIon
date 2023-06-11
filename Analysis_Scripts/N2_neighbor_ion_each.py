#!/usr/bin/env python
# coding: utf-8

# In[52]:


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
each_step   = 500
bond_cutoff = 5.0
ReSmth      = 3.0
box_center  = [0.0,0.0,0.0]
####################change above####################

trj_data    = os.path.join(geo_path,trj_data)
n2__data    = os.path.join(geo_path,n2__data)
trj_lmps    = os.path.join(trj_path,trj_lmps)
n2__lmps    = os.path.join(trj_path,n2__lmps) # 'id type cluster cluster cluster'
Ion_file    = os.path.join(trj_path,'FindIon',Ion_file)
Ion_txt     = np.loadtxt(Ion_file)
Re_txt      = np.loadtxt(os.path.join(trj_path,'Equation.txt'),usecols=8)


# In[55]:


u = mda.Universe(trj_data,trj_lmps, atom_style='id type x y z',format='LAMMPSDUMP')


# In[56]:


n = mda.Universe(n2__data,n2__lmps, atom_style='id type x y z',format='LAMMPSDUMP')


# In[57]:


choose_frames_all = range(trj_skip,len(u.trajectory),trj_step)


# In[58]:


def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)


# In[59]:


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


# In[60]:


def N_O_Dict(u,n,choose_frames,Ion_txt,box_center,bond_cutoff):
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
                for j in ion_o_idx:
                    # find O near N within cutoff
                    # j-1 for match index
                    bond_no = u.select_atoms('index {0}'.format(i),'index {0}'.format(int(j-1))).bond.value()
                    if bond_no <= bond_cutoff:
                        n_o_index[len_nc].append(int(j-1)) 
    return n_o_index


# In[61]:


def Re_N_O_Dict(u,n,choose_frames,Ion_txt,Re,ReSmth,box_center,bond_cutoff):
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
                if len_nc <= Re+ReSmth and len_nc > Re-ReSmth:
                    n_o_index[len_nc] = []
                    for j in ion_o_idx:
                        # find O near N within cutoff
                        # j-1 for match index
                        bond_no = u.select_atoms('index {0}'.format(i),'index {0}'.format(int(j-1))).bond.value()
                        if bond_no <= bond_cutoff:
                            n_o_index[len_nc].append(int(j-1)) 
    return n_o_index


# In[66]:


def min_box_length(u):
    length = []
    for i in range(len(u.trajectory)):
        length.append(u.trajectory[i].dimensions[0]) # iso: dimensions[0] = dimensions[1] = dimensions[2]
    return min(length)


# In[67]:


def atom_counts(bin_centers,choose_frames):
    counts = np.full(bin_centers.size, fill_value=0.0)
    counts_list = counts
    for i in range(len(choose_frames)-1):
        counts_list = np.vstack((counts_list, counts))
    return counts_list


# In[71]:


def no_num_scatter(x,y,idx,savepath='./'):
    fig  = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax   = fig.add_subplot(1, 1, 1)
    ax.scatter(x,y,s=10)
    #ax.legend(fontsize = 12)
    ax.set_xlabel("Distance to the bubble center "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Number of Ions near the Nitrogen", fontsize = 12)
    fig.savefig(os.path.join(savepath, "no_num_scatter_{0}.png".format(idx)), dpi=600, bbox_inches='tight')


# In[73]:


def no_num_hist2d(x,y,idx,binss=[150,15],savepath='./'):
    fig  = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax   = fig.add_subplot(1, 1, 1)
    h = ax.hist2d(x, y, bins=binss,range=[[-0.5,50.5], [-0.25,7.25]],cmap='Blues') 
    fig.colorbar(h[3],shrink=1.0) #orientation='horizontal',
    ax.set_xlabel("Distance to the bubble center "+r"$\ \rm (\AA)$", fontsize = 12)
    ax.set_ylabel("Number of Ions near the Nitrogen", fontsize = 12)
    fig.savefig(os.path.join(savepath, "no_num_hist2d_{0}.png".format(idx)), dpi=600, bbox_inches='tight')


# In[ ]:


for ii in range(int(len(u.trajectory)/each_step)):
    choose_frames = choose_frames_all[ii*each_step+1:(ii+1)*each_step+1]
    n_o_index = Re_N_O_Dict(u,n,choose_frames,Ion_txt,Re_txt[ii],ReSmth,box_center,bond_cutoff)
    
    x = []
    y = []   
    file = open(os.path.join('./', "no_num_{0}.txt".format(ii)), 'w+')
    file.write('%16s%16s\n' %("NO_Distance[A]","No. Ions near N"))
    for i in n_o_index.keys():
        x.append(i)
        y.append(len(n_o_index[i]))
        file.write('%16.8f%16d\n' %(i,len(n_o_index[i])))
    file.close()
                   
    no_num_scatter(x,y,ii)
    no_num_hist2d(x,y,ii)


# In[62]:


#n_o_ratio,n_o_index = N_OO_Dict(u,n,choose_frames,Ion_txt,box_center,bond_cutoff)
#n_o_index = N_O_Dict(u,n,choose_frames,Ion_txt,box_center,bond_cutoff)
#n_o_index = Re_N_O_Dict(u,n,choose_frames,Ion_txt,Re_txt[-1],ReSmth,box_center,bond_cutoff)


# In[63]:


#for i in n_o_ratio.keys():
#    if n_o_ratio[i][1] != 0:
#        print(n_o_ratio[i])


# In[64]:


#print(n_o_ratio)


# In[65]:


#print(n_o_index)


# In[68]:


#delta_r     = 2
#box_length  = 50 # int(min_box_length(u)/2+1)*1.732
#bin_edges   = np.linspace(0, box_length, int(box_length/delta_r)+1)
#bin_centers = bin_edges[:-1] + delta_r/2


# In[69]:


#counts_list = atom_counts(bin_centers,1)
#ratio_sum   = atom_counts(bin_centers,1)


# In[70]:


#for i in n_o_ratio.keys():
#    oo_ratio = n_o_ratio[i][1]/n_o_ratio[i][2]
#    hist, xxx = np.histogram(i, bins=bin_edges)
#    counts_list += hist
#    for j in range(len(hist)):
#        if hist[j] == 1:
#            ratio_sum[j] += oo_ratio
#for i in range(len(counts_list)):
#    ratio_sum[i] = ratio_sum[i]/counts_list[i]

