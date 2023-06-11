#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.offsetbox import AnchoredText
import os
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


# In[5]:


####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/BubbleIon/datafile"
trj_path    = "../"
save_path   = "./"
data_geo    = "D30L55N2_20H_atomic.data"
trj_name    = "bubble_adjust_center.lammpstrj"
delta_r     = 0.2 # bin slice (Angstrom)
trj_step    = 10
trj_skip    = 0
each_step   = 500
box_center  = [0,0,0] # iso
angle       = 140
####################change above####################

data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)


# In[6]:


u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')


# In[7]:


def min_box_length(u):
    length = []
    for i in range(len(u.trajectory)):
        length.append(u.trajectory[i].dimensions[0]) # iso: dimensions[0] = dimensions[1] = dimensions[2]
    return min(length)


# In[8]:


def atom_counts(bin_centers,choose_frames):
    counts = np.full(bin_centers.size, fill_value=0.0)
    counts_list = counts
    for i in range(len(choose_frames)-1):
        counts_list = np.vstack((counts_list, counts))
    return counts_list


# In[9]:


def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)


# In[10]:


def HB_fig(bin_centers,counts_avg,counts_up,counts_down,trj_step,delta_r,save_path,angle,idx):
    fig = plt.figure(figsize=(6.5,6.5), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    
    ax.plot(bin_centers, counts_avg, lw=3, label='Average value', c='firebrick')
    ax.fill_between(bin_centers,counts_up,counts_down,facecolor = 'gold',alpha = 0.6)
    ax.legend(fontsize = 12)
    
    ax.set_xlabel("Radial distance from bubble center "+r"$\ \rm (\AA)$",fontsize = 12)
    ax.set_ylabel("Number of hydrogen bonds per H"+'$\mathregular{_2O}$',fontsize = 12)
    fig.savefig(os.path.join(save_path, "HBperH2O_step{0}_dr{1}_angle{2}_{3}.png".format(trj_step,delta_r,angle,idx)), dpi=150, bbox_inches='tight')    


# In[11]:


box_length  = 28
bin_edges   = np.linspace(0, box_length, int(box_length/delta_r)+1)
bin_centers = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2


# In[12]:


# set hbonds
hbonds = HydrogenBondAnalysis(
    universe           = u,
    donors_sel         = "type 2", # O
    hydrogens_sel      = "type 1", # H
    acceptors_sel      = "type 2", # O
    d_a_cutoff         = 3.5,      # <3.5
    d_h_a_angle_cutoff = angle,      # >140
    update_selections  = False
)

#hbonds.run(
#    start   = trj_skip,                   # skip 10% trajectorys, default:None
#    stop    = None,
#    step    = trj_step,
#    verbose = True
#)

choose_frames_all = range(trj_skip,len(u.trajectory),trj_step)


# In[ ]:


for ii in range(int(len(u.trajectory)/each_step)):
    choose_frames = choose_frames_all[int(ii*each_step/trj_step+1):int((ii+1)*each_step/trj_step+1)]

    hbonds.run(
        start   = int(ii*each_step/trj_step+1),
        stop    = int((ii+1)*each_step/trj_step+1),
        step    = trj_step,
        verbose = True
    )
    
    HB_counts_list   = atom_counts(bin_centers,choose_frames)
    O_counts_list    = atom_counts(bin_centers,choose_frames)
    perO_counts_list = atom_counts(bin_centers,choose_frames)   

    for frame, donor_ix, *_ in hbonds.results.hbonds:
        u.trajectory[frame.astype(int)]
        donor = u.atoms[donor_ix.astype(int)]
        hist, *_ = np.histogram(get_distance(donor.position, box_center), bins=bin_edges)
        # multiply by two as each hydrogen bond involves two water molecules  
        idx = int((frame-int(ii*each_step/trj_step+1))/trj_step)
        HB_counts_list[idx] += hist * 2
    
    for i in hbonds.frames:
        u.trajectory[i]
        idx = int((i-int(ii*each_step/trj_step+1))/trj_step)
        for k in range(len(u.atoms)):
            if u.atoms[k].type =='2': # O
                hist2, xxx = np.histogram(get_distance(u.atoms[k].position, box_center), bins=bin_edges)                
                O_counts_list[idx] += hist2
                
        for l in range(len(perO_counts_list[idx])):
            if O_counts_list[idx][l] == 0:
                perO_counts_list[idx][l] = 0
            else:
                perO_counts_list[idx][l]=HB_counts_list[idx][l]/O_counts_list[idx][l]
    
    counts_avg = np.average(perO_counts_list, axis=0)
    counts_se = np.std(perO_counts_list, axis=0, ddof=1)/np.sqrt(len(choose_frames))
    counts_up   = counts_avg + counts_se
    counts_down = counts_avg - counts_se
    
    HBfile = open(os.path.join(save_path, "HBperH2O_step{0}_dr{1}_angle{2}_{3}.txt".format(trj_step,delta_r,angle,ii)), 'w+')
    HBfile.write('# %16s%16s%16s\n' %('Distance(A)','HBnum','StandErr'))
    for i in range(len(bin_centers)):
        HBfile.write('%16.4f%16.4f%16.4f\n' %(bin_centers[i],counts_avg[i],counts_se[i]))
    HBfile.close()
    
    HB_fig(bin_centers,counts_avg,counts_up,counts_down,trj_step,delta_r,save_path,angle,ii)


# In[ ]:




