#!/usr/bin/env python
# coding: utf-8

# In[91]:


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
save_path   = "./"
trj_lmps    = "bubble_adjust_center.lammpstrj"
n2__lmps    = "bubble_cluster2xyz.lammpstrj"
#Ion_file    = "oh_o_xyz_dump1_cut1.275_bubble_adjust_center.txt"
Ion_file    = "h3o_o_xyz_dump1_cut1.25_bubble_adjust_center.txt"
trj_skip    = 0
trj_step    = 1
each_step   = 500
start_step  = 12500
bond_cutoff = 5.0
box_center  = [0.0,0.0,0.0]
####################change above####################

trj_data    = os.path.join(geo_path,trj_data)
n2__data    = os.path.join(geo_path,n2__data)
trj_lmps    = os.path.join(trj_path,trj_lmps)
n2__lmps    = os.path.join(trj_path,n2__lmps) # 'id type cluster cluster cluster'
Ion_file    = os.path.join(trj_path,'FindIon',Ion_file)
Ion_txt     = np.loadtxt(Ion_file)


# In[93]:


u = mda.Universe(trj_data,trj_lmps, atom_style='id type x y z',format='LAMMPSDUMP')


# In[94]:


n = mda.Universe(n2__data,n2__lmps, atom_style='id type x y z',format='LAMMPSDUMP')


# In[95]:


def get_distance(xyz1, xyz2):
    distance = 0
    for i in range(3):
        distance += np.square((xyz1[i]-xyz2[i]))
    return np.sqrt(distance)


# In[96]:


def Re_O_N_Dict(u,n,choose_frames,Ion_txt,Re,ReSmth,box_center,bond_cutoff):
    n_o_index  = {} # k = distance between N and center of bubble, v = O index list
    for frame_idx in choose_frames:
        u.trajectory[frame_idx]
        n.trajectory[frame_idx]
        
        ion_o_idx = [] # ion O index in a frame
        for k in range(int(Ion_txt[frame_idx][0])):
            ion_o_idx.append(Ion_txt[frame_idx][2+4*k])
        
        for j in ion_o_idx:
            len_oc = get_distance(u.atoms[int(j-1)].position,box_center)
            if len_oc > Re+ReSmth:
                n_o_index[len_oc] = []
                for i in range(len(n.atoms)):
                    bond_no = u.select_atoms('index {0}'.format(i),'index {0}'.format(int(j-1))).bond.value()
                    if bond_no <= bond_cutoff:
                        n_o_index[len_oc].append(i) 
    print(n_o_index)
    return n_o_index


# In[97]:


def Re_O_N_Dict2(u,n,choose_frames,Ion_txt,bond_cutoff):
    n_o_num2  =  0
    #print(choose_frames,len(choose_frames))
    for frame_idx in choose_frames:
        u.trajectory[frame_idx]
        n.trajectory[frame_idx]
        
        ion_o_idx = [] # ion O index in a frame
        n_o_num1  =  0
        for k in range(int(Ion_txt[frame_idx][0])):
            ion_o_idx.append(Ion_txt[frame_idx][2+4*k])
        
        for j in ion_o_idx:
            for i in range(len(n.atoms)):
                bond_no = u.select_atoms('index {0}'.format(i),'index {0}'.format(int(j-1))).bond.value()
                if bond_no <= bond_cutoff:
                    n_o_num1 += 1
        n_o_num1 /= len(ion_o_idx)
        n_o_num2 += n_o_num1
    return n_o_num2/len(choose_frames)


# In[98]:


def no_num_fig(x,y,savepath='./'):
    fig  = plt.figure(figsize=(8,8), dpi=150, facecolor='white')
    ax   = fig.add_subplot(111)
    ax.plot(x,y)
    #ax.legend(fontsize = 12)
    ax.set_xlabel("Time (ns)", fontsize = 12)
    ax.set_ylabel("No. of N2 near Ions", fontsize = 12)
    fig.savefig(os.path.join(savepath, "no_num.png"), dpi=600, bbox_inches='tight')


# In[99]:


#x_time = []
#yN2num = []  
#file = open(os.path.join(trj_path, "Ion_neighor_N2_num.txt"), 'w+')
#file.write('%16s%16s\n' %("time[ns]","N2 near per Ion"))
#
#for ii in range(int(len(u.trajectory)/each_step)):
#    choose_frames = choose_frames_all[ii*each_step+1:(ii+1)*each_step+1]
#    x_time.append((ii+0.5)*each_step/1000+start_step/1000)
#    n_o_index = Re_O_N_Dict(u,n,choose_frames,Ion_txt,Re_txt[ii],ReSmth,box_center,bond_cutoff)
#    
#    total_n_num = 0
#    ions_num = 0
#    for i in n_o_index.keys():
#        total_n_num += len(n_o_index[i])
#        ions_num += 1
#    if ions_num == 0:
#        yN2num.append(0)
#    else:
#        yN2num.append(total_n_num/2/ions_num) # for N 2
#    file.write('%16.8f%16.4f\n' %(x_time[-1],yN2num[-1]))
#               
#file.close()


# In[100]:


choose_frames_all = range(trj_skip,len(u.trajectory),trj_step)


# In[101]:


x_time = []
yN2num = []  
file = open(os.path.join(save_path, "Ion_neighor_N2_num.txt"), 'w+')
file.write('%16s%16s\n' %("time[ns]","N2 near per Ion"))

for ii in range(int(len(u.trajectory)/each_step)):
    choose_frames = choose_frames_all[int(ii*each_step/trj_step+1):int((ii+1)*each_step/trj_step+1)]
    #print(choose_frames)
    x_time.append((ii+0.5)*each_step/1000+start_step/1000)
    n_o_num = Re_O_N_Dict2(u,n,choose_frames,Ion_txt,bond_cutoff)
    yN2num.append(n_o_num)
    file.write('%16.8f%16.4f\n' %(x_time[-1],yN2num[-1]))
               
file.close()


# In[103]:


no_num_fig(x_time,yN2num,save_path)


# In[ ]:




