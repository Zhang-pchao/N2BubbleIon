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
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/BubbleIon/datafile"
trj_path    = "../"
save_path   = "./"
data_geo    = "D30L55N2_20H_atomic.data" #"D30L55N2_20OHcenter_atomic.data"
trj_name    = "bubble_adjust_center.lammpstrj" 
o_xyz       = "h3o_o_xyz_dump1_cut1.25_bubble_adjust_center.txt"
#o_xyz       = "oh_o_xyz_dump1_cut1.275_bubble_adjust_center.txt"
trj_skip    = 0
trj_step    = 10
####################change above####################

data_geo    = os.path.join(geo_path,data_geo)
oid_file    = os.path.join(trj_path,'FindIon',o_xyz)
trj_file    = os.path.join(trj_path,trj_name)

time_start  = 12500 #0 5000
time_step   = 1e-3 # ns
oid_txt     = np.loadtxt(oid_file)

u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
choose_frames = range(trj_skip,len(u.trajectory),trj_step)

# set hbonds
def get_hbonds(u,index,scheme):
    if scheme == "donor":
        hbonds = HydrogenBondAnalysis(
            universe           = u,
            donors_sel         = "index {0}".format(index), # O
            hydrogens_sel      = "type 1", # H
            acceptors_sel      = "type 2", # O
            d_a_cutoff         = 3.5,      # <3.5
            d_h_a_angle_cutoff = 140,      # >140
            update_selections  = False
        )
    else:  # scheme == "acceptor"
        hbonds = HydrogenBondAnalysis(
            universe           = u,
            donors_sel         = "type 2", # O
            hydrogens_sel      = "type 1", # H
            acceptors_sel      = "index {0}".format(index), # O
            d_a_cutoff         = 3.5,      # <3.5
            d_h_a_angle_cutoff = 140,      # >140
            update_selections  = False
        )    
    return hbonds

def get_ion_hb_dict(u,o_id,scheme,frame):    
    hbonds = get_hbonds(u,o_id,scheme)
    hbonds.run(start=frame,stop=frame+1,step=None,verbose=False)
    u.trajectory[frame]
    o_d_or_a = hbonds.results.hbonds.shape[0]    
    return o_d_or_a,hbonds.results.hbonds

def get_neg_hb_dict(u,o_id,scheme,frame):    
    hbonds = get_hbonds(u,o_id,scheme)
    hbonds.run(start=frame,stop=frame+1,step=None,verbose=False)
    u.trajectory[frame]
    o_d_or_a = hbonds.results.hbonds.shape[0]
    return o_d_or_a

def get_fig(ion_donor_dict,ion_accpt_dict,
            ing_donor_dict,ing_accpt_dict,
            save_path,time_start,time_step):
    fig = plt.figure(figsize=(12,6), dpi=150, facecolor='white')    
    
    x = []
    d1 = []
    a1 = []
    d2 = []
    a2 = []
    for i in ion_donor_dict.keys():
        x.append((i+time_start)*time_step)
        d1.append(ion_donor_dict[i])
        a1.append(ion_accpt_dict[i])
        d2.append(ing_donor_dict[i])
        a2.append(ing_accpt_dict[i])
        
    ax = fig.add_subplot(121)    
    ax.plot(x,d1,label='Ion Donor')
    ax.plot(x,a1,label='Ion Acceptor')
    ax.legend(fontsize = 12)
    ax.set_xlabel("Time (ns)",fontsize = 15)
    ax.set_ylabel("No. of Donors or Acceptors",fontsize = 15)
    
    ax2 = fig.add_subplot(122)    
    ax2.plot(x,d2,label='H2O Donor')
    ax2.plot(x,a2,label='H2O Acceptor')
    ax2.legend(fontsize = 12)
    ax2.set_xlabel("Time (ns)",fontsize = 15)
    #ax2.set_ylabel("No. of Donors or Acceptors",fontsize = 15)    
    
    #at = AnchoredText("(a)", prop=dict(size=15), frameon=True, loc='upper right')
    #at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    #ax.add_artist(at) 
    
    fig.savefig(os.path.join(save_path, "HB_DonorAcceptor.png"), dpi=600, bbox_inches='tight')

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path) 

    if not isExists:
        os.makedirs(path) 

def write_dict(path,
               ion_donor_dict,ion_accpt_dict,
               ing_donor_dict,ing_accpt_dict,
               name,time_start,time_step,choose_frames):
    file = open(os.path.join(path, "IonHBnumber_{0}.txt".format(name)), 'w+')
    file.write('# %16s%16s%16s%16s%16s\n' %('Time(ns)','Donors','Acceptors','NEG_Donors','NEG_Acceptors'))
    for i in choose_frames:
        file.write('%16.4f%16.4f%16.4f%16.4f%16.4f\n' %((i+time_start)*time_step,
                                                        ion_donor_dict[i],ion_accpt_dict[i],
                                                        ing_donor_dict[i],ing_accpt_dict[i]))
    file.close()    

ion_donor_dict = {} # find ion as donor
ion_accpt_dict = {}
ing_donor_dict = {} # find neighbor H2O of ion as donor
ing_accpt_dict = {}

#mkdir(save_path)
file = open(os.path.join(save_path, "IonHBnumber.txt"), 'w+')
file.write('# %16s%16s%16s%16s%16s%16s%16s%16s%16s%16s\n' %('Time(ns)','Ions','Ions_d','Ions_a','Ions_dd','Ions_aa',
                                                            'Donors','Acceptors','NEG_Donors','NEG_Acceptors'))
for i in choose_frames:
    ion_donor_dict[i] = 0
    ion_accpt_dict[i] = 0
    ing_donor_dict[i] = 0
    ing_accpt_dict[i] = 0
    for j in range(int(oid_txt[i][0])):
        ion_o_id = int(oid_txt[i][2+4*j]-1)
        ion_o_donors,hbonds_results = get_ion_hb_dict(u,ion_o_id,"donor",   i)
        if ion_o_donors > 0:
            for frame, donor_index, hydrogen_index, acceptor_index, DA_distance, DHA_angle in hbonds_results:
                ion_d_a = acceptor_index.astype(int)
                ion_d_a_d = get_neg_hb_dict(u,ion_d_a,"donor",   i)
                ion_d_a_a = get_neg_hb_dict(u,ion_d_a,"acceptor",i)
                ing_donor_dict[i] += ion_d_a_d
                ing_accpt_dict[i] += ion_d_a_a
        
        ion_o_accpts,hbonds_results = get_ion_hb_dict(u,ion_o_id,"acceptor",i)
        if ion_o_accpts > 0:        
            for frame, donor_index, hydrogen_index, acceptor_index, DA_distance, DHA_angle in hbonds_results:
                ion_a_d = donor_index.astype(int)
                ion_a_d_d = get_neg_hb_dict(u,ion_a_d,"donor",   i)
                ion_a_d_a = get_neg_hb_dict(u,ion_a_d,"acceptor",i) 
                ing_donor_dict[i] += ion_a_d_d
                ing_accpt_dict[i] += ion_a_d_a
        
        ion_donor_dict[i] += ion_o_donors
        ion_accpt_dict[i] += ion_o_accpts
    file.write('%16.4f%16.4f%16.4f%16.4f%16.4f%16.4f' %((i+time_start)*time_step,
                                             int(oid_txt[i][0]),
                                             ion_donor_dict[i],
                                             ion_accpt_dict[i],
                                             ing_donor_dict[i],
                                             ing_accpt_dict[i]
                                            ))    
    ing_donor_dict[i] /= (ion_donor_dict[i]+ion_accpt_dict[i])
    ing_accpt_dict[i] /= (ion_donor_dict[i]+ion_accpt_dict[i])        
    ion_donor_dict[i] /= int(oid_txt[i][0])
    ion_accpt_dict[i] /= int(oid_txt[i][0])
    file.write('%16.4f%16.4f%16.4f%16.4f\n' %(ion_donor_dict[i],
                                              ion_accpt_dict[i],
                                              ing_donor_dict[i],
                                              ing_accpt_dict[i]                       
                                            ))    
file.close()    

get_fig(ion_donor_dict,ion_accpt_dict,
        ing_donor_dict,ing_accpt_dict,
        save_path,time_start,time_step)