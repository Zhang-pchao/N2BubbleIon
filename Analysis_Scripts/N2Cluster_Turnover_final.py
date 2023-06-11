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
data_geo    = "D30L55N2_20H_cluster.data"
trj_name    = "bubble_cluster2index.lammpstrj"
cluster     = "bubble_cluster2xyz.lammpstrj"
trj_skip    = 0
trj_step    = 1
start_time  = 0.0
####################change above####################

data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)
clr_file    = os.path.join(trj_path,cluster) # 'id type cluster cluster cluster'

def trj_info(trj_file):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    atom_index = 0
    i = 0
    for line in lines:
        i += 1
        if "ITEM: NUMBER OF ATOMS" in line:
            atom_index = i
            break
    atom_nums = lines[atom_index].split()[0]
    atom_nums = int(atom_nums)     
    
    xyz_index = []
    j = 0
    for line in lines:
        j += 1
        if "ITEM: ATOMS id" in line:
            xyz_index.append(j)
    
    return atom_nums, xyz_index

def cluster2xyz_trj(atom_nums, xyz_index, trj_file):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    new_trj_file = open(os.path.join(trj_path,"bubble_cluster2xyz.lammpstrj"), "w+")
    
    for i in xyz_index:
        new_trj_file.write(lines[i-9])
        new_trj_file.write(lines[i-8])
        new_trj_file.write(lines[i-7])
        new_trj_file.write(lines[i-6])
        new_trj_file.write(lines[i-5])
        new_trj_file.write(lines[i-4])
        new_trj_file.write(lines[i-3])
        new_trj_file.write(lines[i-2])
        new_trj_file.write("ITEM: ATOMS id type x y z\n")
        for j in range(atom_nums):
            sp = lines[i+j].split()
            new_trj_file.write('%8d%4d%8d%8d%8d\n' %(int(sp[0]),3,int(sp[1]),int(sp[1]),int(sp[1])))
    
    new_trj_file.close()

atom_nums, xyz_index = trj_info(trj_file)

#cluster2xyz_trj(atom_nums, xyz_index[9001:10001], trj_file)
cluster2xyz_trj(atom_nums, xyz_index, trj_file)

u = mda.Universe(data_geo,clr_file, atom_style='id type x y z',format='LAMMPSDUMP')

choose_frames = range(trj_skip,len(u.trajectory),trj_step)

def in_out_cluster(idx,u):
    u.trajectory[idx]
    
    time = u.trajectory.time
    cluster__in = []
    cluster_out = []
    N2_cluster_num = 0
    for i in range(len(u.atoms)):
        if int(u.atoms[i].position[0]) == 1:
            cluster__in.append(u.atoms[i].index)
            N2_cluster_num += 1
        else:
            cluster_out.append(u.atoms[i].index)
    return cluster__in,cluster_out,time,N2_cluster_num

def find_in_out(cluster1,cluster2):
    num = 0
    for i in cluster1:
        if i in cluster2:
            num += 1
    return num/2 # N2

def sum_N2num(list1):
    total = 0
    ele = 0
    while(ele < len(list1)):
        total = total + list1[ele]
        ele += 1
    return total

in__list = []
out_list = []
time     = []
N2_cluster_num = []
for i in choose_frames[:-1]:
    
    cluster__in1,cluster_out1,time1,N2_cluster_num1 = in_out_cluster(i,u)
    cluster__in2,cluster_out2,time2,N2_cluster_num2 = in_out_cluster(i+1,u)
    in__list.append(find_in_out(cluster_out1,cluster__in2))
    out_list.append(find_in_out(cluster__in1,cluster_out2))
    N2_cluster_num.append(N2_cluster_num1/2)
    time.append((time2+1)/1000)

all_in__list = []
all_out_list = []
out_minus_in = []
for i in range(len(in__list)):
    all_in__list.append(sum_N2num(in__list[:i+1]))
for i in range(len(out_list)):
    all_out_list.append(sum_N2num(out_list[:i+1]))
for i in range(len(in__list)):
    out_minus_in.append(all_out_list[i]-all_in__list[i])

def in_out_fig(time,in__list,out_list,all_in__list,all_out_list,out_minus_in,N2_cluster_num,trj_path,start_time):
    fig  = plt.figure(figsize=(13,9), dpi=150, facecolor='white')
    
    ax   = fig.add_subplot(411)
    ax.plot(time,out_list,label="Diffuse from bubble into water")
    ax.plot(time,in__list,label="Diffuse from water into bubble")
    ax.legend(fontsize = 10)
    #ax.set_xlabel("Time (ns)", fontsize = 12)
    ax.set_ylabel("The number of N2", fontsize = 10)
       
    ax   = fig.add_subplot(412)
    ax.plot(time,all_out_list,label="Diffuse from bubble into water")
    ax.plot(time,all_in__list,label="Diffuse from water into bubble")
    ax.legend(fontsize = 10)
    #ax.set_xlabel("Time (ns)", fontsize = 12)
    ax.set_ylabel("Cumulative no. of N2", fontsize = 10)
    
    ax   = fig.add_subplot(413)
    ax.plot(time,out_minus_in,label="[From bubble into water] minus [From water into bubble]")
    ax.legend(fontsize = 10)
    #ax.set_xlabel("Time (ns)", fontsize = 12)
    ax.set_ylabel("Cumulative no. of N2", fontsize = 10) 
    
    ax   = fig.add_subplot(414)
    ax.plot(time,N2_cluster_num)
    #ax.legend(fontsize = 10)
    ax.set_xlabel("Time (ns)", fontsize = 10)
    ax.set_ylabel("no. of N2 in bubble", fontsize = 10)   
    
    file = open(os.path.join(trj_path, "N2Cluster.txt"), 'w+')
    file.write('#%16s%16s%16s%16s%16s%16s%16s\n' %('Time(ns)','each_in','each_out','cumu_in','cumu_out','out_minus_in','N2_num_in_bubble'))
    for i in range(len(time)):
        file.write('%16.4f%16d%16d%16d%16d%16d%16d\n' %(time[i]+start_time,in__list[i],out_list[i],all_in__list[i],all_out_list[i],out_minus_in[i],N2_cluster_num[i]))
    file.close()    
    
    fig.savefig(os.path.join(trj_path, "N2Cluster.png"), dpi=600, bbox_inches='tight')

in_out_fig(time,in__list,out_list,all_in__list,all_out_list,out_minus_in,N2_cluster_num,save_path,start_time)