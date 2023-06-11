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
trj_path    = "../"
trj_name    = "bubble_adjust_center.lammpstrj"
####################change above####################

trj_file    = os.path.join(trj_path,trj_name)

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
        if "ITEM: ATOMS id type" in line:
            xyz_index.append(j)
    
    return atom_nums, xyz_index

atom_nums, xyz_index = trj_info(trj_file)

def onlyxyz_trj(atom_nums, xyz_index, trj_file):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    new_trj_file = open(os.path.join(trj_path,"bubble_only_xyz.lammpstrj"), "w+")
    
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
            new_trj_file.write('%8d%4d%16.6f%16.6f%16.6f\n' %(int(sp[0]),int(sp[1]),float(sp[2]),float(sp[3]),float(sp[4])))
    
    new_trj_file.close()

def shift_xyz(x,box_upper,box_lower,box_length):
    if x > box_upper:
        x = x - box_length
    elif x < box_lower:
        x = x + box_length
    else:
        x = x
    return x

def adjust_trj(atom_nums, xyz_index, trj_path, trj_file, COM_idx=0):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    new_trj_file = open(os.path.join(trj_path,"bubble_adjust_center_real.lammpstrj"), "w+")
    
    k = COM_idx
    for i in xyz_index:
        com_file = os.path.join(trj_path, "COM","COM.{0}".format(k))
        k += 1
        center = np.loadtxt(com_file)
        shift_center = [center[2],center[3],center[4]]
        
        box_lower = [float(lines[i-4].split()[0]),float(lines[i-3].split()[0]),float(lines[i-2].split()[0])]
        box_upper = [float(lines[i-4].split()[1]),float(lines[i-3].split()[1]),float(lines[i-2].split()[1])]
        box_center = [(box_lower[0]+box_upper[0])/2,(box_lower[1]+box_upper[1])/2,(box_lower[2]+box_upper[2])/2]
        box_length = [box_upper[0]-box_lower[0],box_upper[1]-box_lower[1],box_upper[2]-box_lower[2]]
        
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
            
            x = float(sp[2])+box_center[0]-shift_center[0]
            y = float(sp[3])+box_center[1]-shift_center[1]
            z = float(sp[4])+box_center[2]-shift_center[2]                                     
            x = shift_xyz(x,box_upper[0],box_lower[0],box_length[0])
            y = shift_xyz(y,box_upper[1],box_lower[1],box_length[1])
            z = shift_xyz(z,box_upper[2],box_lower[2],box_length[2])                                     
                
            new_trj_file.write('%8d%4d%16.6f%16.6f%16.6f\n' %(int(sp[0]),int(sp[1]),x,y,z))
    
    new_trj_file.close()

adjust_trj(atom_nums, xyz_index[:], trj_path, trj_file, 0)