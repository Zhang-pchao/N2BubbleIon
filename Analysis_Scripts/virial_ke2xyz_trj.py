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
geo_path    = "../"
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

def virial2xyz_trj(atom_nums, xyz_index, trj_file,trj_path):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    new_trj_file = open(os.path.join(trj_path,"bubble_virial2xyz.lammpstrj"), "w+")
    
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
            new_trj_file.write('%8d%4d%16.6f%16.6f%16.6f\n' %(int(sp[0]),int(sp[1]),float(sp[5]),float(sp[7]),float(sp[9])))
    
    new_trj_file.close()

def ke2xyz_trj(atom_nums, xyz_index, trj_file,trj_path):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    new_trj_file = open(os.path.join(trj_path,"bubble_ke2xyz.lammpstrj"), "w+")
    
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
            new_trj_file.write('%8d%4d%16.6f%16.6f%16.6f\n' %(int(sp[0]),int(sp[1]),float(sp[6]),float(sp[8]),float(sp[10])))
    
    new_trj_file.close()

virial2xyz_trj(atom_nums, xyz_index, trj_file,trj_path)
ke2xyz_trj(atom_nums, xyz_index, trj_file,trj_path)