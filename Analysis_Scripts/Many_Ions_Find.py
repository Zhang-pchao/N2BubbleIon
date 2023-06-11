#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib

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
            
    #box_long  = []
    box_ab  = []
    k = 0
    for line in lines:
        k += 1
        if "ITEM: BOX BOUNDS pp pp pp" in line:
            a1 = float(lines[k].split()[0])
            b1 = float(lines[k].split()[1])
            a2 = float(lines[k+1].split()[0])
            b2 = float(lines[k+1].split()[1])
            a3 = float(lines[k+2].split()[0])
            b3 = float(lines[k+2].split()[1])
            #box_long.append([b1-a1, b2-a2, b3-a3])
            box_ab.append([b1, a1, b2, a2, b3, a3])
    
    return atom_nums, xyz_index, box_ab

def write_xyz(atom_nums, xyz_index, box_ab, trj_file, cutoff, scheme,system):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    xyz = []
    pbc_xyz = []
    
    for i in range(atom_nums):
        e_x = lines[xyz_index + i].split()
        atom_index = int(e_x[0])
        atom_type  = int(e_x[1])
        x = float(e_x[2])
        y = float(e_x[3])
        z = float(e_x[4])
        xyz.append([atom_index, atom_type, x, y, z])
    #xyz = np.array(xyz, dtype=np.float16)

    if system == 'sphere':
        pbc_xyz = xyz
    else: # system = bulk, slab
        for j in range(atom_nums):
            e_x = lines[xyz_index + j].split()
            atom_index = int(e_x[0])
            atom_type  = int(e_x[1])
            x1 = xyz[j][2] - (box_ab[0]-box_ab[1])
            x2 = xyz[j][2]
            x3 = xyz[j][2] + (box_ab[0]-box_ab[1])
            y1 = xyz[j][3] - (box_ab[2]-box_ab[3])
            y2 = xyz[j][3]
            y3 = xyz[j][3] + (box_ab[2]-box_ab[3])
            xxx = [x1, x2, x3]
            yyy = [y1, y2, y3]
            
            if system == 'slab':  # for slab at z-axis
                z2 = xyz[j][4]
            else: # bulk
                z1 = xyz[j][4] - (box_ab[4]-box_ab[5])
                z2 = xyz[j][4]
                z3 = xyz[j][4] + (box_ab[4]-box_ab[5])
                zzz = [z1, z2, z3]
            
            for _x in xxx:
                if _x > (box_ab[1]-cutoff-0.15) and _x < (box_ab[0]+cutoff+0.15):
                    for _y in yyy:
                        if _y > (box_ab[3]-cutoff-0.15) and _y < (box_ab[2]+cutoff+0.15):
                            if system == 'slab':  # for slab at z-axis
                                pbc_xyz.append([atom_index, atom_type, _x, _y, z2])
                            else:
                                for _z in zzz:
                                    if _z > (box_ab[5]-cutoff-0.15) and _z < (box_ab[4]+cutoff+0.15):
                                        pbc_xyz.append([atom_index, atom_type, _x, _y, _z])
        #pbc_xyz = np.array(pbc_xyz, dtype=np.float16)
    
    o_index_dict = {}
    o_coord_dict = {}
    h_coord_dict = {}
    oh_label  = False
    h3o_label = False
    o_list = []
        
    for j in range(len(xyz)):        
        if xyz[j][1] == 2.0: # O:2
            o_list = []
            h_xyz  = []
            for i in range(len(pbc_xyz)):
                if pbc_xyz[i][1] == 1.0: # H:1
                    if pbc_xyz[i][2] <= (xyz[j][2]+cutoff) and pbc_xyz[i][2] >= (xyz[j][2]-cutoff):
                        if pbc_xyz[i][3] <= (xyz[j][3]+cutoff) and pbc_xyz[i][3] >= (xyz[j][3]-cutoff):
                            if pbc_xyz[i][4] <= (xyz[j][4]+cutoff) and pbc_xyz[i][4] >= (xyz[j][4]-cutoff):
                                if (np.square((pbc_xyz[i][2]-xyz[j][2]))+np.square((pbc_xyz[i][3]-xyz[j][3]))+np.square((pbc_xyz[i][4]-xyz[j][4]))) <= np.square(cutoff):   
                                    o_list.append(int(pbc_xyz[i][0]))
                                    h_xyz.append([pbc_xyz[i][2],pbc_xyz[i][3],pbc_xyz[i][4]])
        
            #if pbc_xyz[i][1] == 1.0 and xyz[j][1] == 2.0:
            if scheme == 'h3o':
                if len(o_list) == 3:
                    h3o_label = True
                    o_index_dict[int(xyz[j][0])] = o_list                    
                    o_coord_dict[int(xyz[j][0])] = [xyz[j][2],xyz[j][3],xyz[j][4]]
                    h_coord_dict[int(xyz[j][0])] = h_xyz
            else:# scheme == 'oh':
                if len(o_list) == 1:
                    oh_label  = True
                    o_index_dict[int(xyz[j][0])] = o_list
                    o_coord_dict[int(xyz[j][0])] = [xyz[j][2],xyz[j][3],xyz[j][4]]
                    h_coord_dict[int(xyz[j][0])] = h_xyz
        else:
            o_list = []
    
    MDstep = -1   
    if oh_label or h3o_label:
        MDstep = int(lines[xyz_index - 8])
        #print('MD step of containing ions: ', MDstep)
    else: # not find ions
        if scheme == 'h3o':
            o_index_dict[-1] = []
            o_coord_dict[-1] = []
            h_coord_dict[-1] = []
        else: # 'oh'
            o_index_dict[-1] = []
            o_coord_dict[-1] = []
            h_coord_dict[-1] = []        
            
    #print(o_index_dict)
    return o_index_dict,o_coord_dict,MDstep,h_coord_dict


if __name__ == '__main__':

    trj_path = '../'
    save_path= './'
    trj_file = 'bubble_adjust_center.lammpstrj'
    #trj_file = 'bubble_10w.lammpstrj'
    findpart = 3
    trj_start= 5000
    trj_stop = 7501
    step     = 1
    scheme   = 'h3o'
    system   = 'bulk' # slab sphere
    ion_num  = 20
    ion_num_upper = ion_num*2+2
    
    trj_name = trj_file.split('.')[0]
    if scheme == 'h3o':
        cutoff   = 1.25
        indexfile = open(os.path.join(save_path, "{0}_ohidx_dump{1}_cut{2}_{3}_{4}.txt".format(scheme,step,cutoff,trj_name,findpart)), 'w+')
        h_xyzfile = open(os.path.join(save_path, "{0}_h_xyz_dump{1}_cut{2}_{3}_{4}.txt".format(scheme,step,cutoff,trj_name,findpart)), 'w+')
        o_xyzfile = open(os.path.join(save_path, "{0}_o_xyz_dump{1}_cut{2}_{3}_{4}.txt".format(scheme,step,cutoff,trj_name,findpart)), 'w+')
        indexfile.write('# %6s%16s' %('Ions','MD_step'))
        h_xyzfile.write('# %6s%16s' %('Ions','MD_step'))
        o_xyzfile.write('# %6s%16s' %('Ions','MD_step'))
    
        for i in range(ion_num_upper):
            indexfile.write('%8s' %('O'+str(i+1)))
            h_xyzfile.write('%8s%16s%16s%16s' %('O'+str(i+1),'X','Y','Z'))
            o_xyzfile.write('%8s' %('O'+str(i+1)))
            for j in range(3):
                indexfile.write('%8s' %('H'+str(j+1)))
                h_xyzfile.write('%8s%16s%16s%16s' %('H'+str(j+1),'X','Y','Z'))
            o_xyzfile.write('%16s%16s%16s' %('X','Y','Z'))
        indexfile.write('\n')
        h_xyzfile.write('\n')         
        o_xyzfile.write('\n')
    else: # 'oh'
        cutoff   = 1.275
        indexfile = open(os.path.join(save_path, "{0}_ohidx_dump{1}_cut{2}_{3}_{4}.txt".format(scheme,step,cutoff,trj_name,findpart)), 'w+')
        h_xyzfile = open(os.path.join(save_path, "{0}_h_xyz_dump{1}_cut{2}_{3}_{4}.txt".format(scheme,step,cutoff,trj_name,findpart)), 'w+')
        o_xyzfile = open(os.path.join(save_path, "{0}_o_xyz_dump{1}_cut{2}_{3}_{4}.txt".format(scheme,step,cutoff,trj_name,findpart)), 'w+')
        indexfile.write('# %6s%16s' %('Ions','MD_step'))
        h_xyzfile.write('# %6s%16s' %('Ions','MD_step'))        
        o_xyzfile.write('# %6s%16s' %('Ions','MD_step'))
    
        for i in range(ion_num_upper):
            indexfile.write('%8s' %('O'+str(i+1)))
            h_xyzfile.write('%8s%16s%16s%16s' %('O'+str(i+1),'X','Y','Z'))
            o_xyzfile.write('%8s' %('O'+str(i+1)))          
            indexfile.write('%8s' %('H'))
            h_xyzfile.write('%8s%16s%16s%16s' %('H','X','Y','Z'))
            o_xyzfile.write('%16s%16s%16s' %('X','Y','Z'))
        indexfile.write('\n') 
        h_xyzfile.write('\n')
        o_xyzfile.write('\n')
    
    trj_file = os.path.join(trj_path, trj_file)
    atom_nums, xyz_index, box_ab = trj_info(trj_file)   
    
    for i in range(trj_start,trj_stop,step):
        o_index_dict,o_coord_dict,perstep,h_coord_dict = write_xyz(atom_nums, xyz_index[i], box_ab[i], trj_file, cutoff, scheme,system)
        #MD_step.append(perstep)
        
        if -1 in o_index_dict.keys():
            find_ion_num = 0
            indexfile.write('%8d%16d' %(find_ion_num,perstep))
            h_xyzfile.write('%8d%16d' %(find_ion_num,perstep))
            o_xyzfile.write('%8d%16d' %(find_ion_num,perstep))
        else:
            find_ion_num = len(o_index_dict.keys())
            indexfile.write('%8d%16d' %(find_ion_num,perstep))
            h_xyzfile.write('%8d%16d' %(find_ion_num,perstep))
            o_xyzfile.write('%8d%16d' %(find_ion_num,perstep))               
            for ii in o_index_dict.keys():
                indexfile.write('%8d' %ii)
                for jj in o_index_dict[ii]:
                    indexfile.write('%8d' %jj)
    
            for ii in h_coord_dict.keys():
                h_xyzfile.write('%8d' %ii)  # O index
                for jj in o_coord_dict[ii]: # O coord
                    h_xyzfile.write('%16.8f' %jj)
                ll = 0
                for kk in o_index_dict[ii]: # H index and coord                         
                    h_xyzfile.write('%8d%16.8f%16.8f%16.8f' %(kk,h_coord_dict[ii][ll][0],h_coord_dict[ii][ll][1],h_coord_dict[ii][ll][2]))
                    ll += 1
    
            for ii in o_coord_dict.keys():
                o_xyzfile.write('%8d' %ii)
                for jj in o_coord_dict[ii]:
                    o_xyzfile.write('%16.8f' %jj)
        
        for ii in range(ion_num_upper-find_ion_num):
            if scheme == 'h3o':
                indexfile.write('%8d%8d%8d%8d' %(-1,-1,-1,-1))
                h_xyzfile.write('%8d%16d%16d%16d%8d%16d%16d%16d%8d%16d%16d%16d%8d%16d%16d%16d' %(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1))
            else:
                indexfile.write('%8d%8d' %(-1,-1))
                h_xyzfile.write('%8d%16d%16d%16d%8d%16d%16d%16d' %(-1,-1,-1,-1,-1,-1,-1,-1))
            o_xyzfile.write('%8d%16d%16d%16d' %(-1,-1,-1,-1))
                                    
        indexfile.write('\n')
        h_xyzfile.write('\n')   
        o_xyzfile.write('\n')        
    indexfile.close()
    h_xyzfile.close()
    o_xyzfile.close() 