#!/usr/bin/python3
"""Melchor Sanchez-Martinez 2019"""
"""Usage    python get_center_size.py pdb_file.pdb"""

import sys, os

def get_ligand(p):
    leepdb = open(pdb, "r")
    ligand = []
    for line in leepdb: 
        if line.startswith('ATOM') or  line.startswith('HETATM'):
            line=line.split()
            ligand.append(str(line[3]))
    leepdb.close()
    return(ligand[1])

def get_coordinates(p):
    reader=open(p, 'r')
    coordinates=[[], [], []]
    for line in reader:
        #print(line)
        if line.startswith('ATOM') or  line.startswith('HETATM'):
            line=line.split()
            coordinates[0].append(float(line[6]))#6 dependiendo del PDB es posicion 6 en lugar de 5
            coordinates[1].append(float(line[7]))#7 dependiendo del PDB es posicion 7 en lugar de 6
            coordinates[2].append(float(line[8]))#8 dependiendo del PDB es posicion 8 en lugar de 7
    reader.close()
    return coordinates

def get_box_center_and_box_size(coords):
    '''
    Returns the cavity
    '''
    xs = coords[0]
    ys = coords[1]
    zs = coords[2]
    bc = [((min(xs))+(max(xs)))/2, ((min(ys))+(max(ys)))/2, ((min(zs))+(max(zs)))/2]
    bs = [((max(xs)+6)-(min(xs)-6)), (((max(ys)+6)-(min(ys)-6))), ((max(zs)+6)-(min(zs)-6))]

    x = [bc[0]-bs[0]/2, bc[0]+bs[0]/2]
    y = [bc[1]-bs[1]/2, bc[1]+bs[1]/2]
    z = [bc[2]-bs[2]/2, bc[2]+bs[2]/2]
    
    return bc, bs

pdb=sys.argv[1]
coordinates=get_coordinates(pdb)
grid_center,grid_size=get_box_center_and_box_size(coordinates)
ligand_code = get_ligand(pdb)

print("Binding site of " + str(ligand_code))
print("Grid center: X " + str(grid_center[0]) +" Y "+ str(grid_center[1]) +" Z "+ str(grid_center[2]))
print("Grid size: X " + str(grid_size[0]) +" Y "+ str(grid_size[1]) +" Z "+ str(grid_size[2]))


with open('./grid_center_and_size', 'w')as fw:
    fw.write("Binding site of " + str(ligand_code) + '\n')
    fw.write("Grid center: X " + str(grid_center[0]) +" Y "+ str(grid_center[1]) +" Z "+ str(grid_center[2]) + '\n')
    fw.write("Grid size: X " + str(grid_size[0]) +" Y "+ str(grid_size[1]) +" Z "+ str(grid_size[2]) +'\n')
fw.close()