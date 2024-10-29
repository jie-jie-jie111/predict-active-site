# !/usr/bin/env python
# -*-coding:utf-8 -*-

# File       : descriptors_generation.py
# Description：1.This file is used to calculate quantum descriptors which are applied to training
#              2.Consolidate all descriptors to generate descriptor files

import os
import RF
import numpy as np
import re
from rdkit import Chem
import csv
import Template
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
from rdkit.Chem import rdPartialCharges
from openbabel import pybel



def get_DDEC_BOs(molecule_container,bond_orders_file):
    total_atoms = molecule_container.GetNumAtoms()
    fo = open(bond_orders_file, 'r')
    data = fo.readlines()
    fo.close()
    try:
        # Find aromatic carbon atoms
        all_aromatic_carbons = [x.GetIdx() for x in molecule_container.GetAromaticAtoms()]
        for i in range(total_atoms):
            temp = (data[i + 2]).split()
            # print(i, temp)
            current_atom = molecule_container.GetAtomWithIdx(i)
            # print(current_atom.GetSymbol())
            if current_atom.GetSymbol() == temp[0]:
                current_atom.SetProp('Bond_orders_sum', temp[-1])
                if current_atom.GetIdx() in all_aromatic_carbons:
                    #
                    my_hydrogens = []
                    for n in current_atom.GetNeighbors():
                        # print(n.GetSymbol())
                        if n.GetSymbol() == 'H':
                            my_hydrogens.append(str(n.GetIdx() + 1))
                    # print(my_hydrogens)
                    if len(my_hydrogens) > 0:
                        counter = 0
                        for line in data:
                            if re.search(' Printing BOs for ATOM # (\s)* ' + str(i + 1), line):
                                # print(line)
                                break
                            counter = counter + 1
                        # print( current_atom.GetSymbol(), str(i+1), counter, my_hydrogens)
                        line = data[counter]
                        c = 0
                        while not (line.startswith(' ======================')):
                            # print(line)
                            if line.startswith(' Bonded to the (  0,   0,   0) translated image of atom number'):
                                temp = line.split()
                                if temp[12] in my_hydrogens:
                                    # print(temp[20])
                                    current_atom.SetProp("CHB", temp[20])
                            c = c + 1
                            line = data[counter + c]

    except:
        print("*CHB Error")
    print("CHB finished")
    return molecule_container


# Calculate Fukui index
def get_FMO_fukuis(molecule_container,overlap_file,log_file):
    total_atoms = molecule_container.GetNumAtoms()
    print(total_atoms)
    fo = open(overlap_file, 'r')
    data = fo.readlines()
    fo.close()
    # Creating a 2D Coverage Matrix
    overlap_matrix = np.zeros((total_atoms, total_atoms), dtype=float)
    total_integrals = int(data[0])
    for line in data[3:(total_integrals + 3)]:
        # print('-', line)
        t = line.split()
        atom1 = int(t[0]) - 1
        atom2 = int(t[1]) - 1
        overlap = float(t[-1])
        # print(overlap)
        # if 0 <= atom1 < total_atoms and 0 <= atom2 < total_atoms:
        #     overlap_matrix[atom1, atom2] = overlap

        overlap_matrix[atom1, atom2] = overlap
        #     print(atom1, atom2, overlap)

    fo = open(log_file, 'r')
    data = fo.readlines()
    fo.close()
    # Find HOMO
    for i in range(len(data)):
        line = data[i]
        m = re.search('(\d+) alpha electrons(\s)+(\d+) beta electrons', line)
        # print(m)
        if m:
            # print('>', line, m.groups())
            if m.groups()[0] == m.groups()[2]:
                HOMO_nbr = int(m.groups()[0])
                # print('HOMO :', HOMO_nbr)
            else:
                print('*alpha and beta electrons different')
        elif line.startswith(' Alpha  occ. eigenvalues --'):
            HOMO_energy = line.split()[-1]
            # print('HOMO E:', HOMO_energy)
        elif re.search('Molecular Orbital Coefficients: ', line):
            break
    # Initialising the HOMO density
    HOMO_density = np.zeros((total_atoms))
    p = re.compile('(Eigenvalues --)(.)+ ' + HOMO_energy)
    for i in range(len(data)):
        m = p.search(data[i])
        if m:
            l = (data[i]).split()
            # print(l)
            index = 7 - l.index(HOMO_energy)
            print("!"+str(index))
            line_number = i
            # print("line_num :"+str(line_number))
            # print("*"+str(HOMO_nbr))
            # print(data[i-2])
            # print("#"+str(data[i - 2]).split()[-index])
            if not ((data[i - 2]).split()[-index+1] == str(HOMO_nbr)):
                print('HOMO Error')
                return -1
            break
    # print(line_number)
    counter = line_number
    for line in data[(counter + 1):]:
        counter = counter + 1
        l = line.split()

        if len(l) == 9:
            current_atom = int(l[1]) - 1
            # print(current_atom)
            HOMO_density[current_atom] = 0
        HOMO_density[current_atom] = HOMO_density[current_atom] + float(l[-index]) * float(l[-index])
        if (data[counter + 3]).startswith('     Eigenvalues --'):
            break
    for i in range(total_atoms):
        current_atom = molecule_container.GetAtomWithIdx(i)
        fukui = HOMO_density[i]
        for j in range(total_atoms):
            fukui = fukui + HOMO_density[i] * HOMO_density[j] * overlap_matrix[i, j]
            # print(HOMO_density[i], HOMO_density[j], overlap_matrix[i,j])
        fukui = "{:6.6f}".format(fukui)
        current_atom.SetProp("Fukui", fukui)
        # print(i+1, current_atom.GetSymbol(),fukui)
    print("Calculate Fukui finish")
    return molecule_container


# Generate descriptor files
def make_feature_file(filename,molecule_container):
    # Makes a csv file for this molecule where each atom is an entry
    # if file already exists, adds to what's there
    feature_file = 'csv file'
    my_features = 'Atom#, Atom, ' + (', '.join(script.ALL_PROPERTIES)) + '\n'

    molecule_container = Chem.RemoveHs(molecule_container)
    total_atoms = molecule_container.GetNumAtoms()

    my_descriptors = []
    for atom in molecule_container.GetAtoms():
        d = {'IDX': atom.GetIdx(), 'Symbol': atom.GetSymbol()}

        for prop in script.ALL_PROPERTIES:
            d[prop] = 'NA'

        my_descriptors.append(d)
    print("Finding equivalent atoms")
    Chem.AssignStereochemistry(molecule_container, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    found_equivalent = []

    for atom in molecule_container.GetAtoms():
        # print (atom.GetIdx(), atom.GetProp('_CIPRank'))
        if atom.GetSymbol() == 'C':
            for a in molecule_container.GetAtoms():
                if atom.GetIdx() != a.GetIdx():
                    if atom.GetProp('_CIPRank') == a.GetProp('_CIPRank'):
                        # Find equivalent atoms and save index values
                        found_equivalent.append([atom.GetIdx(), a.GetIdx()])
    for a in found_equivalent:
        a.sort()
    found_equivalent = set(map(tuple, found_equivalent))

    for i, j in found_equivalent:
        atom1 = molecule_container.GetAtomWithIdx(i)
        atom2 = molecule_container.GetAtomWithIdx(j)
        for prop in script.ALL_PROPERTIES:
            if atom1.HasProp(prop) and atom2.HasProp(prop):
                temp1 = float(atom1.GetProp(prop))
                temp2 = float(atom2.GetProp(prop))
                temp = (temp1 + temp2) / 2.0
                temp = "{:6.6f}".format(temp)
                atom1.SetProp(prop, temp)
                atom2.SetProp(prop, temp)

    for i in range(total_atoms):
        current_atom = molecule_container.GetAtomWithIdx(i)
        my_features = my_features + str(current_atom.GetIdx()) + ', ' + current_atom.GetSymbol()

        for prop in script.ALL_PROPERTIES:
            if current_atom.HasProp(prop):
                temp = current_atom.GetProp(prop)
            else:
                temp = 'NA'
            my_features = my_features + ', ' + temp

        my_features = my_features + ' \n'

    fw = open(feature_file, 'w')
    fw.write(my_features)
    fw.close()
    # print(my_features)
    print("Write in feature file.")
    return feature_file


def compile_data_dictionary(my_features_csv, my_molecule, filefolder):
    for root, dirs, files in os.walk(filefolder):
        # print("*",root)
        # print(files)
        for file in files:
            if file.split(".")[1] == "sdf":
                filename = file.split(".")[0]
                descriptors_data = {'idx': []}

                for d in script.ALL_PROPERTIES:
                    descriptors_data[d] = []

                aromatic = [a.GetIdx() for a in my_molecule.GetAromaticAtoms()]
                fw = open('csv files', 'a')
                # fw = open( '../all_testing_data.csv' , 'a')

                with open(my_features_csv) as csvfile:
                    readCSV = csv.reader(csvfile, delimiter=',')
                    l = 0
                    next(readCSV)  # This is needed to skip the header line
                    for row in readCSV:
                        atom_index = int(row[0])
                        atom = my_molecule.GetAtomWithIdx(atom_index)
                        if (atom_index in aromatic) and (atom.GetSymbol() == 'C'):
                            my_neihgbors = [x.GetSymbol() for x in atom.GetNeighbors()]
                            # print(atom_index, atom.GetSymbol(), my_neihgbors)
                            if ('H' in my_neihgbors):
                                fw.write(','.join(
                                    [filename, str(atom_index), row[2], row[3], row[4], row[5], row[6], row[7]]))
                                fw.write('\n')

                                descriptors_data['idx'].append(float(atom_index))
                                descriptors_data['Q'].append(float(row[2]))
                                descriptors_data['CHB'].append(float(row[3]))
                                descriptors_data['BO'].append(float(row[4]))
                                descriptors_data['ASA'].append(float(row[5]))
                                descriptors_data['Fukui'].append(float(row[6]))
                                descriptors_data['AN'].append(float(row[7]))
                                l = l + 1

                X = np.zeros((len(script.ALL_PROPERTIES) + 1, l))


                for prop in script.ALL_PROPERTIES:
                    # print (prop,  MAP_paths.ALL_PROPERTIES.index(prop))
                    X[script.ALL_PROPERTIES.index(prop):] = descriptors_data[prop]

                X[(len(script.ALL_PROPERTIES)):] = descriptors_data['idx']
                X = X.transpose()

    return X


def extractDescriptors(my_smi,filefolder):
    # 1. convert sdf file into molecule object
    filename = (my_smi.split(".")[0]).split("\\")[-1]
    my_molecule = (Chem.SDMolSupplier(my_smi))[0]
    my_molecule = Chem.AddHs(my_molecule)
    for root, dirs, files in os.walk(filefolder):
        # print("*",root)
        # print(files)
        for file in files:
            if "DDEC6_even_tempered_net_atomic_charges" in file.split(".")[0]:
                # print("-",file)
                net_atomic_file = root + "\\" + file
                # print(net_atomic_file)
                # my_molecule = get_DDEC_Qs(my_molecule, net_atomic_file)
            elif "DDEC6_even_tempered_bond_orders" in file.split(".")[0]:
                bond_file = root + "\\" + file
                # my_molecule = get_DDEC_BOs(my_molecule, bond_file)
            elif "overlap_populations" in file.split(".")[0]:
                overlap_file = root + "\\" + file
            else:
                pass
        # print(dirpath)
        my_molecule = get_DDEC_Qs(my_molecule, net_atomic_file)
        my_molecule = get_DDEC_BOs(my_molecule, bond_file)
        my_molecule = get_FMO_fukuis(my_molecule, overlap_file, log_file)
    # 3. compile data together
    my_features_file = make_feature_file(filename,my_molecule)
    my_data = compile_data_dictionary(my_features_file, my_molecule,filefolder)

    return my_data




