# !/usr/bin/env python
# -*-coding:utf-8 -*-

# File       : run_code.py
# Description：provide some example codes to run our project


import os
# import RF
# import pretty_output
import numpy as np
import re
from rdkit import Chem
# import csv
import Gaussian
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import rdFreeSASA
from rdkit.Chem import rdPartialCharges


# Function to generate Gaussian input file from SMILES
def smiles_to_gjf(smiles, output_file, option, name, charge=0, multiplicity=1):
    # Convert SMILES to molecule object
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
    AllChem.UFFOptimizeMolecule(mol)  # Optimize the geometry
    
    # Get 3D coordinates
    conf = mol.GetConformer()
    
    # Write the Gaussian input file
    with open(output_file, 'w') as f:
        # print(gjf_template)
        # g_input = (gjf_template.substitute(TITLE=name, CHARGE=charge,MULT=multiplicity,XYZ=info.replace(",","\n")))
        # print(g_input)
        # g_file = open(str(name)+'.gjf','w')
        
        xyz_info = []
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            line = f"{symbol}  {pos.x:.5f}  {pos.y:.5f}  {pos.z:.5f}"
            # print(line)
            xyz_info.append(line)
        print(xyz_info)
        info = " "
        for lines in xyz_info:
            if info != " ":
                info = info + "," + lines
            else:
                info = lines
        gjf_template = Gaussian.ALL_OPTIONS.get(option)
        g_input = (gjf_template.substitute(TITLE=name, CHARGE=charge,MULT=multiplicity,XYZ=info.replace(",","\n")))
        # f.write("\n")
        f.write(g_input)
        f.close()

# Example usage
# output_file = ".gjf"
# smiles_to_gjf(smi, output_file,"OPT"," ")

# Function to generate DDEC6 input file
def generate_ddec(name):
    # for name,value in coordinfo.items():
    # 1.job_control file 
    name1 = 'job_control.txt'
    template1 = Gaussian.ALL_OPTIONS.get('DDEC_JOB_CONTROL')
    job_control_input= (template1.substitute(TITLE=name,PATH_TO_CHARGEMOL=Gaussian.DDEC6_Path))
    job_control_file = open(str(name1),'w')
    job_control_file.write(job_control_input)
    job_control_file.close()
# 2.job_submit
    template2 = Gaussian.ALL_OPTIONS.get('DDEC_JOB_SUBMIT')
    job_submit_input = (template2.substitute(TITLE=name))
    job_submit_file = open(str(name2),'w')
    job_submit_file.write(job_submit_input)
    job_submit_file.close()

    return job_control_file,job_submit_file

# Example usage
# generate_ddec(' ')
