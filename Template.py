# !/usr/bin/env python
# -*-coding:utf-8 -*-

# File       : Template.py
# Description：Provide some template files

from string import Template

# Job Path
Gaussian_Path = "  "
DDEC6_Path = "  "
DDEC6_ATOM_DENSITY = "  "

# 1.Gaussian gjf input file
OPT = Template("""%Chk=$TITLE.chk
#p  OPT B3LYP/6-31G(d)  pop=regular

 $TITLE

$CHARGE $MULT
$XYZ

--Link1--
%Chk=$TITLE.chk
%NoSave
#p B3LYP/6-31G(d) GFINPUT IOP(6/7=3) 6d density=current output=wfx Geom=allCheck Guess=Read

$TITLE.wfx
 
0 1

$TITLE.wfx


""")


# 2.submit gjf file
GAUSSIAN_SUBMIT = Template("""#!/bin/bash
#SBATCH -J=$TITLE          	 
#SBATCH -o out.%J     			
#SBATCH -e err.%J      			
#SBATCH -N 1                    
#SBATCH -n 24                   
#SBATCH -p  cpu48c

cd  $SLURM_SUBMIT_DIR   

g16 < $TITLE.gjf >  $TITLE.log    
""")

# 3.ddec6 control file
DDEC_JOB_CONTROL = Template("""<atomic densities directory complete path>
$PATH_TO_CHARGEMOL/atomic_densities/
</atomic densities directory complete path>

<input filename>
$TITLE.wfx
</input filename>

<charge type>
DDEC6
</charge type>

<compute BOs>
.true.
</compute BOs>

""")

# 4.ddec6 submit file
DDEC_JOB_SUBMIT = Template( """#!/bin/bash
#SBATCH --job-name=$TITLE 
#SBATCH --output=$TITLE.out 
#SBATCH --partition=cpu
#SBATCH --share  
#SBATCH --nodes 1  
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
export OMP_NUM_THREADS=4
module load gcc

echo "run complete on `hostname`: `date`" 1>&2

    """)


ALL_OPTIONS = { 'OPT': OPT,'GAUSSIAN_SUBMIT':GAUSSIAN_SUBMIT,'DDEC_JOB_CONTROL':DDEC_JOB_CONTROL,\
               'DDEC_JOB_SUBMIT': DDEC_JOB_SUBMIT}


RANDOM_FOREST_MODEL = "RF_model_eq.sav"

MINI_Iterations = 200
MAX_CONF=20

CONVERT_AU2KCAL_MOL = 627.5