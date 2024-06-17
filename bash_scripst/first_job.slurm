#1/bin/bash

#SBATCH --time=12:00:00
#SBATCH --partition=sched_mit_ejossou
#SBATCH -N 1
#SBATCH --mem=2G
#SBATCH -o output_myfirstjob_%j.txt
#SBATCH -e error_myfirstjob_%j.txt
#SBATCH --mail-type= BEGIN,END,FAIL
#SBATCH --mail-user=trem@mit.edu

module load mit/matlab
module load gcc/12.2.0-x86_64
module load python/3.10.8-x86_64


export VASP_PP_PATH=/nobackup1/users/trem/PBE/



#python file for alloy analysis
source myenv/bin/activate
#export PYTHONPATH="${PYTHONPATH}:/home/trem/.local/lib/python3.11/site-packages:"
python alloy_pyscript.py