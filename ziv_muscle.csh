#!/bin/env conda run -n noam_env_min
#$ -N ziv_ribosome
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o "/gpfs0/biores/users/mishmarlab/Ziv/stdout/"
#$ -e "/gpfs0/biores/users/mishmarlab/Ziv/stderr/"

python "/gpfs0/biores/users/mishmarlab/Ziv/ribozom_protein_blast/mucsle.py"