#!/bin/env conda run -n noam_env_min
#$ -N ziv_ribosome
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o "/gpfs0/biores/users/mishmarlab/Ziv/stdout/"
#$ -e "/gpfs0/biores/users/mishmarlab/Ziv/stderr/"

python "/gpfs0/biores/users/mishmarlab/Ziv/ribozom_protein_blast/blast-ribozom-protein.py --file_path=/gpfs0/biores/users/mishmarlab/Ziv/ribozom_protein_blast/protein-ribozom.txt"
python "/gpfs0/biores/users/mishmarlab/Ziv/ribozom_protein_blast/check_results.py"
python "/gpfs0/biores/users/mishmarlab/Ziv/ribozom_protein_blast/create-multy-gene.py"