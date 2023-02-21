import pickle
from pathlib import Path
import modules
import pandas as pd
from Bio import SeqIO
import json
import os

from consts import CSV_FILES_NAME, HITS_FILES_NAME


def generate_fasta_files_list():
    directory = HITS_FILES_NAME
    files = Path(directory).glob('*.fasta')
    files_list = list(files)
    return files_list


def generate_csv_files_list():
    directory = CSV_FILES_NAME
    files = Path(directory).glob('*.csv')
    files_list = list(files)
    return files_list


def read_seq_from_fasta(organism_id, file_name):
    with open(os.path.dirname(os.path.realpath(__file__)) + '/hits-files/' + file_name + '_hits.fasta', 'r') as file:
        sequences = [i for i in SeqIO.parse(file, 'fasta')]
        for seq in sequences:
            if(seq.id.split('|')[1] == organism_id):
                return (seq.seq)


def multy_gene(organism_name, gene_maping):
    multy_gene = ""
    gene_maping[organism_name] = {}
    csv_list = generate_csv_files_list()
    for file in csv_list:
        df = pd.read_csv(file)
        # organism_list = df[['organism', 'query_id']]
        line = df.loc[df['organism'] == organism_name]
        new_df = pd.DataFrame(line)
        if not new_df.empty:
            organism_id = new_df.iloc[:,0].values[0]
            protein_name = file.name.split('.')[0]
            cur_seq = read_seq_from_fasta(organism_id, protein_name)
            gene_maping[organism_name][protein_name] = len(multy_gene)
            multy_gene += str(cur_seq)

    with open('multy_gene.fasta', 'a') as write_gene:
        write_gene.write(organism_name + '\n')
        write_gene.write(multy_gene + '\n')


if __name__ == '__main__':
    gene_maping = {}
    with open('multy_gene.fasta' ,'w'):
        pass
    with open('common_organisms', 'rb') as f:
        common_organisms = pickle.load(f)
        for organism in common_organisms:
            multy_gene(organism, gene_maping)

    with open('gene_maping.json' , 'w') as file:
        file.write(json.dumps(gene_maping, indent=4))
