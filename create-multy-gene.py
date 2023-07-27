import pickle
from pathlib import Path
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
    # files = [Path('../bio_projects_files/Q9Y399.csv'), Path('../bio_projects_files/Q9P0M9.csv')]
    files_list = list(files)
    return files_list


def read_seq_from_fasta(organism_id, file_name):
    with open(os.path.dirname(os.path.realpath(__file__)) + '/hits-files/' + file_name + '_hits.fasta', 'r') as file:
        sequences = [i for i in SeqIO.parse(file, 'fasta')]
        for seq in sequences:
            if(seq.id.split('|')[1] == organism_id):
                return (seq.seq)


def human_multy_gene(gene_maping):
    protein_list = generate_csv_files_list()
    human_multy_gene = ""
    gene_maping["human"] = {}
    for protein in protein_list:
        protein_name = protein.name.split('.')[0]
        with open(os.path.dirname(os.path.realpath(__file__)) + '/seq-files/' + protein_name + '.fasta',
                  'r') as file:
            lines = file.readlines()
            i = 1
            gen_seq = ""
            while(lines[i] != '\n'):
                gen_seq += lines[i].strip()
                i += 1
            human_multy_gene += gen_seq
        gene_maping["human"][protein_name] = len(human_multy_gene)

    with open('/Users/zivseker/Desktop/Projects/bio-project/multy_gene.fasta', 'a') as write_gene:
        write_gene.write('>' + "human" + '\n')
        write_gene.write(human_multy_gene + '\n')


def multy_gene(organism_name, gene_maping):
    multy_gene = ""
    gene_maping[organism_name] = {}
    csv_list = generate_csv_files_list()
    for file in csv_list:
        df = pd.read_csv(file)
        print(df)
        # organism_list = df[['organism', 'query_id']]
        organism_name_temp = organism_name.replace('_',' ')
        line = df.loc[df['organism'] == organism_name_temp]
        new_df = pd.DataFrame(line)
        print(new_df)
        if not new_df.empty:
            organism_id = new_df.iloc[:, 0].values[0]
            organism_id = organism_id.replace(' ', '_')
            protein_name = file.name.split('.')[0]
            cur_seq = read_seq_from_fasta(organism_id, protein_name)
            if cur_seq == None:
                raise Exception("cur seq in empty " + organism_name + ' ' + protein_name+  ' '+ organism_id)
            gene_maping[organism_name][protein_name] = len(multy_gene)
            if 'o' in str(cur_seq):
                print(cur_seq)
                raise Exception("O in " + organism_name + ' ' + protein_name)
            multy_gene += str(cur_seq)
    # rrnL = get_organism_rna(organism_name, "rrnL_fasta.fasta" )
    # rrnS = get_organism_rna(organism_name, "rrnS_fasta.fasta" )
    # gene_maping[organism_name]["rrnL"] = len(multy_gene)
    # multy_gene += rrnL
    # gene_maping[organism_name]["rrnS"] = len(multy_gene)
    # multy_gene += rrnS


    with open('/Users/zivseker/Desktop/Projects/bio-project/multy_gene.fasta', 'a') as write_gene:
        organism_name = organism_name.replace(' ', '_')
        write_gene.write('>' + organism_name + '\n')
        # f.write(bytes('Text', 'utf-8'))
        write_gene.write(multy_gene + '\n')
        # write_gene.write(multy_gene + '\n')


def get_organism_rna(organism_name, file_name):
    with open(os.path.dirname(os.path.realpath(__file__)) + '/' + file_name, 'r') as file:
        lines = file.readlines()
        counter = 0
        for line in lines:
            if organism_name in line:
                return lines[counter + 1]
            counter += 1
        return "not found"


if __name__ == '__main__':
    # csv_list = generate_csv_files_list()
    # for file in csv_list:
    #     with open(file, 'r') as f:
    #         lines = f.read()
    #         lines = lines.replace('_', ' ')
    #         file.write_text(lines)

    gene_maping = {}
    with open('./multy_gene.fasta', 'w'):
        pass
    # human_multy_gene(gene_maping)
    with open('common_organisms', 'rb') as common_organisms:
        # common_organisms = pickle.load(f)
        organisms = common_organisms.readlines()
        for organism in organisms:
            organism = organism.decode("utf-8").strip()
            multy_gene(organism, gene_maping)

    with open('../gene_maping.json', 'w') as file:
        file.write(json.dumps(gene_maping, indent=4))

