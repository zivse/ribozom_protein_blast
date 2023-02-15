from pathlib import Path
import modules
import pandas as pd
from Bio import SeqIO

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


def read_seq_from_fasta(organism_id):
    files_list = generate_fasta_files_list()
    for file in files_list:
        sequences = [i for i in SeqIO.parse(file, 'fasta')]
        for seq in sequences:
           # print(seq.id.split('|')[1])
            if(seq.id.split('|')[1] == organism_id):
                return (seq.seq)




def multy_gene(organism_name):
    multy_gene = ""
    csv_list = generate_csv_files_list()
    for file in csv_list:
        df = pd.read_csv(file)
        organism_list = df[['organism', 'query_id']]
        #print(organism_list.get(organism_name))
        #print(organism_list)
        line = df.loc[df['organism'] == organism_name]
        new_df = pd.DataFrame(line)
      #  print(line)
        if not new_df.empty:
            organism_id = new_df.iloc[:,0].values[0]
            cur_seq = read_seq_from_fasta(organism_id)
            multy_gene += str(cur_seq)
        write_gene = open(multy_gene.fasta, 'w')
        write_gene.write(multy_gene)
        #print(multy_gene)


if __name__ == '__main__':
    #read_seq_from_fasta()
    # files_list = generate_files_list();
    # for file in files_list:
    #     print(file.name)

    with open('common_organisms', 'w') as f:
        for organism in list:
            multy_gene(organism)


