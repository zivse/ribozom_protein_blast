import os
from pathlib import Path
import pandas as pd
import re

from Bio import Entrez, SeqIO


def check_csv():
    #check the names of all the proteins in the csv files
    directory = 'csv-files'
    files = Path(directory).glob('*')
    #go over all the files in csv-files
    for file in files:
        #using df to check the names of all the proteins
        df = pd.read_csv(file)
        protein_names = df[['protein_name']]
        original_name = protein_names.iloc[0].values[0].strip()
        numpy_proteins = protein_names.to_numpy()
        for protein in numpy_proteins:
            predicted = 'PREDICTED: '+original_name
            if protein[0].strip() != original_name and protein[0].strip() != predicted:
                delete_protein_from_csv(protein[0], file)
                delete_protein_from_hits_files(protein[0], file)
        break
    #remove duplicate from csv
    remove_duplicate()


def delete_protein_from_csv(protein_to_delete, file_name):
    delete_protein = pd.read_csv(file_name)
    delete_protein.drop(delete_protein.index[(delete_protein["protein_name"] == protein_to_delete)], axis=0, inplace=True)
    delete_protein.to_csv(file_name, sep=',', index=False)
    #delete the protein from the csv fike and then send it to be deleted from the fasta


def delete_protein_from_hits_files(protein_to_delete, file_name):
    protein_id = os.path.basename(file_name)
    file_name_array = protein_id.split('.')
    fasta_file_name = 'hits-files/' + file_name_array[0] + '_hits.fasta'
    with open(fasta_file_name, 'r') as fasta_file:
        fasta_read = fasta_file.read()
        new_fasta = re.sub(r"^(.*?%s(.|\n)*?)>ref" % protein_to_delete, '>ref', fasta_read, 0, re.M)
    with open(fasta_file_name, 'w') as fasta_file:
        fasta_file.write(new_fasta)


def protein_from_animal():
    common_organisms = animals_list()
    directory = 'csv-files'
    files = Path(directory).glob('*')
    # go over all the files in csv-files
    for file in files:
        df = pd.read_csv(file)
        protein_and_organism_names = df[['protein_name', 'organism']]
        numpy_proteins = protein_and_organism_names.to_numpy()
        #check if organism in common organisms list
        for protein in numpy_proteins:
            #if not delete it
            if protein[1] not in common_organisms:
                delete_protein_from_csv(protein[0], file)
                delete_protein_from_hits_files(protein[0], file)
            break


def animals_list():
    directory = 'csv-files'
    files = Path(directory).glob('*')
    files_list = list(files)
    first_file = files_list[0]
    df = pd.read_csv(first_file)
    df = df.drop(df.index[:1])
    #take the list from the first csv
    organisms_list = df[['organism']]
    numpy_organisms = organisms_list.to_numpy()
    not_common_organisms = []
    common_organisms = []
    #for each organisms in the lists check that appear in all the other csvs
    for organism in numpy_organisms:
        for file in files_list:
            check_df = pd.read_csv(file)
            if organism[0] not in set(check_df['organism']):
                not_common_organisms.append(organism[0])
                break
        common_organisms.append(organism[0])
    df = pd.read_csv('table_of_organisms.csv')
    organism_from_table_list = df[['organism']]
    #check that all the organisms in the lists appear in the table
    numpy_organism_from_table_list = organism_from_table_list.to_numpy()
    for organism in common_organisms:
        if organism not in numpy_organism_from_table_list:
            common_organisms.remove(organism)
            #check the break
            break
    return common_organisms


def remove_duplicate():
    directory = 'csv-files'
    files = Path(directory).glob('*')
    files_list = list(files)
    #for each csv find the duplicate and select the one with the lowest e-value
    for file in files_list:
        df = pd.read_csv(file)
        df1 = df.loc[df.groupby('organism', sort=False)['e-value'].idxmin()]
        df1.to_csv(file, sep=',', index=False)


def record_check(ID, type = 'gb'):

  """
  Check whether org genbank exists, if not download it, based on RefSeq ID.

    Parameters
    ----------
    ID : str
        ID of the organism.
    type : str
        Type of the record.
        Default is 'gb'.
        Other options are 'fasta'.
  """
  PATH = os.getcwd()
  suffix = '.gbk' if type == 'gb' else '.fasta'
  loc = 'genbank_DB' if type == 'gb' else 'fasta_DB'
  if not os.path.isdir(loc):
    os.mkdir(loc) # create directory if it does not exist
  filename = os.path.join(PATH, loc, ID + suffix)
  if not os.path.isfile(filename):
    print('Downloading ' + ID + '...')
    net_handle = Entrez.efetch(db = 'nucleotide', id = ID, rettype = type, retmode = 'text')
    out_handle = open(filename, 'w')
    out_handle.write(net_handle.read())
    net_handle.close()
    out_handle.close()

def rrnL_csv():
    PATH = os.getcwd()
    common_organism = animals_list()
    df = df = pd.read_csv("table_of_organisms.csv")
    for organism in animals_list:
        organism_id = df.loc[df['organism'] == organism]['RefSeq']
        organism_locs = df.loc[df['organism'] == organism]['Gene_locations']
        organism_gene_order = df.loc[df['organism'] == organism]['Gene_order']
        record_check(
            organism_id)  # Use the above function to create the record for the organism if it does not exist (better to save handles if you have enough memory)
        record = SeqIO.read(os.path.join(PATH, 'genbank_DB', organism_id + '.gbk'),
                            'genbank')  # Load the record file (genbank format)
        seq = record.seq  # Isolate the mtDNA sequence (entire sequence)
        name = 'rrnL'  # The gene name
        ind = organism_gene_order.index(name)  # The gene location based on the list of gene order
        loc = organism_locs[ind]  # The gene location based on the list of gene locations
        start, end, strand = loc.split(':')  # Split the start, end and strand
        print(
            f'rrnL starts at {start} and ends at {end} in the {strand} strand')  # Print the start and end location of the gene
        if strand == -1:
            cur_seq = seq[int(start):int(end)].reverse_complement()  # Reverse complement ONLY IF strand is -1
        with open('test_fasta.fasta', 'w') as fasta:  # Write the sequence to a fasta file
            fasta.write(f'>{organism_id}_{name}\n{cur_seq}\n')  # This is just for one organism, need to do for all

# if __name__ == '__main__':
#     # animals_list()
#     #check_csv()
#     protein_from_animal()
