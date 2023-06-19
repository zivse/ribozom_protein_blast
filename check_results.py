import ssl  # monkey patch for BioPython 1.68 & 1.69
ssl._create_default_https_context = ssl._create_unverified_context


import os
from pathlib import Path
import pandas as pd
import pickle
from Bio import Entrez, SeqIO, Seq
from ast import literal_eval # Used to convert the list column to a real list
Entrez.email = 'zivse@post.bgu.ac.il' # Enter your email address here
Entrez.api_key = '016d35b4600f9c5d1d5ced586898c3ff3a09' # Enter your API key here
from consts import HITS_FILES_NAME, CSV_FILES_NAME


def check_csv():
    """check that all the protein in the csv have right names and delete the ones that dont"""
    #check the names of all the proteins in the csv files
    files = generate_csv_files_list()
    #go over all the files in csv-files
    for file in files:
        #using df to check the names of all the proteins
        df = pd.read_csv(file)
        protein_names = df[['protein_name']]
        original_name = protein_names.iloc[0].values[0].strip()
        numpy_proteins = protein_names.to_numpy()
        # for protein in numpy_proteins:
        #     protein_0 = protein[0]
        #     predicted = 'PREDICTED: '+original_name
        #     if protein_0.strip() != original_name and protein_0.strip() != predicted:
        #         delete_protein_from_csv(protein_0, file)

    remove_duplicate_organisms_from_csv_files()


def delete_protein_from_csv(protein_to_delete, file_name):
    delete_protein = pd.read_csv(file_name)
    delete_protein.drop(delete_protein.index[(delete_protein["protein_name"] == protein_to_delete)], axis=0, inplace=True)
    delete_protein.to_csv(file_name, sep=',', index=False)


def delete_protein_from_csv_by_organism(organism_to_delete, file_name):
    delete_protein = pd.read_csv(file_name)
    delete_protein.drop(delete_protein.index[(delete_protein["organism"] == organism_to_delete)], axis=0, inplace=True)
    delete_protein.to_csv(file_name, sep=',', index=False)


def escape_characters(string):
    """
    escapes the '[' and ']' characters
    """
    return string.replace('[', '\[').replace(']', '\]')


def animals_list():
    """returns a list of strings representing the common organizems and saving them to a file"""
    files_list = generate_csv_files_list()
    first_file = files_list[0]
    df = pd.read_csv(first_file)
    #take the list from the first csv
    organisms_list = df[['organism']]
    numpy_organisms = organisms_list.to_numpy()
    common_organisms, not_common_organisms = generate_common_and_noncommon_organisms(numpy_organisms, files_list)
    # common_organisms = compare_common_with_table(common_organisms)
    for not_common_org in not_common_organisms:
        if not_common_org in common_organisms:
            common_organisms.remove(not_common_org)
    with open('common_organisms', 'w') as f:  # Save the list to a file
        for organism in common_organisms:
            print(organism)
            f.write(organism + '\n')
    return common_organisms


def compare_common_with_table(common_organisms):
    """check that all the organisms in the lists appear in the table"""
    df = pd.read_csv('../bio_projects_files/table_of_organisms.csv')
    organism_from_table_list = df[['organism']]
    numpy_organism_from_table_list = organism_from_table_list.to_numpy()
    for organism in common_organisms.copy():
        if organism not in numpy_organism_from_table_list:
            common_organisms.remove(organism)

    return common_organisms


def generate_common_and_noncommon_organisms(numpy_organisms, files_list):
    """for each organism check if exist in all the proteins, and add to the common/ noncommon lists"""
    not_common_organisms = []
    common_organisms = []
    # for each organisms in the lists check that appear in all the other csvs
    for organism in numpy_organisms:
        for file in files_list:
            check_df = pd.read_csv(file)
            if organism[0] not in set(check_df['organism']):
                not_common_organisms.append(organism[0])
                break
        common_organisms.append(organism[0])
    return common_organisms, not_common_organisms


def protein_from_animal():
    """delete all the information related to organisms that dont exist in the common organisms list """
    common_organisms = animals_list()
    files = generate_csv_files_list()
    # go over all the files in csv-files
    for file in files:
        df = pd.read_csv(file)
        protein_and_organism_names = df[['protein_name', 'organism']]
        numpy_proteins = protein_and_organism_names.to_numpy()
        #check if organism in common organisms list
        for protein in numpy_proteins:
            #if not delete it
            if protein[1] not in common_organisms:
                delete_protein_from_csv_by_organism(protein[1], file)


def remove_duplicate_organisms_from_csv_files():
    """"#for each csv find the duplicate and select the one with the lowest e-value"""
    files_list = generate_csv_files_list()
    for file in files_list:
        df = pd.read_csv(file)
        df1 = df.loc[df.groupby('organism', sort=False)['e-value'].idxmin()]
        df1.to_csv(file, sep=',', index=False)


def generate_csv_files_list():
    directory = CSV_FILES_NAME
    files = Path(directory).glob('*.csv')
    files_list = list(files)
    return files_list


def generate_fasta_files_list():
    directory = HITS_FILES_NAME
    files = Path(directory).glob('*.fasta')
    files_list = list(files)
    return files_list


def record_check( ID, type = 'gb'):
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
    # print('Downloading ' + ID + '...')
    net_handle = Entrez.efetch(db = 'nucleotide', id = ID, rettype = type, retmode = 'text')
    out_handle = open(filename, 'w')
    out_handle.write(net_handle.read())
    net_handle.close()
    out_handle.close()


def find_protein_seq(PATH,organism_id):
    record_check(organism_id)
    record = SeqIO.read(os.path.join(PATH, 'genbank_DB', organism_id + '.gbk'),
                        'genbank')  # Load the record file (genbank format)
    seq = record.seq  # Isolate the mtDNA sequence (entire sequence)
    return seq


def gene_csv(gene):
    """
    Receive a gene name, isolate its sequence from all organisms in the common_organisms list, and save all sequences into a fasta file.

    Parameters
    ----------
    gene : str
        Name of the gene.
    """
    PATH = os.getcwd()
    if os.path.exists('common_organisms'): # Load the list from file if it exists
        with open('common_organisms', 'rb') as f:
            common_organisms = pickle.load(f)  # If the list doesn't exist, create it

    else:
        common_organisms = animals_list()
    df = pd.read_csv("../bio_projects_files/table_of_organisms.csv")
    for c in ['Gene_locations', 'Gene_order']: # Convert the column values to lists
        df[c] = df[c].apply(literal_eval)
    if "Homo sapiens" not in common_organisms:
        common_organisms.insert(0, "Homo sapiens")
    fasta_lines = ''
    for organism_name in common_organisms:
        try: # Try to find the gene in the organism
            organism_id = df.loc[df['organism'] == organism_name]['RefSeq'].iloc[0]  # Get the RefSeq ID of the organism, .iloc[0] is required because otherwise it returns a series
            organism_locs = df.loc[df['organism'] == organism_name]['Gene_locations'].iloc[0]
            organism_gene_order = df.loc[df['organism'] == organism_name]['Gene_order'].iloc[0]
        except IndexError: # If the organism is not in the table, skip it
            print('ERROR: Could not find ' + organism_name + ' in the table')
            continue
        seq = find_protein_seq(PATH, organism_id)
        try: seq[0]
        except Seq.UndefinedSequenceError:
            print('ERROR: Could not find ' + organism_name + ' in the genbank DB bad sequence!')
            continue
        name = gene  # The gene name
        try:
            ind = organism_gene_order.index(name)  # The gene location based on the list of gene order
        except ValueError: # If the gene is not in the organism, skip it
            print('Could not find ' + name + ' in ' + organism_name)
            continue
        loc = organism_locs[ind]  # The gene location based on the list of gene locations
        print(loc)
        start, end, strand = loc.split(':')  # Split the start, end and strand
        print(f'{gene} starts at {start} and ends at {end} in the {strand} strand')  # Print the start and end location of the gene
        if strand == -1:
            cur_seq = seq[int(start):int(end)].reverse_complement()  # Reverse complement ONLY IF strand is -1
        else:
            cur_seq = seq[int(start):int(end)]
        fasta_lines += f'>{organism_name}_{name}\n{cur_seq}\n'  # This is just for one organism, need to do for all
    with open(f'{gene}_fasta.fasta', 'w') as fasta:  # Write the sequence to a fasta file. Ziv - If you want to write multiple sequences, you need to append to the file or move it out of the loop and write all at once
        fasta.write(fasta_lines)


def remove_proteins_from_hits_if_not_in_ids(protein_ids, hits_file):
    new_content = []
    with open(hits_file) as fh:
        for line1 in fh:
            line2 = next(fh)
            for protein_id in protein_ids:
                if protein_id[0] in line1:
                    new_content.append(line1)
                    new_content.append(line2)
    with open(hits_file, 'w') as fasta_file:
        content_string = ''.join(new_content)
        fasta_file.write(content_string)


def sync_hits_files_with_csv_files():
    files_list = generate_csv_files_list()
    for file in files_list:
        df = pd.read_csv(file)
        query_ids = df[['query_id']].to_numpy()
        protein_file_path = os.path.basename(file)
        file_name_array = protein_file_path.split('.')
        fasta_file_name = HITS_FILES_NAME + '/' + file_name_array[0] + '_hits.fasta'
        remove_proteins_from_hits_if_not_in_ids(query_ids, fasta_file_name)


if __name__ == '__main__':
    # gene_csv('rrnL')
    # gene_csv('rrnS')
    # animals_list()
    # check_csv()
     remove_duplicate_organisms_from_csv_files()
     protein_from_animal()
    # print(animals_list())
    # sync_hits_files_with_csv_files()
    # print(generate_files_list())
    #df = pd.read_csv(pathlib.PosixPath('csv-files/O15235.csv'))

