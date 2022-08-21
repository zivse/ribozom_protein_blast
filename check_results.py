import os
from pathlib import Path
import pandas as pd
import re


def check_csv():
    #check the names of all the proteins in the csv files
    directory = 'csv-files'
    files = Path(directory).glob('*')
    #go over all the files in csv-files
    for file in files:
        # delete_protein_from_hits_files('28S ribosomal protein S10', file)
        # break
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
        for protein in protein_and_organism_names:
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
            break
    return common_organisms


if __name__ == '__main__':
    # animals_list()
    #check_csv()
    protein_from_animal()
