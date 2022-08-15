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
    #delete from fasta - doesnt work yet
    protein_id = os.path.basename(file_name)
    file_name_array = protein_id.split('.')
    fasta_file_name = 'hits-files/' + file_name_array[0] + '_hits.fasta'
    with open(fasta_file_name, 'r') as fasta_file:
        fasta_read = fasta_file.read()
        new_fasta = re.sub(r"^(.*?%s(.|\n)*?)>ref" % protein_to_delete, '>ref', fasta_read, 0, re.M)
    with open(fasta_file_name, 'w') as fasta_file:
        fasta_file.write(new_fasta)


if __name__ == '__main__':
    check_csv()
