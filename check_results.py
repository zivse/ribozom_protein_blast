import os
from pathlib import Path
import pandas as pd
import re


def check_csv():
    # print('1')
    directory = 'csv-files'
    files = Path(directory).glob('*')
    for file in files:
        delete_protein_from_hits_files('28S ribosomal protein S10', file)
        break
        df = pd.read_csv(file)
        protein_names = df[['protein_name']]
        original_name = protein_names.iloc[0].values[0].strip()
        numpy_proteins = protein_names.to_numpy()
        for protein in numpy_proteins:
            predicted = 'PREDICTED: '+original_name
            if protein[0].strip() != original_name and protein[0].strip() != predicted:
                delete_protein_from_csv(protein[0], file)
                break
            break


def delete_protein_from_csv(protein_to_delete, file_name):
    delete_protein = pd.read_csv(file_name)
    delete_protein.drop(delete_protein.index[(delete_protein["protein_name"] == protein_to_delete)], axis=0, inplace=True)
    delete_protein.to_csv(file_name, sep=',', index=False)
    delete_protein_from_hits_files(protein_to_delete, file_name)


def delete_protein_from_hits_files(protein_to_delete, file_name):
    protein_id = os.path.basename(file_name)
    file_name_array =protein_id.split('.')
    fasta_file_name = 'hits-files/' + file_name_array[0] + '_hits.fasta'
    # print(fasta_file_name)
    with open(fasta_file_name, 'r') as fasta_file:
        new_fasta = re.sub("^>ref.*{%s}*.\n"%protein_to_delete, fasta_file)
        print(new_fasta)


if __name__ == '__main__':
    check_csv()
