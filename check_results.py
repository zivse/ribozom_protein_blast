import os
from pathlib import Path
import pandas as pd


def check_csv():
    directory = 'csv-files'
    files = Path(directory).glob('*')
    for file in files:
        df = pd.read_csv(file)
        protein_names = df[['protein_name']]
        original_name = protein_names.iloc[0].values[0].strip()
        numpy_proteins = protein_names.to_numpy()

        for protein in numpy_proteins:
            predicted = 'PREDICTED: '+original_name
            if protein[0].strip() != original_name and protein[0].strip() != predicted:
                print(protein[0].strip())
                delete_protein_from_csv(protein[0], file)


def delete_protein_from_csv(protein_to_delete, file_name):
    delete_protein = pd.read_csv(file_name)
    delete_protein.drop(delete_protein.index[(delete_protein["protein_name"] == protein_to_delete)], axis=0, inplace=True)
    delete_protein.to_csv(file_name, sep=',',index=False)
    delete_protein_from_hits_files(protein_to_delete, file_name)


def delete_protein_from_hits_files(protein_to_delete, file_name):
    protein_id = os.path.basename(file_name)
    parts =protein_id.split('.')
    # print(parts)


if __name__ == '__main__':
    check_csv()