import json
from pathlib import Path
import pandas
import plotly.express as px

PSICOV_PATH = '/Users/zivseker/Desktop/Projects/bio-project/psicov-files'
FILTER_VALUE = 0.8


def get_all_file_paths():
    files = Path(PSICOV_PATH).glob('*')
    return list(files)


def get_file_name_from_path(path):
    return path.split('/')[-1]


def get_metadata_from_file_name(file_name):
    splitted_file_name = file_name.split('_')
    left_protein_name = splitted_file_name[0]
    right_protein_start_index = splitted_file_name[1]
    right_protein_name = splitted_file_name[2].split('.')[0]
    return left_protein_name, int(right_protein_start_index), right_protein_name


def get_psicov_data(file_path):
    return pandas.read_csv(file_path, delimiter=" ")


def get_num_of_valid_rows(filter_value, middle_index, df):
    counter = 0
    for _, row in df.iterrows():
        left_index = row[0]
        right_index = row[1]
        value = row[4]
        if left_index > middle_index or right_index < middle_index:
            continue

        if value < filter_value:
            continue

        counter += 1

    return counter


def get_protein_index(proteins, protein_name):
    sorted_proteins = sorted(proteins)
    for index, protein in enumerate(sorted_proteins):
        if protein_name == protein:
            return index
    return -1


def main_logic():
    proteins = set()
    file_paths = get_all_file_paths()

    if not file_paths:
        print("no files found")
        return

    for file_path in file_paths:
        file_name = get_file_name_from_path(str(file_path))
        left_protein_name, _, right_protein_name = get_metadata_from_file_name(file_name)
        proteins.add(left_protein_name)
        proteins.add(right_protein_name)

    counting_matrix = [[0 for _ in range(len(proteins))] for _ in range(len(proteins))]

    for file_path in file_paths:
        file_name = get_file_name_from_path(str(file_path))
        df = get_psicov_data(file_path)
        left_protein_name, right_protein_start_index, right_protein_name = get_metadata_from_file_name(file_name)
        left_protein_index = get_protein_index(proteins, left_protein_name)
        right_protein_index = get_protein_index(proteins, right_protein_name)
        num_of_valid_rows = get_num_of_valid_rows(FILTER_VALUE, right_protein_start_index, df)
        counting_matrix[left_protein_index][right_protein_index] = num_of_valid_rows

    return counting_matrix, proteins


if __name__ == '__main__':
    counting_matrix, proteins = main_logic()
    protein_to_index = {protein: get_protein_index(proteins, protein) for protein in proteins}
    with open("protein_to_index_mapping.txt", "w") as f:
        f.write(json.dumps(protein_to_index, indent=4))
    fig = px.imshow(counting_matrix, color_continuous_scale='RdBu_r', origin='lower')
    fig.show()




















