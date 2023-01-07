from Bio.Align.Applications import MuscleCommandline
from pathlib import Path
import os
from consts import XML_FOLDER_NAME, MUSCLE_EXE, MUSCLE_FILES_NAME, HITS_FILES_NAME


def generate_files_list():
    files = Path(HITS_FILES_NAME).glob('*')
    return list(files)


if __name__ == '__main__':
    files_list = generate_files_list()
    for file in files_list:
        out_file = f'{MUSCLE_FILES_NAME}/{os.path.basename(file)}'
        muscle_cline = MuscleCommandline(MUSCLE_EXE, input=file, out=out_file)
        print(muscle_cline)
        muscle_cline()



