import ssl  # monkey patch for BioPython 1.68 & 1.69
ssl._create_default_https_context = ssl._create_unverified_context
from Bio.Blast import NCBIWWW

# sequence_data = open("sequence.fasta").read()
# print(sequence_data)
protein = "MRPS28"
with open('protein-ribozom2.txt') as f:
    proteins = f.readlines()
for protein in proteins:
    result_handle = NCBIWWW.qblast("blastp", "nt", protein)
    with open('results.xml', 'a') as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)
    print(f"done with {protein}")

