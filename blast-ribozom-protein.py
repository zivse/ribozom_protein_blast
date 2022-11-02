import ssl  # monkey patch for BioPython 1.68 & 1.69
ssl._create_default_https_context = ssl._create_unverified_context
import os
# TODO: change name of the variable below to be more specific
PATH = os.getcwd()
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
Entrez.email = 'zivse@post.bgu.ac.il' # Enter your email address here
Entrez.api_key = '016d35b4600f9c5d1d5ced586898c3ff3a09' # Enter your API key here
import pandas as pd

"""
TODO: move file names to constants in a new file named `consts.py`
TODO: add function documentation to all functions
"""


def handle_entrez(protein_id):
    """Fetch the handle for the protein, based on protein ID as a fasta file"""
    handle = Entrez.efetch(db="protein", id=protein_id, rettype='fasta', retmode='text')
    # Save handle as a fasta file
    output = open(os.path.join(PATH, 'seq-files', f'{protein_id}.fasta'), 'w')
    output.write(handle.read())
    output.close()
    # Close the handle
    handle.close()


def blast_results_to_xml(record, protein_id):
    result_handle = NCBIWWW.qblast(program="blastp", database="refseq_protein", sequence=record.seq, hitlist_size=500)
    # Save the results of the search as an XML file MAKE SURE YOU SAVE SINCE RUNTIME IS LONG
    with open(os.path.join(PATH, 'xml-files', protein_id + '.xml'), 'w') as f:
        f.write(result_handle.read())
    result_handle.close()


def get_blast_from_xml(protein_id):
    """Parse the BLASTp results which were saved into xml file as a BLAST record object"""
    with open(os.path.join(PATH, 'xml-files', protein_id + '.xml'), 'r') as result_handle:
        blast_records = NCBIXML.read(result_handle)
    return blast_records


def results_to_csv(blast_records, protein_id):
    """iterate over all hits, print information, generate a report dataframe and a fasta of the protein sequence of all hits"""
    report = {i: [] for i in ['query_id', 'protein_name', 'organism', 'e-value',
                              'score']}  # initialize a dictionary to store the report data
    sqeuence_fasta = open(os.path.join(PATH, 'hits-files', protein_id + '_hits' + '.fasta'), 'w')  # initialize a fasta file to store the protein sequences of all hits
    for alignment in blast_records.alignments:  # iterate over all hits
        for hsp in alignment.hsps:  # iterate over all hsps
            if hsp.expect < 0.05:  # if the e-value is less than 0.05
                report['query_id'].append(alignment.title.split('|')[1])  # add the query id to the report dictionary
                report['protein_name'].append(
                    alignment.title.split('|')[2].split(',')[0])  # add the protein name to the report dictionary
                report['organism'].append(
                    alignment.title.split('[')[1].split(']')[0])  # add the organism to the report dictionary
                report['e-value'].append(hsp.expect)  # add the e-value to the report dictionary
                report['score'].append(hsp.score)  # add the alignment score to the report dictionary
                sqeuence_fasta.write(
                    '>' + alignment.title + '\n' + hsp.sbjct + '\n')  # add the protein sequence to the fasta file
    report_df = pd.DataFrame(report)  # convert the report dictionary to a dataframe
    report_df.to_csv(os.path.join(PATH, 'csv-files', protein_id + '.csv'), index=False)  # save the report dataframe to a csv file
    sqeuence_fasta.close()  # close the fasta file handle


def handle_protein(protein_id):
    """run BLASTp search for the protein and return the hits"""
    handle_entrez(protein_id)
    # Read the record after it has been saved as a file.
    record = SeqIO.read(os.path.join('seq-files', f'{protein_id}.fasta'), 'fasta')
    # Run the BLASTp search against a database of all RefSeq proteins
    # Currently set to run against the RefSeq protein database and return atleast 500 hits
    blast_results_to_xml(record, protein_id)
    blast_records = get_blast_from_xml(protein_id)
    results_to_csv(blast_records, protein_id)


with open('protein-ribozom.txt') as f:
    proteins = f.readlines()
    for protein_id in proteins:
        handle_protein(protein_id.replace('\n', ''))

