import ssl  # monkey patch for BioPython 1.68 & 1.69

from consts import XML_FOLDER_NAME, HITS_FILES_NAME, SEQ_FILES_NAME, CSV_FILES_NAME

ssl._create_default_https_context = ssl._create_unverified_context
import os
import argparse

PATH = os.getcwd()
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
Entrez.email = 'zivse@post.bgu.ac.il' # Enter your email address here
Entrez.api_key = '016d35b4600f9c5d1d5ced586898c3ff3a09' # Enter your API key here
import pandas as pd


def handle_entrez(protein_id):
    """Fetch the handle for the protein, based on protein ID as a fasta file"""
    handle = Entrez.efetch(db="protein", id=protein_id, rettype='fasta', retmode='text')
    # Save handle as a fasta file
    output = open(os.path.join(PATH, SEQ_FILES_NAME, f'{protein_id}.fasta'), 'w')
    output.write(handle.read())
    output.close()
    # Close the handle
    handle.close()


def blast_results_to_xml(record, protein_id):
    result_handle = NCBIWWW.qblast(program="blastp", database="refseq_protein", sequence=record.seq, hitlist_size=1000)
    # Save the results of the search as an XML file MAKE SURE YOU SAVE SINCE RUNTIME IS LONG
    with open(os.path.join(PATH, XML_FOLDER_NAME, protein_id + '.xml'), 'w') as f:
        my_xml = result_handle.read()
        my_xml = my_xml.replace("CREATE_VIEW", "")
        # create_view_index = my_xml.find("CREATE_VIEW")
        # if create_view_index != -1:
        #     my_xml = my_xml[:create_view_index] + my_xml[create_view_index + 11:]
        f.write(my_xml)
    result_handle.close()


def get_blast_from_xml(protein_id):
    """Parse the BLASTp results which were saved into xml file as a BLAST record object"""
    with open(os.path.join(PATH, XML_FOLDER_NAME, protein_id + '.xml'), 'r') as result_handle:
        blast_records = NCBIXML.read(result_handle)
    return blast_records


def results_to_csv(blast_records, protein_id):
    """iterate over all hits, print information, generate a report dataframe and a fasta of the protein sequence of all hits"""
    report = {i: [] for i in ['query_id', 'protein_name', 'organism', 'e-value',
                              'score']}  # initialize a dictionary to store the report data
    sqeuence_fasta = open(os.path.join(PATH, HITS_FILES_NAME, protein_id + '_hits' + '.fasta'), 'w')  # initialize a fasta file to store the protein sequences of all hits
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
    report_df.to_csv(os.path.join(PATH, CSV_FILES_NAME, protein_id + '.csv'), index=False)  # save the report dataframe to a csv file
    sqeuence_fasta.close()  # close the fasta file handle


def handle_protein(protein_id):
    """run BLASTp search for the protein and return the hits"""
    handle_entrez(protein_id)
    # Read the record after it has been saved as a file.
    record = SeqIO.read(os.path.join(SEQ_FILES_NAME, f'{protein_id}.fasta'), 'fasta')
    # Run the BLASTp search against a database of all RefSeq proteins
    # Currently set to run against the RefSeq protein database and return atleast 500 hits
    blast_results_to_xml(record, protein_id)
    blast_records = get_blast_from_xml(protein_id)
    results_to_csv(blast_records, protein_id)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_path', help='foo help')
    args = parser.parse_args()
    file_path = args.file_path
    with open(file_path) as f:
        proteins = f.readlines()
        for protein_id in proteins:
            handle_protein(protein_id.replace('\n', ''))