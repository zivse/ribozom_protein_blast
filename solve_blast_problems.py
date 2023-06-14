import ssl  # monkey patch for BioPython 1.68 & 1.69

from consts import XML_FOLDER_NAME, HITS_FILES_NAME, SEQ_FILES_NAME, CSV_FILES_NAME

ssl._create_default_https_context = ssl._create_unverified_context
import os
import argparse
from pathlib import Path
PATH = os.getcwd()
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
Entrez.email = 'zivse@post.bgu.ac.il' # Enter your email address here
Entrez.api_key = '016d35b4600f9c5d1d5ced586898c3ff3a09' # Enter your API key here
import pandas as pd


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
    all_organisms = []
    alraedy_done_organisms = []
    for alignment in blast_records.alignments:  # iterate over all hits
        organisms = get_organisms_from_title(alignment.title, all_organisms)
        all_organisms.extend(organisms)
        titles = split_title(alignment.title)
        hsp = get_max_score_hsp(alignment.hsps)
        if hsp.expect < 0.05:  # if the e-value is less than 0.05
            for organism in organisms:
                if alignment.title.split('|')[2].split(',')[0] == 'ref':
                    start = 2
                else:
                    start = 0
                report['query_id'].append(alignment.title.split('|')[1 + start])  # add the query id to the report dictionary
                report['protein_name'].append(
                alignment.title.split('|')[2+start].split('[')[0].replace(' ', '_').replace(',', '_'))  # add the protein name to the report dictionary
                report['organism'].append(organism)  # add the organism to the report dictionary
                report['e-value'].append(hsp.expect)  # add the e-value to the report dictionary
                report['score'].append(hsp.score)  # add the alignment score to the report dictionary
            for title in titles:
                if title.split('[')[1].split(']')[0] not in alraedy_done_organisms:
                    sqeuence_fasta.write(
                        '>' + title + '\n' + hsp.sbjct + '\n')  # add the protein sequence to the fasta file
                    # print(hsp.sbjct)
                    alraedy_done_organisms.append(title.split('[')[1].split(']')[0])

    report_df = pd.DataFrame(report)  # convert the report dictionary to a dataframe
    report_df.to_csv(os.path.join(PATH, CSV_FILES_NAME, protein_id + '.csv'), index=False)  # save the report dataframe to a csv file
    sqeuence_fasta.close()  # close the fasta file handle


def get_organisms_from_title(title, all_organisms):
    title_array = title.split('>ref')
    title_organisms_array = []
    for ref in title_array:
        if 'PREDICTED' in ref or 'predicted' in ref:
            break
        organism = ref.split('[')[1].split(']')[0]
        organism = organism.replace(' ', '_')
        words = organism.split()
        # if(len(words) > 3 ):
        #     organism = words[0] +  + words[1] + " " + words[2]
        if organism not in set(title_organisms_array) and organism not in set(all_organisms):
            title_organisms_array.append(organism)
        # print(organism)
    return title_organisms_array


def split_title(title):
    title = title[:-1]
    split_title = []
    title_array = title.split('>')
    for organism in title_array:
        if 'PREDICTED' in organism or 'predicted' in organism:
            break
        if organism != '':
            if ']' not in organism:
                organism = organism + ']'
            split_title.append(organism)
    return split_title


def get_max_score_hsp(hsps):
    score = hsps[0].score
    max_hsp = hsps[0]
    for hsp in hsps:
        if hsp.score > score:
            score = hsp.score
            max_hsp = hsp
    return max_hsp


if __name__ == '__main__':
    file_path = '/Users/zivseker/Desktop/Projects/bio-project/protein-ribozom.txt'
    with open(file_path) as f:
        proteins = f.readlines()
        for protein_id in proteins:
            protein_id = protein_id.strip()
            blast_records = get_blast_from_xml(protein_id)
            results_to_csv(blast_records, protein_id)
    # blast_records = get_blast_from_xml('O75394')
    # results_to_csv(blast_records, 'O75394')
