#### modified to fit JB's new probound script
from Bio import SeqIO
import time 
import subprocess
from Bio.Seq import Seq
from Bio import pairwise2
import numpy as np
import pandas as pd
import glob
import os
import seaborn as sns
import itertools
import tempfile
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_chunk(fasta_chunk, tfbs_model):
    df_all_affinities = pd.DataFrame()
    # ProBound jar path: set $PROBOUND_JAR, or pass --jar on the CLI.
    # See https://github.com/RubeLab/ProBound for building the jar.
    probound_jar = os.environ.get('PROBOUND_JAR',
        '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/util_scripts/ProBoundTools/ProBoundTools/target/ProBound-jar-with-dependencies.jar')
    probound_default_cmd = f'java -cp {probound_jar} proBoundTools/App '
    
    for record in fasta_chunk:
        header = str(record.id)
        print("Working on {} sequence...".format(header))
        seq = str(record.seq).upper()
        for index, row in tfbs_model.iterrows():
            tfbs = row['TF_name']
            tfbs_ID = row['TFID']
            tfbs_size = int(row['model_size'])
            print('\tAligning to {} sites.'.format(tfbs))
            print('\tTFBS size = {}.'.format(tfbs_size))
            kmer_length = tfbs_size

            # Create the dataframe
            df_kmer = pd.DataFrame(columns=["start_pos", "end_pos", "kmer_seq"])

            for i in range(0, len(seq)):
                start_pos = i
                end_pos = i + kmer_length
                kmer_seq = seq[start_pos:end_pos]
                if len(kmer_seq) < kmer_length:
                    break
                else:
                    tmp_kmer = pd.DataFrame({"start_pos": [start_pos], "end_pos": [end_pos], "kmer_seq": [kmer_seq]})
                    df_kmer = pd.concat([df_kmer, tmp_kmer], ignore_index=True)
            
            with tempfile.NamedTemporaryFile() as temp_substring:
                with tempfile.NamedTemporaryFile() as temp_probound:
                    df_kmer['kmer_seq'].to_csv(temp_substring.name, sep='\t', index=False, header=False)
                    cmd_str = probound_default_cmd + "-c 'loadMotifCentralModel(%s).buildConsensusModel().addNScoring().inputTXT(%s).bindingModeScores(/dev/stdout,profile)' > %s" % (tfbs_ID, temp_substring.name, temp_probound.name)
                    subprocess.check_call(cmd_str, shell=True)
                    probound_affinities = pd.read_csv(temp_probound.name, sep='\t', header=None)
                    probound_affinities.columns = ["substring","affinity_FOR","affinity_REV"]
                    probound_affinities_for = dict(zip(probound_affinities['substring'], probound_affinities['affinity_FOR']))
                    probound_affinities_rev = dict(zip(probound_affinities['substring'], probound_affinities['affinity_REV']))
                    
            df_kmer['affinity_FOR'] = df_kmer['kmer_seq'].map(probound_affinities_for)
            df_kmer['affinity_REV'] = df_kmer['kmer_seq'].map(probound_affinities_rev)
            df_kmer['TF_name'] = tfbs
            df_kmer['substring_l'] = tfbs_size
            df_kmer['CRE'] = header
            
            df_all_affinities = pd.concat([df_all_affinities, df_kmer], ignore_index=True)
    return df_all_affinities

def process_fasta(fasta, tfbs_model_path, n_threads, output):
    tfbs_model = pd.read_csv(tfbs_model_path, sep='\t')
    tfbs_model = tfbs_model.sort_values(by=['model_size'])
    tfbs_max_aff = dict(zip(tfbs_model['TF_name'], tfbs_model['max_aff']))
    
    # ProBound jar path: set $PROBOUND_JAR, or pass --jar on the CLI.
    # See https://github.com/RubeLab/ProBound for building the jar.
    probound_jar = os.environ.get('PROBOUND_JAR',
        '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/util_scripts/ProBoundTools/ProBoundTools/target/ProBound-jar-with-dependencies.jar')
    probound_default_cmd = f'java -cp {probound_jar} proBoundTools/App '
    
    # Read the fasta sequences
    fasta_sequences = list(SeqIO.parse(fasta, 'fasta'))
    chunk_size = len(fasta_sequences) // n_threads

    if chunk_size <= 0:
        df_all_affinities = process_chunk(fasta_sequences, tfbs_model)
        
    else:
        # Split sequences into chunks for parallel processing
        fasta_chunks = [fasta_sequences[i:i + chunk_size] for i in range(0, len(fasta_sequences), chunk_size)]

        df_all_affinities = pd.DataFrame()

        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = [executor.submit(process_chunk, chunk, tfbs_model) for chunk in fasta_chunks]

            for future in as_completed(futures):
                df_all_affinities = pd.concat([df_all_affinities, future.result()], ignore_index=True)
    
    ####normalize affinity
    df_all_affinities_norm = pd.DataFrame()
    for tfbs in df_all_affinities['TF_name'].unique():
        tmp = df_all_affinities[df_all_affinities['TF_name'] == tfbs]
        tmp['norm_affinity_FOR'] = tmp['affinity_FOR']/tfbs_max_aff[tfbs]
        tmp['norm_affinity_REV'] = tmp['affinity_REV']/tfbs_max_aff[tfbs]
        df_all_affinities_norm = pd.concat([df_all_affinities_norm, tmp], ignore_index = True)
    
    # Add 'TFBS_ori' column
    df_all_affinities_norm["TFBS_orientation"] = df_all_affinities_norm.apply(lambda row: "+" if row["norm_affinity_FOR"] >= row["norm_affinity_REV"] else "-", axis=1)

    # Add 'norm_affinity' column
    df_all_affinities_norm["norm_affinity"] = df_all_affinities_norm.apply(
        lambda row: row["norm_affinity_FOR"] if row["norm_affinity_FOR"] >= row["norm_affinity_REV"] else row["norm_affinity_REV"], 
        axis=1
    )
    
    df_all_affinities_norm.to_csv(output, sep='\t', index=False, compression='gzip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',required=True,help="fasta file")
    parser.add_argument('-m','--model', required=True, help="tfbs model file")
    parser.add_argument('-o','--output', required=True, help="output file name")
    parser.add_argument('-t','--threads', type=int, default = 1, help="number of threads")
    args = parser.parse_args()
    process_fasta(tfbs_model_path = args.model, fasta = args.fasta, n_threads = args.threads, output = args.output)
