import os
import argparse
import sys
import pandas as pd
import numpy as np
import glob
import tempfile
import subprocess
from Bio import SeqIO
from Bio import Seq
import pyBigWig
import pandas as pd
import numpy as np
import deepdish as dd
import os
os.environ['CUDA_VISIBLE_DEVICES'] ="0"
import pyfaidx
import random
import pickle as pkl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tensorflow as tf
import argparse
import json
import chrombpnet.training.utils.losses as losses
from chrombpnet.training.utils.data_utils import get_seq as get_seq
import chrombpnet.training.utils.one_hot as one_hot
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import concurrent.futures
from evolvability import load_model_wrapper
from iterative_evo_mutations import iterative_evo_mutations, update_evo_locations

AFFINITY_SCRIPT = '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/util_scripts/TFBS_affinity_hits/get_TFBS_hits_v2_20221019.R'

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
PWM_SCHEMA = ["MOTIF_NAME", "MOTIF_PWM_FWD"]
model = None
inputlen = None
regions_seqs = None

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def get_footprint_for_motif(seqs, motif, model, inputlen, batch_size):
    '''
    Returns footprints for a given motif. Motif is inserted in both the actual sequence and reverse complemented version.
    seqs input is already assumed to be one-hot encoded. motif is in sequence format.
    '''
    midpoint=inputlen//2

    w_mot_seqs = seqs.copy()
    w_mot_seqs[:, midpoint-len(motif)//2:midpoint-len(motif)//2+len(motif)] = one_hot.dna_to_one_hot([motif])

    # midpoint of motif is the midpoint of sequence
    pred_output=model.predict(w_mot_seqs, batch_size=batch_size, verbose=True)
    footprint_for_motif_fwd = softmax(pred_output[0])*(np.exp(pred_output[1])-1)

    # reverse complement the sequence
    w_mot_seqs_revc = w_mot_seqs[:, ::-1, ::-1]
    pred_output_rev=model.predict(w_mot_seqs_revc, batch_size=batch_size, verbose=True)
    footprint_for_motif_rev = softmax(pred_output_rev[0])*(np.exp(pred_output_rev[1])-1)

    # add fwd sequence predictions and reverse sesquence predictions (not we flip the rev predictions)
    counts_for_motif = np.exp(pred_output_rev[1]) - 1 + np.exp(pred_output[1]) - 1
    footprint_for_motif_tot = footprint_for_motif_fwd+footprint_for_motif_rev[:,::-1]
    footprint_for_motif =  footprint_for_motif_tot / footprint_for_motif_tot.sum(axis=1)[:, np.newaxis]

    return footprint_for_motif.mean(0), counts_for_motif.mean(0)

def process_motif(row):
    motif = row["MOTIF_NAME"]
    motif_to_insert_fwd = row["MOTIF_PWM_FWD"]
#     print("inserting motif: ", motif)
#     print(motif_to_insert_fwd)
    motif_footprint, motif_counts = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, model, INPUT_SEQLEN, 64)
    return motif, [motif_footprint, motif_counts]



def run_dms(wt_seq):
    ### run wt seqs
    wt_footprints = {}
    wt_inputs = pd.DataFrame({'MOTIF_NAME':['REF'],'MOTIF_PWM_FWD':[wt_seq]})
    future_to_row = {}
    for index, row in wt_inputs.iterrows():
        future = executor.submit(process_motif, row)
        future_to_row[future] = row

    # Collect the results as they become available
    for future in concurrent.futures.as_completed(future_to_row):
        row = future_to_row[future]
        motif, result = future.result()
        wt_footprints[motif] = result
        
    wt_inputs["motif_footprint"] = [wt_footprints[x][-1][0] for x in wt_inputs['MOTIF_NAME'].tolist()]
    wt_inputs['norm_footprint'] =  wt_inputs["motif_footprint"]/footprints_at_motifs['control'][-1][0]
    
    
    mut_seqs = generate_in_silico_mutatagenesis(wt_seq, 0, 320)
    mut_seqs.columns = ['MOTIF_NAME','MOTIF_PWM_FWD']
    
    dms_footprints = {}
    # Submit jobs for each row in the dataframe
    future_to_row = {}
    for index, row in mut_seqs.iterrows():
        future = executor.submit(process_motif, row)
        future_to_row[future] = row

    # Collect the results as they become available
    for future in concurrent.futures.as_completed(future_to_row):
        row = future_to_row[future]
        motif, result = future.result()
        dms_footprints[motif] = result
    
    mut_seqs["motif_footprint"] = [dms_footprints[x][-1][0] for x in mut_seqs['MOTIF_NAME'].tolist()]
    mut_seqs['norm_footprint'] =  mut_seqs["motif_footprint"]/footprints_at_motifs['control'][-1][0]
    mut_seqs['log2_FC'] = np.log2(mut_seqs['norm_footprint']/wt_inputs['norm_footprint'].tolist()[0])
    dms_seqs = pd.concat([wt_inputs, mut_seqs])
    
    return dms_seqs

def run_snp_marginal_footprint(ref_seq, target_seq, snp_motifs):
    ### run wt seqs
    wt_footprints = {}
    wt_inputs = pd.DataFrame({'MOTIF_NAME':['REF','TARGET'],'MOTIF_PWM_FWD':[ref_seq, target_seq]})
    future_to_row = {}
    for index, row in wt_inputs.iterrows():
        future = executor.submit(process_motif, row)
        future_to_row[future] = row

    # Collect the results as they become available
    for future in concurrent.futures.as_completed(future_to_row):
        row = future_to_row[future]
        motif, result = future.result()
        wt_footprints[motif] = result
        
    wt_inputs["motif_footprint"] = [wt_footprints[x][-1][0] for x in wt_inputs['MOTIF_NAME'].tolist()]
    wt_inputs['norm_footprint'] =  wt_inputs["motif_footprint"]/footprints_at_motifs['control'][-1][0]
    target_input = wt_inputs[wt_inputs['MOTIF_NAME'] == 'TARGET']
    
    
    snp_footprints = {}
    # Submit jobs for each row in the dataframe
    future_to_row = {}
    for index, row in snp_motifs.iterrows():
        future = executor.submit(process_motif, row)
        future_to_row[future] = row

    # Collect the results as they become available
    for future in concurrent.futures.as_completed(future_to_row):
        row = future_to_row[future]
        motif, result = future.result()
        snp_footprints[motif] = result
    
    snp_motifs["motif_footprint"] = [snp_footprints[x][-1][0] for x in snp_motifs['MOTIF_NAME'].tolist()]
    snp_motifs['norm_footprint'] =  snp_motifs["motif_footprint"]/footprints_at_motifs['control'][-1][0]
    snp_motifs = pd.concat([wt_inputs, snp_motifs])
    
    return snp_motifs

def run_wt_table(wt_inputs):
    ### run wt seqs
    wt_footprints = {}
    future_to_row = {}
    for index, row in wt_inputs.iterrows():
        future = executor.submit(process_motif, row)
        future_to_row[future] = row

    # Collect the results as they become available
    for future in concurrent.futures.as_completed(future_to_row):
        row = future_to_row[future]
        motif, result = future.result()
        wt_footprints[motif] = result
        
    wt_inputs["motif_footprint"] = [wt_footprints[x][-1][0] for x in wt_inputs['MOTIF_NAME'].tolist()]
    wt_inputs['norm_footprint'] =  wt_inputs["motif_footprint"]/footprints_at_motifs['control'][-1][0]
    
    return wt_inputs

def run_snp_table(snps_df):
    # run snp seqs
    snp_footprints = {}
    # Submit jobs for each row in the dataframe
    future_to_row = {}
    for index, row in snps_df.iterrows():
        future = executor.submit(process_motif, row)
        future_to_row[future] = row

    # Collect the results as they become available
    for future in concurrent.futures.as_completed(future_to_row):
        row = future_to_row[future]
        motif, result = future.result()
        snp_footprints[motif] = result
    
    snps_df["motif_footprint"] = [snp_footprints[x][-1][0] for x in snps_df['MOTIF_NAME'].tolist()]
    snps_df['norm_footprint'] =  snps_df["motif_footprint"]/footprints_at_motifs['control'][-1][0]
    
    return snps_df

ancestral_annotations = pd.read_csv('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/ancestral_reconstructions/mouseLineage_ancestral_path_list.txt', sep = '\t')

twist_meta = pd.read_csv('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/cactus_CREs/Twist_ortholog_CREs.csv')
twist_meta['species'] = [x.split('|')[-1] for x in twist_meta['name'].tolist()]
twist_meta = twist_meta[twist_meta['species'].isin(ancestral_annotations['AncID'].tolist())]
twist_max = twist_meta[twist_meta['target_type'] == 'max_tiles']

phyloP_scores = pd.read_csv('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/syn_CREs/evolvability_optimization/DMS_300bp_tiles_1bp_size_phyloP.bed', sep = '\t', header = None)
phyloP_scores.columns = ['chr','start','end','name','phyloP']
phyloP_fasta = dict()
for record in SeqIO.parse('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/syn_CREs/evolvability_optimization/DMS_300bp_tiles_1bp_size.fa','fasta'):
    header = record.id.split("::")[0]
    phyloP_fasta[header] = str(record.seq).upper()

phyloP_scores['seq'] = phyloP_scores['name'].map(phyloP_fasta)
phyloP_scores['index'] = [int(x.split('_')[-1]) - 1 for x in phyloP_scores['name']]
phyloP_scores['CRE'] = ['_'.join(x.split('_')[:-1]) for x in phyloP_scores['name']]


genome_fasta = pyfaidx.Fasta('/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.fa')
#Dependent variables
BASE_DIR='/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet'
MODEL_DIR=BASE_DIR + '/chrombpnet_models/endo/models/'


##################################################################
# Load model
##################################################################
tf.keras.backend.clear_session()
model_path = MODEL_DIR + 'chrombpnet_nobias.h5'
model = load_model_wrapper(model_path = model_path)
INPUT_SEQLEN=2114
OUTPUT_SEQLEN=1000

splits_dict = json.load(open('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/splits/fold_0.json'))
chroms_to_keep = set(splits_dict["test"])

regions_df = pd.read_csv('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/cactus_CREs/bg_beds/endo_negatives_10ksubset.bed', sep='\t', names=NARROWPEAK_SCHEMA)
regions_subsample = regions_df[(regions_df["chr"].isin(chroms_to_keep))]
regions_subsample = regions_subsample.sample(n = 1000, random_state = 42)
print(regions_subsample.head())
print('getting background sequences')
regions_seqs = get_seq(regions_subsample, genome_fasta, INPUT_SEQLEN)

from evolvability import get_seq, load_model_wrapper, minimize_next_generation, maximize_next_generation, generate_in_silico_mutatagenesis

footprints_at_motifs = {}

avg_response_at_tn5 = []

motif = "control"
motif_to_insert_fwd = ""
print("inserting motif: ", motif)
print(motif_to_insert_fwd)
motif_footprint, motif_counts = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, model, INPUT_SEQLEN, 64)
footprints_at_motifs[motif]=[motif_footprint,motif_counts]

num_threads = 8
executor = concurrent.futures.ThreadPoolExecutor(max_workers=num_threads)

tiles_df = pd.read_csv("/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/syn_CREs/evolvability_optimization/JB_ortholog_CREs_corrected.txt", sep = '\t')

pairs = [['Gata4_chr14_5729','Mus_musculus','Rattus_norvegicus'],
        ['Gata4_chr14_5729','Mus_musculus','Homo_sapiens'],
        ['Gata4_chr14_5729','Mus_musculus','fullTreeAnc55'],
        ['Gata4_chr14_5729','Mus_musculus','fullTreeAnc41'],
        ['Gata4_chr14_5729','Mus_musculus','fullTreeAnc239'],
         ['Gata4_chr14_5729','Rattus_norvegicus','fullTreeAnc239'],
         ['Gata4_chr14_5729','Homo_sapiens','fullTreeAnc239'],
        ['Epas1_chr17_10063','Mus_musculus','fullTreeAnc239'],
         ['Epas1_chr17_10063','Rattus_norvegicus','fullTreeAnc239'],
         ['Epas1_chr17_10063','Homo_sapiens','fullTreeAnc239'],
        ['Sparc_chr11_7211','Mus_musculus','fullTreeAnc239'],
         ['Sparc_chr11_7211','Rattus_norvegicus','fullTreeAnc239'],
         ['Sparc_chr11_7211','Homo_sapiens','fullTreeAnc239'],
        ['Bend5_chr4_8201','Mus_musculus','fullTreeAnc239'],
         ['Bend5_chr4_8201','Rattus_norvegicus','fullTreeAnc239'],
        ['Lama1_chr17_7784','Mus_musculus','fullTreeAnc239'],
        ['Lama1_chr17_7784','Rattus_norvegicus','fullTreeAnc239']]


snp_marginal_fp = pd.DataFrame()

from pairwise_extract_func import formatSeqFromMutDataframe, alignment2mutDataframe, align_seqs, extract_seqs, linkAlignRef2PhyloP, count_nucleotides_with_underscore, addphyloP2MutDf
from Levenshtein import distance as levenshtein_distance

n_iter = 25

for pair in pairs:
    print("working on ")
    print(pair)
    ref_seq, target_seq = extract_seqs(pair, tiles_df)
    phyloP_CRE = phyloP_scores[phyloP_scores['CRE'] == pair[0]]
    standard_df = pd.DataFrame({'start':[0,0],
                               'end':[len(ref_seq),len(target_seq)],
                                'Event':['ref','target'],
                               'ref_nuc':[np.nan, np.nan],
                               'target_nuc':[np.nan, np.nan]})
    standard_df['seq'] = [ref_seq, target_seq]
    standard_df['id'] = standard_df['Event'] + '_' + standard_df.index.astype(str)
    standard_df['MOTIF_NAME'] = standard_df['id']
    standard_df['MOTIF_PWM_FWD'] = standard_df['seq']
    standard_df = run_wt_table(standard_df)
    standard_df = standard_df.drop(['id','MOTIF_NAME','MOTIF_PWM_FWD'], axis = 1)
    
    aligned_seq_a, aligned_seq_b = align_seqs(ref_seq, target_seq)
    mut_df = alignment2mutDataframe(aligned_seq_a, aligned_seq_b)
    mut_df = formatSeqFromMutDataframe(aligned_seq_b, mut_df)
    if pair[1] == 'Mus_musculus':
        aligned_seq_a_idx = count_nucleotides_with_underscore(aligned_seq_a)
        aligned_phyloP_a = linkAlignRef2PhyloP(aligned_seq_a, phyloP_CRE)
        mut_df = addphyloP2MutDf(mut_df, aligned_seq_a_idx, aligned_phyloP_a)
    
    total_random_selects = pd.DataFrame()
    
    ### add snps based on random order
    orig_seqs = standard_df
    orig_seqs['step'] = 0
    for i in range(0, 50): ## sample ten times
        print("random iteration step :{}".format(i))
        random_selects = mut_df.sample(n = 1, random_state = i)
        current_seq = random_selects['seq'].tolist()[0]
        edit_distance = levenshtein_distance(ref_seq, current_seq)
        random_selects['edit_distance_from_ref'] = edit_distance
        random_selects['step'] = 1
        
        size = random_selects.shape[0]
        step = 2
        while size != 0:
            print("\tstep :{}".format(step))
            aligned_seq_a, aligned_seq_b = align_seqs(ref_seq, current_seq)
            iter_df = alignment2mutDataframe(aligned_seq_a, aligned_seq_b)
            iter_df = formatSeqFromMutDataframe(aligned_seq_b, iter_df)
            if pair[1] == 'Mus_musculus':
                aligned_seq_a_idx = count_nucleotides_with_underscore(aligned_seq_a)
                aligned_phyloP_a = linkAlignRef2PhyloP(aligned_seq_a, phyloP_CRE)
                iter_df = addphyloP2MutDf(iter_df, aligned_seq_a_idx, aligned_phyloP_a)
            if iter_df.shape[0] == 0:
                print(iter_df)
                break
            sample = iter_df.sample(n = 1, random_state = i)
            current_seq = sample['seq'].tolist()[0]
            edit_distance = levenshtein_distance(ref_seq, current_seq)
            sample['step'] = step
            sample['edit_distance_from_ref'] = edit_distance
            random_selects = pd.concat([random_selects, sample])
            size = iter_df.shape
            step += 1
        
        if ref_seq == current_seq:
            print('alignment equal now')
        else:
            print('error not same sequence')
            print(ref_seq)
            print(current_seq)
        random_selects['iter'] = i+1
        total_random_selects = pd.concat([total_random_selects, random_selects])
    
    total_random_selects['id'] = ['SNP_'+ str(x) for x in range(0, total_random_selects.shape[0])]
    total_random_selects['MOTIF_NAME'] = total_random_selects['id']
    total_random_selects['MOTIF_PWM_FWD'] = total_random_selects['seq']
    total_random_selects = run_snp_table(total_random_selects)
    total_random_selects = total_random_selects.drop(['id','MOTIF_NAME','MOTIF_PWM_FWD'], axis = 1)
    
    snp_output = pd.concat([standard_df, total_random_selects])
    snp_output[['CRE','reference','target']] = pair
    snp_marginal_fp = pd.concat([snp_marginal_fp,snp_output])

snp_marginal_fp.to_csv('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/syn_CREs/evolvability_optimization/ancestral_lineage_random_walk_snps_marginal_fp_v2.txt', sep = '\t', index = False)
    

