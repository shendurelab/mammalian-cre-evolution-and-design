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

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# Small inputs ship next to this script under data/.
# Two large inputs (mm10 fasta, chromBPNet nobias model) are NOT shipped;
# set MM10_FASTA and CHROMBPNET_NOBIAS_MODEL in the environment, or edit
# the defaults below.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
# reuse CRE tile FASTA that already ships with pCREs_v2
CRE_FASTA = os.path.join(
    SCRIPT_DIR, '..', 'pCREs_v2', 'data',
    'DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa',
)
FOLD_JSON = os.path.join(DATA_DIR, 'fold_0.json')
BG_REGIONS_BED = os.path.join(DATA_DIR, 'endo_negatives_10ksubset.bed')

MM10_FASTA = os.environ.get(
    'MM10_FASTA',
    '/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.fa',
)
CHROMBPNET_NOBIAS_MODEL = os.environ.get(
    'CHROMBPNET_NOBIAS_MODEL',
    '/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/chrombpnet_models/endo/models/chrombpnet_nobias.h5',
)
OUTPUT_TSV = os.environ.get('OUTPUT_TSV', 'compact_marginal_fp.txt')

for p, label in [(CRE_FASTA, 'CRE tile FASTA'),
                 (FOLD_JSON, 'fold_0.json'),
                 (BG_REGIONS_BED, 'endo_negatives_10ksubset.bed'),
                 (MM10_FASTA, 'mm10 FASTA (set $MM10_FASTA)'),
                 (CHROMBPNET_NOBIAS_MODEL, 'chromBPNet nobias model (set $CHROMBPNET_NOBIAS_MODEL)')]:
    if not os.path.exists(p):
        sys.exit("ERROR: {} not found: {}".format(label, p))

def load_model_wrapper(model_path):
    # read .h5 model
    custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_path, compile=False)
    print("got the model")
    model.summary()
    return model


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

genome_fasta = pyfaidx.Fasta(MM10_FASTA)


##################################################################
# Load model
##################################################################
tf.keras.backend.clear_session()
model = load_model_wrapper(model_path = CHROMBPNET_NOBIAS_MODEL)
INPUT_SEQLEN=2114
OUTPUT_SEQLEN=1000

splits_dict = json.load(open(FOLD_JSON))
chroms_to_keep = set(splits_dict["test"])

regions_df = pd.read_csv(BG_REGIONS_BED, sep='\t', names=NARROWPEAK_SCHEMA)
regions_subsample = regions_df[(regions_df["chr"].isin(chroms_to_keep))]
regions_subsample = regions_subsample.sample(n = 1000, random_state = 42)
print(regions_subsample.head())
print('getting background sequences')
regions_seqs = get_seq(regions_subsample, genome_fasta, INPUT_SEQLEN)


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

seq_dict = dict()

for record in SeqIO.parse(CRE_FASTA,'fasta'):
    header= str(record.id)
    cre_id = header.split("_tile_")[1]
    seq_dict[cre_id] = str(record.seq)
    

def one_bp_deletions(string):
    deletions = []
    del_names = []
    for i in range(len(string)):
        deletion = string[:i] + string[i+1:]
        deletions.append(deletion)
        del_names.append('del_'+str(i+1)) #### note this is 1-indexed
    return del_names, deletions

def extract_seqs(pair, df):
    """
    takes a list pair containing [CRE, ref, target] and extracts sequence accordingly from df
    
    returns ref and target seq
    """
    cre_df = df[df['CRE'] == pair[0]]
    ref_seq = cre_df[cre_df['species'] == pair[1]]['sequence'].tolist()[0]
    
    return ref_seq



pairs = [['Gata4_chr14_5729','Mus_musculus'],
        ['Epas1_chr17_10063','Mus_musculus'],
        ['Sparc_chr11_7211','Mus_musculus'],
        ['Bend5_chr4_8201','Mus_musculus'],
        ['Lama1_chr17_7784','Mus_musculus']]


snp_marginal_fp = pd.DataFrame()

for cre in seq_dict:
    print("working on {}".format(cre))
    ref_seq = seq_dict[cre]
    standard_df = pd.DataFrame({'id':['ref'],
                               'length':[len(ref_seq)],
                                'num_del':[0]})
    standard_df['seq'] = [ref_seq]
    standard_df['MOTIF_NAME'] = standard_df['id']
    standard_df['MOTIF_PWM_FWD'] = standard_df['seq']
    standard_df = run_wt_table(standard_df)
    standard_df = standard_df.drop(['MOTIF_NAME','MOTIF_PWM_FWD'], axis = 1)
    
    n_del = len(ref_seq)
    current_seq = ref_seq
    
    snp_output = standard_df
    
    for i in range(0, n_del):
        if i == n_del:
            pass
        else:
            del_names, deletions = one_bp_deletions(current_seq)
            del_df = pd.DataFrame({'id':del_names,
                                   'length':[len(x) for x in deletions],
                                  'seq':deletions})
            del_df['MOTIF_NAME'] = del_df['id']
            del_df['MOTIF_PWM_FWD'] = del_df['seq']
            del_df = run_snp_table(del_df)
            del_df = del_df.drop(['MOTIF_NAME','MOTIF_PWM_FWD'], axis = 1)
            del_df = del_df.sort_values(by = ['norm_footprint'], ascending = False)

            best_rank = del_df.head(n = 1)
            best_rank['num_del'] = i+1

            current_seq = best_rank['seq'].tolist()[0]
            snp_output = pd.concat([snp_output, best_rank])

    snp_output[['CRE']] = cre
    snp_marginal_fp = pd.concat([snp_marginal_fp,snp_output])

snp_marginal_fp.to_csv(OUTPUT_TSV, sep = '\t', index = False)
    

