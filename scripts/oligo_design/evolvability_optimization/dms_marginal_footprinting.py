import pyBigWig
import pandas as pd
import numpy as np
import deepdish as dd
import os
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
from Bio import SeqIO


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



genome_fasta = pyfaidx.Fasta('/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.fa')
#Dependent variables
# Default chromBPNet nobias model ships in the repo at
# data/chrombpnet_models/endo/chrombpnet_nobias.h5. Override via
# $CHROMBPNET_NOBIAS_MODEL to use a different model.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
CHROMBPNET_NOBIAS_MODEL = os.environ.get(
    'CHROMBPNET_NOBIAS_MODEL',
    os.path.join(_REPO_ROOT, 'data', 'chrombpnet_models', 'endo', 'chrombpnet_nobias.h5'),
)


##################################################################
# Load model
##################################################################
tf.keras.backend.clear_session()
model_path = CHROMBPNET_NOBIAS_MODEL
model = load_model_wrapper(model_path = model_path)
INPUT_SEQLEN=2114
OUTPUT_SEQLEN=1000


splits_dict = json.load(open('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/splits/fold_0.json'))
chroms_to_keep = set(splits_dict["test"])

regions_df = pd.read_csv('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/cactus_CREs/bg_beds/endo_negatives_10ksubset.bed', sep='\t', names=NARROWPEAK_SCHEMA)
regions_subsample = regions_df[(regions_df["chr"].isin(chroms_to_keep))]
regions_subsample = regions_subsample.sample(n = 500, random_state = 42)
print(regions_subsample.head())
print('getting background sequences')
regions_seqs = get_seq(regions_subsample, genome_fasta, INPUT_SEQLEN)

from evolvability import get_seq, load_model_wrapper, minimize_next_generation, maximize_next_generation, generate_in_silico_mutatagenesis

name = []
seq = []
cre = []
for record in SeqIO.parse('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/syn_CREs/evolvability_optimization/DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa','fasta'):
    header= str(record.id)
    cre_id = header.split("_tile_")[1]
    name.append(header)
    cre.append(cre_id)
    seq.append(str(record.seq))
    
max_mouse_tiles_input = pd.DataFrame({'MOTIF_NAME':name,'MOTIF_PWM_FWD':seq,'CRE':cre})

footprints_at_motifs = {}

avg_response_at_tn5 = []

motif = "control"
motif_to_insert_fwd = ""
print("inserting motif: ", motif)
print(motif_to_insert_fwd)
motif_footprint, motif_counts = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, model, INPUT_SEQLEN, 64)
footprints_at_motifs[motif]=[motif_footprint,motif_counts]

num_threads = 16
executor = concurrent.futures.ThreadPoolExecutor(max_workers=num_threads)

dms = pd.DataFrame()

### run wt seqs
wt_footprints = {}
max_mouse_tiles_input['MOTIF_NAME'] = max_mouse_tiles_input['MOTIF_NAME'] + '(WT)'
future_to_row = {}
for index, row in max_mouse_tiles_input.iterrows():
    future = executor.submit(process_motif, row)
    future_to_row[future] = row

# Collect the results as they become available
for future in concurrent.futures.as_completed(future_to_row):
    row = future_to_row[future]
    motif, result = future.result()
    wt_footprints[motif] = result

max_mouse_tiles_input["motif_footprint"] = [wt_footprints[x][-1][0] for x in max_mouse_tiles_input['MOTIF_NAME'].tolist()]
max_mouse_tiles_input['norm_footprint'] =  max_mouse_tiles_input["motif_footprint"]/footprints_at_motifs['control'][-1][0]
wt_cre_norm_footprint = dict(zip(max_mouse_tiles_input['CRE'], max_mouse_tiles_input['norm_footprint']))

dms = pd.concat([dms, max_mouse_tiles_input])

## do deep mutational scans
for index, row in max_mouse_tiles_input.iterrows():
    wt_seq = row['MOTIF_PWM_FWD']
    mut_seqs = generate_in_silico_mutatagenesis(wt_seq, 0, 300)
    mut_seqs.columns = ['MOTIF_NAME','MOTIF_PWM_FWD']
    cre_name = row['CRE']
    mut_seqs['CRE'] = cre_name
    
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
    mut_seqs['log2_FC'] = np.log2(mut_seqs['norm_footprint']/wt_cre_norm_footprint[cre_name])
    
    dms = pd.concat([dms, mut_seqs])
    
    
    
    

print(dms.head())
dms.to_csv('/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/syn_CREs/evolvability_optimization/dms_max_300bptile_selected_CREs_marginal_footprint_v3.txt', sep = '\t', index = False)

print('Done')
