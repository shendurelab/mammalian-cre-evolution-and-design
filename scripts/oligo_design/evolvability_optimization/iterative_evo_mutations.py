import os
import argparse
import sys
import pandas as pd
import numpy as np
import glob
import tempfile
import subprocess
    
def iterative_evo_mutations(seq, mut_df):
    '''
    given mut_df containing all possible snps/indels for seq
    generate all possible mutations from the list of snps/indels
    
    return a df containing all possible mutations of seqs by mut id
    '''
    
    motif_seqs = []
    
    id_snp = mut_df['id'].tolist()
    for index, row in mut_df.iterrows():
        tmp_seq = seq
        event = row['Event']
        current_start = int(row['PatternStart'])
        current_end = int(row['PatternEnd'])
        direction = row['direction']
        if event == 'mismatch':
            print(event)
            input_seq = tmp_seq[current_start:current_end]
            if input_seq != row['PatternSubstring']:
                print('error')
                print(current_start,current_end)
                print(input_seq)
                print(row['PatternSubstring'])
                print(tmp_seq)
                break
            tmp_seq = tmp_seq[:current_start] + row['SubjectSubstring'] + tmp_seq[current_end:]
            motif_seqs.append(tmp_seq)
        else:
            if event == 'deletion': ## handle deletion events
                print(event)
                print(current_start,current_end)
                print(row['PatternSubstring'])
                substract_char_len = len(row['PatternSubstring'])
                input_seq = tmp_seq[current_start:current_end]
                if input_seq != row['PatternSubstring']:
                    print('error')
                    print(current_start,current_end)
                    print(input_seq)
                    print(row['PatternSubstring'])
                    print(tmp_seq)
                    break
                tmp_seq = tmp_seq[:current_start] + tmp_seq[current_end:]
                motif_seqs.append(tmp_seq)
            else: ## handle insertion events
                print(event)
                print(current_start,current_end)
                print(row['SubjectSubstring'])
                add_char_len = len(row['PatternSubstring'])
                tmp_seq = tmp_seq[:current_start] + row['SubjectSubstring'] + tmp_seq[current_start:]
                motif_seqs.append(tmp_seq)
    
    snps_df = pd.DataFrame({'MOTIF_NAME':id_snp, 'MOTIF_PWM_FWD':motif_seqs})
    return snps_df


def update_evo_locations(mut_df, event_df):
    """
    update locations if deletion or insertion events happen
    """
    initial_starts = {int(k):int(k) for k in mut_df['PatternStart'].tolist()}
    initial_ends = {int(k):int(k) for k in mut_df['PatternEnd'].tolist()}
    
    event = event_df['Event'].tolist()[0]
    length_event = len(event_df['PatternSubstring'].tolist()[0])
    event_start = int(event_df['PatternStart'].tolist()[0])
    event_end = int(event_df['PatternEnd'].tolist()[0])
    
    if event == "deletion":
        ### update index pos to reflect substracing characters to string
        for k in initial_starts:
            if k >= event_start:
                initial_starts[k] = initial_starts[k] - length_event
        for k in initial_ends:
            if k >= event_end:
                initial_ends[k] = initial_ends[k] - length_event
                
    else:
        ### update index pos to reflect adding new characters to string
        for k in initial_starts:
            if k > event_start:
                initial_starts[k] = initial_starts[k] + length_event
        for k in initial_ends:
            if k > event_start:
                initial_ends[k] = initial_ends[k] + length_event
    print(event)
    print(event_start, event_end)
    mut_df['PatternStart'] = mut_df['PatternStart'].map(initial_starts)
    mut_df['PatternEnd'] = mut_df['PatternEnd'].map(initial_ends)
    if (mut_df['PatternStart'] > mut_df['PatternEnd']).any():
        print(mut_df[mut_df['PatternStart'] > mut_df['PatternEnd']])
#     mut_df['PatternStart'] = [initial_starts[x] if x in initial_starts else x for x in mut_df['PatternStart'].tolist()]
#     mut_df['PatternEnd'] = [initial_ends[x] if x in initial_ends else x for x in mut_df['PatternEnd'].tolist()]
    
    return mut_df