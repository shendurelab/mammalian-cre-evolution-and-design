import pandas as pd
import os
import numpy as np
import glob
import math
from sequence_align.pairwise import hirschberg, needleman_wunsch

# Function to merge rows with adjacent intervals and different events
def merge_adjacent_mismatch_intervals(df):
    """
    takes a mismatch 
    """
    merged_rows = []
    i = 0
    
    while i < len(df) - 1:
        # Check if the end of the current row is equal to the start of the next row
        if (df.iloc[i]['end'] == df.iloc[i+1]['start']) and (df.iloc[i]['Event'] == df.iloc[i+1]['Event']):
            # Merge the two rows
            merged_row = {
                'start': df.iloc[i]['start'],
                'end': df.iloc[i+1]['end'],
                'Event': df.iloc[i]['Event'],
                'ref_nuc': df.iloc[i]['ref_nuc'] + df.iloc[i+1]['ref_nuc'],
                'target_nuc': df.iloc[i]['target_nuc'] + df.iloc[i+1]['target_nuc']
            }
            merged_rows.append(merged_row)
            # Skip the next row
            i += 2
        else:
            # If the intervals are not adjacent or events are different, keep the current row unchanged
            merged_rows.append(df.iloc[i].to_dict())
            i += 1
    
    # Add the last row if it wasn't merged
    if i == len(df) - 1:
        merged_rows.append(df.iloc[i].to_dict())
    
    # Convert the list of dictionaries to a DataFrame
    merged_df = pd.DataFrame(merged_rows)
    
    return merged_df

def generate_target_mut_id(mut_df, reverse = False):
    mut_ids = []
    if reverse:
        for index, row in mut_df.iterrows():
            if row['Event'] == 'mismatch':
                mutid = row['ref_nuc'] + str(row['ref_start']) + row['target_nuc']
                mut_ids.append(mutid)
            elif row['Event'] == 'deletion':
                ref_pos = '_'.join([str(x) for x in list(row['ref_start'])])
                mutid = ref_pos + 'ins' + row['target_nuc']
                mut_ids.append(mutid)
            else:
                if len(row['target_nuc']) < 2:
                    mutid = str(row['ref_start']) + 'del' + row['ref_nuc']
                else:
                    ref_pos = str(row['ref_start']) + "_" + str(row['ref_end'])
                    mutid = ref_pos + 'del' + row['ref_nuc']
                mut_ids.append(mutid)
    else:
        for index, row in mut_df.iterrows():
            if row['Event'] == 'mismatch':
                mutid = row['target_nuc'] + str(row['target_start']) + row['ref_nuc']
                mut_ids.append(mutid)
            elif row['Event'] == 'insertion':
                ref_pos = '_'.join([str(x) for x in list(row['target_start'])])
                mutid = ref_pos + 'ins' + row['ref_nuc']
                mut_ids.append(mutid)
            else:
                if len(row['ref_nuc']) < 2:
                    mutid = str(row['target_start']) + 'del' + row['target_nuc']
                else:
                    ref_pos = str(row['target_start']) + "_" + str(row['target_end'])
                    mutid = ref_pos + 'del' + row['target_nuc']
                mut_ids.append(mutid)
    mut_df['target_mut_id'] = mut_ids
    return(mut_df)

def formatSeqFromMutDataframe(target_seq, mut_df):
    mut_seqs = []
    
    compare_seq = ''.join(target_seq)
    compare_seq = compare_seq.replace("_","")
    
    for index, row in mut_df.iterrows():
        start = int(row['start'])
        end = int(row['end'])
        mutation = row['Event']
        ref_nuc = row['ref_nuc']
        target_nuc = row['target_nuc']
        
        if mutation == 'insertion':
            if start == 0:
                mut_seq = [ref_nuc] + target_seq
            else:
                mut_seq = target_seq[:start] + [ref_nuc] + target_seq[start:]
            mut_seq = ''.join(mut_seq)
            mut_seq = mut_seq.replace("_","")
            if len(mut_seq) != len(compare_seq) + len(ref_nuc):
                print('wrong length in insertion')
                print(row)
                print(len(mut_seq))
                print(len(compare_seq) + len(ref_nuc))
                print(mut_seq)
                print(compare_seq)
                break
        elif mutation == 'deletion':
            mut_seq = target_seq[:start] + target_seq[end:]
            mut_seq = ''.join(mut_seq)
            mut_seq = mut_seq.replace("_","")
            
            if len(mut_seq) != len(compare_seq) - len(ref_nuc):
                print('wrong length in deletion')
                print(row)
                print(len(mut_seq))
                print(len(compare_seq) - len(ref_nuc))
                print(mut_seq)
                print(compare_seq)
                break
        else:
            mut_seq = target_seq[:start] + [ref_nuc] + target_seq[end:]
            mut_seq = ''.join(mut_seq)
            mut_seq = mut_seq.replace("_","")
            
            if len(mut_seq) != len(compare_seq):
                print('wrong length in mismatch')
                print(row)
                print(len(mut_seq))
                print(len(compare_seq) - len(ref_nuc))
                print(mut_seq)
                print(compare_seq)
                break
        mut_seqs.append(mut_seq)
    
    mut_df['seq'] = mut_seqs
    return mut_df
        

def addphyloP2MutDf(mut_df, aligned_seq_a_idx, aligned_phyloP_a):
    phyloP_col = []
    for index, row in mut_df.iterrows():
        start, end = row['ref_start'], row['ref_end']
#         width = end - start
        if row['Event'] != 'deletion':
            if math.isnan(end):
                start = aligned_seq_a_idx.index(start)
                phyloP_score = aligned_phyloP_a[start]
            elif math.isnan(start):
                end = aligned_seq_a_idx.index(end)
                phyloP_score = aligned_phyloP_a[end]
            else:
                start, end = aligned_seq_a_idx.index(start), aligned_seq_a_idx.index(end)
                phyloP_score = np.mean(aligned_phyloP_a[start:end])
            phyloP_col.append(phyloP_score)
        else:
            start = row['ref_start']
            start, end = aligned_seq_a_idx.index(start[0]), aligned_seq_a_idx.index(start[-1])
            ### for deletions take the average from left-right of the deletion
            phyloP_score = aligned_phyloP_a[end]
            phyloP_col.append(phyloP_score)
    mut_df['phyloP'] = phyloP_col
    return mut_df
    
def linkAlignRef2PhyloP(aligned_seq_a, phyloP_df):
    """
    
    """
    phyloP_dict = {(row['seq'], row['index']): row['phyloP'] for idx, row in phyloP_df.iterrows()}
    phyloP_df = phyloP_df.sort_values(by = ['index'])
    phyloP_seq = ''.join(phyloP_df['seq'].tolist())
    aligned_seq = ''.join([x for x in aligned_seq_a if x != '_'])
    if aligned_seq == phyloP_seq:
        phyloP_ind = []
        seq_i = 0
        for i, seq in enumerate(aligned_seq_a):
            if seq == '_':
                phyloP_ind.append(np.nan)
            else:
                phyloP_score = phyloP_dict[(seq, seq_i)]
                phyloP_ind.append(phyloP_score)
                seq_i += 1
        return phyloP_ind
    else:
        print(phyloP_seq)
        print(aligned_seq)

        
def linkAlignRef2DMS(aligned_seq_a, DMS_df):
    """
    
    """
    dms_dict = {(row['seq'], row['index']): row['log2FC'] for idx, row in DMS_df.iterrows()}
    DMS_df = DMS_df.sort_values(by = ['index'])
    DMS_seq = ''.join(DMS_df['seq'].tolist())
    aligned_seq = ''.join([x for x in aligned_seq_a if x != '_'])
    if aligned_seq == DMS_seq:
        dms_ind = []
        seq_i = 0
        for i, seq in enumerate(aligned_seq_a):
            if seq == '_':
                dms_ind.append(np.nan)
            else:
                dms_score = dms_dict[(seq, seq_i)]
                dms_ind.append(dms_score)
                seq_i += 1
        return dms_ind
    else:
        print(aligned_seq)
        print(aligned_seq)

def count_nucleotides_with_underscore(lst):
    counter = 0
    result = []
    
    for item in lst:
        if item == "_":
            result.append(np.nan)
        else:
            result.append(counter)
            counter += 1
    return result
        
def align_count_nucleotides_with_underscore(lst):
    result = []
    count = -1  # Start at -1 to make it zero-based
    
    for i, char in enumerate(lst):
        if char != "_":
            count += 1
            result.append(count)
        else:
            # Find the nearest left position
            left = count if count >= 0 else None
            
            # Find the nearest right position
            right = None
            for j in range(i + 1, len(lst)):
                if lst[j] != "_":
                    right = count + 1
                    break
            
            result.append((left, right))
    
    return result
        
def alignment2mutDataframe(aligned_seq_a, aligned_seq_b):
    """
    given two lists of alignments - function will parse and identify all snps/indels between two sequences.
    
    a = ref seq
    b = target seq
    
    insertion/deletion cases are directed at b seq
    """
    # Initialize variables to track positions and intervals
    positions = []
    start = None
    prev_underscore_a = False
    prev_underscore_b = False
    ref_nuc_position = align_count_nucleotides_with_underscore(aligned_seq_a)
    ref_align_position = [i for i in range(0, len(aligned_seq_a))]
    ref_nuc_position = dict(zip(ref_align_position, ref_nuc_position))
    target_nuc_position = align_count_nucleotides_with_underscore(aligned_seq_b)
    target_align_position = [i for i in range(0, len(aligned_seq_b))]
    target_nuc_position = dict(zip(target_align_position, target_nuc_position))
    

    # Iterate through both lists simultaneously
    for i, (item_a, item_b) in enumerate(zip(aligned_seq_a, aligned_seq_b)):
        # Check if item_a is "_"
        if item_a == "_":
            if not prev_underscore_a:
                start = i
                prev_underscore_a = True
        # Check if item_b is "_"
        elif item_b == "_":
            if not prev_underscore_b:
                start = i
                prev_underscore_b = True
        # If neither item is "_", check for mismatches
        elif item_a != item_b:
            if prev_underscore_a or prev_underscore_b:
                positions.append([start, i, "deletion" if prev_underscore_a else "insertion", 
                                  ''.join(aligned_seq_a[start:i]), ''.join(aligned_seq_b[start:i])])
                prev_underscore_a = False
                prev_underscore_b = False
            positions.append([i, i+1, "mismatch", ''.join(aligned_seq_a[i:i+1]), ''.join(aligned_seq_b[i:i+1])])
        # If items match and an interval is ongoing, end it
        elif prev_underscore_a or prev_underscore_b:
            positions.append([start, i, "deletion" if prev_underscore_a else "insertion", 
                              ''.join(aligned_seq_a[start:i]), ''.join(aligned_seq_b[start:i])])
            prev_underscore_a = False
            prev_underscore_b = False
            
    # If an interval is ongoing at the end, end it
    if prev_underscore_a:
        positions.append((start, len(aligned_seq_a)-1, "insertion"))
    elif prev_underscore_b:
        positions.append((start, len(aligned_seq_a)-1, "deletion"))
    
#     print(positions)
    if len(positions) == 0:
        print(positions)
        return pd.DataFrame()
    elif len(positions[0]) == 3:
        positions = pd.DataFrame(positions, columns=['start','end','Event'])
        positions['ref_start'] = positions['start'].map(ref_nuc_position)
        positions['ref_end'] = positions['end'].map(ref_nuc_position)
        positions['target_start'] = positions['start'].map(target_nuc_position)
        positions['target_end'] = positions['end'].map(target_nuc_position)
    else:
        positions = pd.DataFrame(positions, columns=['start','end','Event','ref_nuc','target_nuc'])
        mismatch_df = positions[positions['Event'] == 'mismatch']
        mismatch_df = merge_adjacent_mismatch_intervals(mismatch_df)
        positions = positions[positions['Event'] != 'mismatch']
        positions = pd.concat([positions, mismatch_df])
        positions['ref_start'] = positions['start'].map(ref_nuc_position)
        positions['ref_end'] = positions['end'].map(ref_nuc_position)
        positions['target_start'] = positions['start'].map(target_nuc_position)
        positions['target_end'] = positions['end'].map(target_nuc_position)
        
    positions = positions.sort_values(by='start')
    
    return positions

def calculate_identity(sequenceA, sequenceB):
    """
    Returns the percentage of identical characters between two sequences.
    Assumes the sequences are aligned.
    """

    sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
    gapless_sa = [sa[i] for i in range(len(sa)) if (sa[i] != '_' and sb[i] != '_')]
    gapless_sb = [sb[i] for i in range(len(sb)) if (sa[i] != '_' and sb[i] != '_')] 
    gapless_sl = sum([1 for i in range(sl) if (sa[i] != '_' and sb[i] != '_')])
    if gapless_sl == 0:
        return(0)
    else:
        matches = [gapless_sa[i] == gapless_sb[i] for i in range(gapless_sl)]
        blast_id = (100 * sum(matches)) / gapless_sl

        return blast_id

def align_seqs(ref_seq, target_seq, align_type = 'global'):
    """
    given two sequences - function will return the two alignments in list format.
    
    default global alignment
    """
    if align_type == 'global':
        aligned_seq_a, aligned_seq_b = needleman_wunsch(
            list(ref_seq),
            list(target_seq),
            match_score=1.0,
            mismatch_score=-1.0,
            indel_score=-1.0,
            gap="_",
        )
    else:
        aligned_seq_a, aligned_seq_b = hirschberg(
            list(ref_seq),
            list(target_seq),
            match_score=2.0,
            mismatch_score=-1.0,
            indel_score=-2.0,
            gap="_",
        )

    return aligned_seq_a, aligned_seq_b
    
    
def extract_seqs(pair, df):
    """
    takes a list pair containing [CRE, ref, target] and extracts sequence accordingly from df
    
    returns ref and target seq
    """
    cre_df = df[df['CRE'] == pair[0]]
    ref_seq = cre_df[cre_df['species'] == pair[1]]['sequence'].tolist()[0]
    target_seq = cre_df[cre_df['species'] == pair[2]]['sequence'].tolist()[0]
    
    return ref_seq, target_seq
