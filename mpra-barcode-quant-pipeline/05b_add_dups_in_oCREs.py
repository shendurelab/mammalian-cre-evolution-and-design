#!/usr/bin/env python3
"""
Step 5b: Add duplicate CRE sequences to the count table.

Some oCRE sequences appear as duplicates in the reference (same sequence,
different tile IDs). These duplicates weren't recovered in sequencing because
they are identical, but their counts should be attributed to all duplicate
CRE IDs.

Usage:
    python 05b_add_dups_in_oCREs.py <count_table_wo_zeroes.txt.gz> <dup_ref.txt> <output_prefix>

Inputs:
    count_table_wo_zeroes.txt.gz : Count table from Step 5 (without zeroes)
    dup_ref.txt                  : Duplicate CRE reference file (tab-delimited:
                                   tile_name <tab> [list of duplicate tile names])

Outputs:
    {output_prefix}_duplicate_id_collapsed.txt                : CRE_id to quantification_id mapping
    {output_prefix}_count_table_wo_zeroes_w_dups.txt.gz       : Extended count table with duplicates
"""

import sys
import ast
import pandas as pd

if len(sys.argv) < 4:
    print("Usage: python 05b_add_dups_in_oCREs.py <count_table.txt.gz> <dup_ref.txt> <output_prefix>")
    sys.exit(1)

count_file = sys.argv[1]
dup_ref_file = sys.argv[2]
out_prefix = sys.argv[3]

print(f"Loading count table: {count_file}")
df = pd.read_csv(count_file, sep='\t', low_memory=False)
print(f"  Rows: {df.shape[0]}")

missing_dup_count = pd.DataFrame()
still_did_not_sequence = []

print(f"Reading duplicate reference: {dup_ref_file}")
with open(dup_ref_file, 'r') as fh:
    for line in fh:
        dup_read = line.strip().split('\t')
        if dup_read[0] == 'tile_name':
            continue
        literal_dup = ast.literal_eval(dup_read[1])
        dup_read = [dup_read[0]] + literal_dup
        temp_df = df[df['CRE_id'].isin(dup_read)]
        if temp_df.shape[0] == 0:
            print(f"  Did not sequence: {dup_read}")
            still_did_not_sequence += dup_read
        else:
            uniq_cre_id = temp_df['CRE_id'].unique()
            if len(uniq_cre_id) < len(dup_read):
                missing = [x for x in dup_read if x not in uniq_cre_id]
                for x in missing:
                    new_df = temp_df.copy()
                    new_df['CRE_id'] = x
                    new_df['quantification_id'] = uniq_cre_id[0]
                    missing_dup_count = pd.concat([missing_dup_count, new_df])

print(f"  Unsequenced duplicates: {len(still_did_not_sequence)}")

# Write duplicate ID mapping
if not missing_dup_count.empty:
    dup_id = missing_dup_count[['CRE_id', 'quantification_id']].drop_duplicates()
    dup_id.to_csv(f"{out_prefix}_duplicate_id_collapsed.txt", sep='\t', index=False)
    print(f"  Wrote: {out_prefix}_duplicate_id_collapsed.txt")

    # Prepare and concatenate
    missing_dup_count = missing_dup_count.drop(['quantification_id'], axis=1)
    missing_dup_count['is_duplicate'] = True
    df_w_dups = pd.concat([df, missing_dup_count])
else:
    df_w_dups = df

print(f"  Original rows: {df.shape[0]}")
print(f"  After adding duplicates: {df_w_dups.shape[0]}")

out_file = f"{out_prefix}_count_table_wo_zeroes_w_dups.txt.gz"
df_w_dups.to_csv(out_file, sep='\t', compression='gzip',
                  na_rep='NULL', index=False)
print(f"  Wrote: {out_file}")
