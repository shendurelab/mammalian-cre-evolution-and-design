import os
import pandas as pd
from Bio import SeqIO
import random
import time
import subprocess
from SequencePremutationTools import SequencePremutationTools, dinuclShuffle

random.seed(42) ## set random seed

# all input files live in ./data/ next to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')

def create_fasta_list(fasta):
    fasta_list = list()
    for record in SeqIO.parse(fasta,'fasta'):
        fasta_list.append(str(record.seq))
    return fasta_list

def create_fasta(sequence_dict, output):
    with open(output, 'w') as fh:
        for header in sequence_dict:
            fh.write(">" + header +"\n")
            fh.write(sequence_dict[header]+"\n")

def create_pwm(sequence_dict, output):
    count = 0
    with open(output, 'w') as fh:
        for header in sequence_dict:
            seq = sequence_dict[header]
            fh.write(header +"\t" + seq+ "\ttile" +str(count) +"\n")
            count += 1

fasta_dict = dict()
for record in SeqIO.parse(os.path.join(DATA_DIR, 'DMS_300bp_tiles_15bp_symmetric_extension_max_activity_270bp_subtile_20231205.fa'),'fasta'):
    header= str(record.id)
    cre_id = header.split("_tile_")[1]
    fasta_dict[cre_id] = str(record.seq)


seq_dict = dict()
# create_fasta(fasta_dict, 'max_endo_chromBPnet_320bp_subtiles.fa')
endo_tfbs_hits = pd.read_csv(os.path.join(DATA_DIR, "TFBS_ProBound_Gata6_Sox17_Foxa2_Klf4_ParEndo_maxTileDMS_20240409.txt"), sep = '\t')
endo_tfbs_hits['start_pos'] = endo_tfbs_hits['TFBS_start'].astype(int) - 1
endo_tfbs_hits['end_pos'] = endo_tfbs_hits['TFBS_end'].astype(int)
endo_tfbs_hits = endo_tfbs_hits[endo_tfbs_hits['CRE'].isin(fasta_dict.keys())]

endo_tfbs_minus1_hits = pd.read_csv(os.path.join(DATA_DIR, "TFBS_flanks_minus1_ProBound_Gata6_Sox17_Foxa2_Klf4_ParEndo_maxTileDMS_20240409.txt"), sep = '\t')
endo_tfbs_minus1_hits['start_pos'] = endo_tfbs_minus1_hits['TFBS_start'].astype(int) - 1
endo_tfbs_minus1_hits['end_pos'] = endo_tfbs_minus1_hits['TFBS_end'].astype(int)
endo_tfbs_minus1_hits = endo_tfbs_minus1_hits[endo_tfbs_minus1_hits['CRE'].isin(fasta_dict.keys())]

endo_tfbs_minus2_hits = pd.read_csv(os.path.join(DATA_DIR, "TFBS_flanks_minus2_ProBound_Gata6_Sox17_Foxa2_Klf4_ParEndo_maxTileDMS_20240409.txt"), sep = '\t')
endo_tfbs_minus2_hits['start_pos'] = endo_tfbs_minus2_hits['TFBS_start'].astype(int) - 1
endo_tfbs_minus2_hits['end_pos'] = endo_tfbs_minus2_hits['TFBS_end'].astype(int)
endo_tfbs_minus2_hits = endo_tfbs_minus2_hits[endo_tfbs_minus2_hits['CRE'].isin(fasta_dict.keys())]

### use old background sequences from pCRE v1
bg_df = pd.read_csv(os.path.join(DATA_DIR, "max_endo_chromBPnet_320bp_subtile_shuffled_v2.tsv"), sep = '\t', header = None)
bg_df = bg_df[bg_df[0].str.startswith('background')]
bg_full_list = [x.lower() for x in bg_df[1].tolist()]
bg_list = [x[10:310].lower() for x in bg_full_list]

### add 100 new ones
bg_seq_list = create_fasta_list(os.path.join(DATA_DIR, "endo_bg_2k.fa"))
bg_seq_list = [x for x in bg_seq_list if x[:320].lower() not in bg_full_list]
new_bg_list = [x for x in bg_seq_list if len(x) >= 300]
new_bg_list = random.sample(new_bg_list, k=100)
new_bg_list = [x[:300].lower() for x in new_bg_list]
bg_list = bg_list + new_bg_list
shuffle_bg_list = [dinuclShuffle(x.upper()).lower() for x in bg_list]

print(len(bg_list))

select_bg_tiles = ['background_seq|background|8','background_seq|background|28',
                  'background_seq|background|23','background_seq|background|90',
                  'background_seq|background|19','background_seq|background|25']
select_bg_dict = {'background_seq|background|8':'permissive','background_seq|background|28':'permissive',
                  'background_seq|background|23':'intermediate','background_seq|background|90':'intermediate',
                  'background_seq|background|19':'recalcitrant','background_seq|background|25':'recalcitrant'}

select_bg_df = bg_df[bg_df[0].isin(select_bg_tiles)]
select_bg_list = select_bg_df[1].tolist()
select_bg_list = [x[10:310].lower() for x in select_bg_list]
print(len(select_bg_list))

order_dict = dict()

for cre in fasta_dict:
    print(cre)
    check_dict = dict()
    template_seq = fasta_dict[cre].upper()
    seq_dict[cre + "|original_seq"] = template_seq
    temp_tfbs = endo_tfbs_hits[endo_tfbs_hits['CRE'] == cre]
    temp_tfbs = temp_tfbs.sort_values(by = ['start_pos','end_pos'])
    spt = SequencePremutationTools(template_seq, temp_tfbs, len(template_seq))
    spt.process_intervals()
    print('\tMaking dishuffles')
    for i in range(0, 200):
        shuffle_seq = spt.dinuclShuffle_seq()
        shuffle_flank_seq = spt.shuffleFlankingTFBS()
        shuffle_TFBS_seq = spt.shuffleTFBS()
        seq_dict[cre + "|dishuffle_all|shuffle_" + str(i)] = shuffle_seq
        seq_dict[cre + "|dishuffle_TFBS_full|shuffle_" + str(i)] = shuffle_TFBS_seq
    
    print('\tMaking thripsis')
    for i in range(0,200):
        shuffle_segment_break5 = spt.shuffle_RandomSegments(num_cuts = 5,random_seed = i)
        seq_dict[cre + "|CRE_thripsis_5breaks_TFBS_full|ss_" + str(i)] = shuffle_segment_break5
        shuffle_segment_break10 = spt.shuffle_RandomSegments(num_cuts = 10, random_seed = i)
        seq_dict[cre + "|CRE_thripsis_10breaks_TFBS_full|ss_" + str(i)] = shuffle_segment_break10
        shuffle_segment_break20 = spt.shuffle_RandomSegments(num_cuts = 20, random_seed = i)
        seq_dict[cre + "|CRE_thripsis_20breaks_TFBS_full|ss_" + str(i)] = shuffle_segment_break20
    
    print('\tMaking reconstitutions and depositions')
    for i in range(0, len(bg_list)):
        bg_seq = bg_list[i].lower()
        shuffle_bg_seq = shuffle_bg_list[i].lower()
        TFBS_implant_seq = spt.TFBS_implantation(bg_seq)
        seq_dict[cre + "|reconstitution_TFBS_full|bg_" + str(i + 1)] = TFBS_implant_seq
        check_dict[cre + "|reconstitution_TFBS_full|bg_" + str(i + 1)] = TFBS_implant_seq
        shuffle_TFBS_implant_seq = spt.shuffle_TFBS_implantation(bg_seq, 101)
        seq_dict[cre + "|random_deposition_all_bkg_TFBS_full|bg_" + str(i + 1) + "|fixed_deposition_1"] = shuffle_TFBS_implant_seq
        check_dict[cre + "|random_deposition_all_bkg_TFBS_full|bg_" + str(i + 1) + "|fixed_deposition_1"] = shuffle_TFBS_implant_seq
        shuffle_TFBS_implant_seq = spt.shuffle_TFBS_implantation(bg_seq, 102)
        seq_dict[cre + "|random_deposition_all_bkg_TFBS_full|bg_" + str(i + 1) + "|fixed_deposition_2"] = shuffle_TFBS_implant_seq
        check_dict[cre + "|random_deposition_all_bkg_TFBS_full|bg_" + str(i + 1) + "|fixed_deposition_2"] = shuffle_TFBS_implant_seq
        ### make shuffled background reconstitutions and depostions
        TFBS_implant_seq = spt.TFBS_implantation(shuffle_bg_seq)
        seq_dict[cre + "|reconstitution_TFBS_full|shuffled_bg_" + str(i + 1)] = TFBS_implant_seq
        shuffle_TFBS_implant_seq = spt.shuffle_TFBS_implantation(shuffle_bg_seq, 101)
        seq_dict[cre + "|random_deposition_all_bkg_TFBS_full|shuffled_bg_" + str(i + 1)+ "|fixed_deposition_1"] = shuffle_TFBS_implant_seq
        check_dict[cre + "|random_deposition_all_bkg_TFBS_full|shuffled_bg_" + str(i + 1)+ "|fixed_deposition_1"] = shuffle_TFBS_implant_seq
        shuffle_TFBS_implant_seq = spt.shuffle_TFBS_implantation(shuffle_bg_seq, 102)
        seq_dict[cre + "|random_deposition_all_bkg_TFBS_full|shuffled_bg_" + str(i + 1)+ "|fixed_deposition_2"] = shuffle_TFBS_implant_seq
        check_dict[cre + "|random_deposition_all_bkg_TFBS_full|shuffled_bg_" + str(i + 1)+ "|fixed_deposition_2"] = shuffle_TFBS_implant_seq
    
    print('\tMaking more depositions')
    for i in range(0, len(select_bg_list)):
        bg_seq = select_bg_list[i].lower()
        bg_type = select_bg_dict[select_bg_tiles[i]]
        for seed in range(0, 100):
#             name = cre + "|random_deposition_fixed_bkg_TFBS_full|"+ bg_type +"_bg_" + str(i + 1) +"|deposition_" + str(seed+1)
#             if name == 'Bend5_chr4_8201|random_deposition_fixed_bkg_TFBS_full|permissive_bg_1|deposition_2':
            shuffle_TFBS_implant_seq = spt.shuffle_TFBS_implantation(bg_seq, seed)
#                 print(shuffle_TFBS_implant_seq)
#                 print(bg_type)
            seq_dict[cre + "|random_deposition_fixed_bkg_TFBS_full|"+ bg_type +"_bg_" + str(i + 1) +"|deposition_" + str(seed+1)] = shuffle_TFBS_implant_seq
            check_dict[cre + "|random_deposition_fixed_bkg_TFBS_full|"+ bg_type +"_bg_" + str(i + 1) +"|deposition_" + str(seed+1)] = shuffle_TFBS_implant_seq

    print('\tMaking +1 flanking reconstitutions')
    ### minus 1 flank
    temp_tfbs = endo_tfbs_minus1_hits[endo_tfbs_minus1_hits['CRE'] == cre]
    temp_tfbs = temp_tfbs.sort_values(by = ['start_pos','end_pos'])
    spt = SequencePremutationTools(template_seq, temp_tfbs, len(template_seq))
    spt.process_intervals()
            
    for i in range(0, len(bg_list)):
        bg_seq = bg_list[i].lower()
        shuffle_bg_seq = shuffle_bg_list[i].lower()
        TFBS_implant_seq = spt.TFBS_implantation(bg_seq)
        seq_dict[cre + "|reconstitution_TFBS_minus1_flank|bg_" + str(i + 1)] = TFBS_implant_seq
        ### make shuffled background reconstitutions
        TFBS_implant_seq = spt.TFBS_implantation(shuffle_bg_seq)
        seq_dict[cre + "|reconstitution_TFBS_minus1_flank|shuffled_bg_" + str(i + 1)] = TFBS_implant_seq
        
    print('\tMaking +2 flanking reconstitutions')
    ### minus 2 flank
    temp_tfbs = endo_tfbs_minus2_hits[endo_tfbs_minus2_hits['CRE'] == cre]
    temp_tfbs = temp_tfbs.sort_values(by = ['start_pos','end_pos'])
    spt = SequencePremutationTools(template_seq, temp_tfbs, len(template_seq))
    spt.process_intervals()
        
    for i in range(0, len(bg_list)):
        bg_seq = bg_list[i].lower()
        shuffle_bg_seq = shuffle_bg_list[i].lower()
        TFBS_implant_seq = spt.TFBS_implantation(bg_seq)
        seq_dict[cre + "|reconstitution_TFBS_minus2_flank|bg_" + str(i + 1)] = TFBS_implant_seq
        ### make shuffled background reconstitutions
        TFBS_implant_seq = spt.TFBS_implantation(shuffle_bg_seq)
        seq_dict[cre + "|reconstitution_TFBS_minus2_flank|shuffled_bg_" + str(i + 1)] = TFBS_implant_seq
        
    
    for seq_name in check_dict:
        seq = check_dict[seq_name]
        for tfbs in temp_tfbs['substring'].tolist():
            if tfbs not in seq:
                print('error missing tfbs: {}'.format(tfbs))
                print(seq_name)
                print(seq)
    

for i in range(0,len(bg_list)):
    seq_dict["bg_" + str(i + 1) + "|background_seq"] = bg_list[i].lower()
    seq_dict["shuffled_bg_" + str(i + 1) + "|shuffled_background_seq"] = shuffle_bg_list[i].lower()

create_pwm(seq_dict, 'JB_tiles_300bp_subtile_shuffled_v4.tsv')
