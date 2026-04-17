"""
Reproducible generation of extended-data figure 2c/d for paper 1.

Uses a trained endoderm chromBPNet model to compute DeepLIFT attribution
scores for the Epas1 (mouse vs African dog) and Gata4 (mouse vs rat) oCREs,
and plots sequence-logo tracks annotated with functional TFBS boxes.

Usage (from either the repo root or scripts/):
    python scripts/make_ext_data_fig2cd.py

Environment:
    conda activate /net/shendure/vol10/projects/tli/ledidi     # (tangermeme + bpnetlite + ledidi)

Inputs (all resolved relative to REPO_ROOT):
    data/ext_data_fig2cd/all_oCRE_cactusv2_reoriented_20240417_pwm.txt
    data/ext_data_fig2cd/tfbs_tables/Epas1_oCRE_func_TFBS_mapping.txt
    data/ext_data_fig2cd/tfbs_tables/Rattus_norvegicus_Gata4_CRE_func_TFBS_mapping.txt
    data/ext_data_fig2cd/chrombpnet_models/endo/endo_bias.h5
    data/ext_data_fig2cd/chrombpnet_models/endo/chrombpnet_nobias.h5
    data/mouse_TFBS_impact_and_affinity_CRM.txt
    data/motif_databases/CIS-BP_Mus_musculus.meme

Outputs (written to REPO_ROOT/figures/):
    ext_data_fig2c.pdf
    ext_data_fig2d.pdf
"""

import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn
import tensorflow as tf
import torch
from bpnetlite import ChromBPNet
from tangermeme.plot import plot_logo
from tangermeme.utils import random_one_hot
from tangermeme.ersatz import substitute
from tangermeme.deep_lift_shap import deep_lift_shap
from tangermeme.seqlet import recursive_seqlets
from tangermeme.io import read_meme
from tangermeme.annotate import annotate_seqlets

os.environ.setdefault('CUDA_VISIBLE_DEVICES', '0')
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

inputlen = 2114

# ---- Resolve repo paths robustly (runnable from REPO_ROOT or scripts/) ----
def find_repo_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in [here, here.parent, Path.cwd(), Path.cwd().parent]:
        if (candidate / "data").is_dir() and (candidate / "scripts").is_dir():
            return candidate.resolve()
    raise RuntimeError("Could not locate repo root (expected data/ and scripts/ directories).")

REPO_ROOT = Path(os.environ.get("REPO_ROOT", find_repo_root()))
DATA_DIR  = REPO_ROOT / "data"
FIG_DATA  = DATA_DIR / "ext_data_fig2cd"
OUT_DIR   = REPO_ROOT / "figures"
OUT_DIR.mkdir(exist_ok=True)

# Set global style parameters
plt.rcParams.update({
    # --- Base font settings ---
    'font.family': 'Helvetica',
    'font.size': 6,
    'axes.titlesize': 9,       # plot title
    'axes.labelsize': 7,       # axis titles
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'legend.title_fontsize': 6,

    # --- Line and axis appearance ---
    'axes.linewidth': 0.25,    # thin black axes lines
    'grid.linewidth': 0.25,
    'lines.linewidth': 0.5,
    'axes.edgecolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',

    # --- Strip/Facet style equivalents ---
    'axes.facecolor': 'white',
    'axes.labelcolor': 'black',
    'figure.facecolor': 'white',

    # --- Legend position and spacing ---
    'legend.loc': 'lower center',  # similar to ggplot bottom legend
    'legend.frameon': False,
    'legend.handlelength': 1.2,
    'legend.borderaxespad': 0.5,

    # --- Layout and padding ---
    'figure.autolayout': True,
    'axes.titlepad': 10,       # spacing below title
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5
})
seaborn.set_style("whitegrid", {'axes.edgecolor': 'black', 'axes.linewidth': 0.25})

class Wrapper(torch.nn.Module):
    def __init__(self, model):
        super().__init__()
        self.model = model 
    def forward(self, X):
        return self.model(X)[1]

final_TFs = ["Jun_Atf3","Foxa2","Gata4/6","Sox17","Klf4","Hnf1b"]
# --- add TFBS colored boxes ---
# define colors for each TF (or use a palette)
tf_colors = {
    "Gata4/6": "#EE7600",     # darkorange2
    "Sox17":  "#1E90FF",        # dodgerblue1
    "Jun_Atf3": "#EE3B3B",      # firebrick2
    "Foxa2": "#228B22",        # forestgreen
    "Klf4": "#B452CD",         # mediumorchid3
    "Hnf1b": "#FFC125"         # goldenrod1
}

# --- Read the file ---
mouse_DMS_tfbs = pd.read_csv(DATA_DIR / "mouse_TFBS_impact_and_affinity_CRM.txt", sep='\t')

# --- Add orientation numeric code ---
mouse_DMS_tfbs["TF_ori_pos"] = mouse_DMS_tfbs["TFBS_orientation"].apply(lambda x: 1 if x == '+' else -1)

# --- Combine Gata4 and Gata6 ---
gata_tfbs = mouse_DMS_tfbs[mouse_DMS_tfbs["TF_name"].isin(["Gata4", "Gata6"])].copy()

# Correct orientation issue for Gata6
gata_tfbs.loc[gata_tfbs["TF_name"] == "Gata6", "TF_ori_pos"] *= -1

# Recompute TFBS_orientation from TF_ori_pos
gata_tfbs["TFBS_orientation"] = gata_tfbs["TF_ori_pos"].apply(lambda x: '+' if x == 1 else '-')

# Merge Gata4 and Gata6 into "Gata4/6" and keep max affinity per (TFBS_start, TFBS_end, CRE)
gata_tfbs["TF_name"] = "Gata4/6"
gata_tfbs = (
    gata_tfbs.sort_values("norm_affinity", ascending=False)
    .groupby(["TFBS_start", "TFBS_end", "CRE"], as_index=False)
    .first()
)

# --- Replace Gata4 and Gata6 entries with merged Gata4/6 entries ---
mouse_DMS_tfbs = pd.concat([
    mouse_DMS_tfbs[~mouse_DMS_tfbs["TF_name"].isin(["Gata4", "Gata6"])],
    gata_tfbs
], ignore_index=True)

# Filter for final TF list
mouse_DMS_tfbs = mouse_DMS_tfbs[mouse_DMS_tfbs["TF_name"].isin(final_TFs)]

# --- Keep only functional TFs ---
mouse_func_tfbs = mouse_DMS_tfbs[mouse_DMS_tfbs["functional_TF"]]
    
oCRE_seq_dataframe = pd.read_csv(FIG_DATA / "all_oCRE_cactusv2_reoriented_20240417_pwm.txt",
                                 sep='\t', header=None)
oCRE_seq_dataframe[[ "CRE", "species" ]] = oCRE_seq_dataframe[0].str.split("__", expand=True)
oCRE_seq_dataframe.columns = ['id','seq','CRE','species']

mouse_seq_dataframe = oCRE_seq_dataframe[oCRE_seq_dataframe['species'] == 'Mus_musculus']

motifs = read_meme(str(DATA_DIR / "motif_databases" / "CIS-BP_Mus_musculus.meme"))
endo_TFs = ['Sox17','Sox4','Sox10','Sox13','Sox2','Sox9','Sox15','Sox6', ### sox family
            'Gata4','Gata1','Gata6','Gata2','Gata5', ### Gata family
            'Foxa2','Foxa3','Foxp1','Foxp2','Foxo1', ### Fox family
            'Klf4','Klf5','Klf6','Klf3','Klf1', ### Klf family
            'Hnf1b',
            'Jun','Atf3','Atf2','Atf1','Creb1','Fosl1','Fosl2','Fos','Fosb','Crem' #### Atf family
           ]
endo_motifs = {k: v for k, v in motifs.items() if any(tf in k for tf in endo_TFs)}
motif_names = list(endo_motifs.keys())
family_names = {'Gata4/6':['M00165_2.00 Gata6', 'M00167_2.00 Gata5',
                           'M08125_2.00 (Gata2)_(Homo_sapiens)_(DBD_1.00)', 'M08126_2.00 Gata4', 'M08127_2.00 Gata1'],
                'Hnf1b':['M00408_2.00 Hnf1b'],
               'Klf4':['M00242_2.00 (Klf11)_(Homo_sapiens)_(DBD_0.97)','M00787_2.00 Klf12',
                       'M02888_2.00 (Klf16)_(Homo_sapiens)_(DBD_0.94)','M02903_2.00 (Klf13)_(Homo_sapiens)_(DBD_0.99)', 
                       'M04400_2.00 (Klf6)_(Homo_sapiens)_(DBD_1.00)', 'M04418_2.00 (Klf3)_(Homo_sapiens)_(DBD_1.00)',
                      'M08086_2.00 (Klf5)_(Homo_sapiens)_(DBD_0.96)', 'M08099_2.00 Klf1',"M08021_2.00 Klf4",
                      'M08246_2.00 (Klf10)_(Homo_sapiens)_(DBD_0.99)', 'M08323_2.00 (Klf15)_(Homo_sapiens)_(DBD_1.00)',
                       'M08398_2.00 (Klf14)_(Homo_sapiens)_(DBD_0.97)'],
               'Sox17':['M08040_2.00 Sox17','M00826_2.00 Sox10','M00202_2.00 Sox15', 
                        'M00208_2.00 Sox21', 'M00212_2.00 Sox13', 'M00213_2.00 Sox4',
                       'M00826_2.00 Sox10','M03501_2.00 (Sox9)_(Homo_sapiens)_(DBD_1.00)','M08041_2.00 Sox2','M08169_2.00 Sox6'],
                'Foxa2':['M00792_2.00 Foxp2', 'M00793_2.00 Foxp1','M04831_2.00 (Foxa3)_(Homo_sapiens)_(DBD_1.00)',
                        'M08121_2.00 Foxa2', 'M08122_2.00 Foxo1'],
               'AP-1':['M01804_2.00 Atf1', 'M01805_2.00 Fosl1', 'M01806_2.00 Creb1', 
                       'M01807_2.00 Atf3', 'M01808_2.00 Atf2', 'M01809_2.00 Fosl2', 
                       'M01816_2.00 Jun', 'M01821_2.00 Crem', 'M01822_2.00 Jund',
                      'M04272_2.00 (Fosb)_(Homo_sapiens)_(DBD_1.00)','M08069_2.00 (Fos)_(Homo_sapiens)_(DBD_1.00)',
                      ]}
motif_to_family = {
    motif: family
    for family, motifs in family_names.items()
    for motif in motifs
}


# Initialize variables
bias_model = str(FIG_DATA / "chrombpnet_models" / "endo" / "endo_bias.h5")
chrombpnet_model = str(FIG_DATA / "chrombpnet_models" / "endo" / "chrombpnet_nobias.h5")

model = ChromBPNet.from_chrombpnet(bias_model, chrombpnet_model, 'endo')
model = model.to(device)
wrapper = Wrapper(model.accessibility)

# Generate one-hot encoded sequence
X_random = random_one_hot((100, 4, 2114), random_state = 42).type(torch.float32)

# predict mouse attr
mouse_EPAS1_seq = mouse_seq_dataframe[mouse_seq_dataframe['CRE'] == 'Epas1_chr17_10063']['seq'].tolist()[0]
X = substitute(X_random, mouse_EPAS1_seq)

# Define sequence parameters
seq_length = 2114
middle = seq_length // 2
mouse_string_length = len(mouse_EPAS1_seq)
half_length = mouse_string_length // 2
mouse_start = middle - half_length - 2
mouse_end = middle + half_length + 2

# DeepLIFT for Fine-tuned Model
mouse_attr = deep_lift_shap(wrapper, X, device='cuda', random_state=0)
mouse_attr_avg = mouse_attr.mean(dim=0, keepdim=True)
seqlets = recursive_seqlets(mouse_attr_avg.sum(dim = 1), threshold = 0.1, max_seqlet_len = 10)
motif_idxs, motif_pvalues = annotate_seqlets(X[0:1], seqlets, endo_motifs)
mouse_seqlet_annotations = seqlets.copy()
mouse_seqlet_annotations['example_idx'] = [motif_names[idx] for idx in motif_idxs]
mouse_seqlet_annotations['example_idx'] = mouse_seqlet_annotations['example_idx'].map(motif_to_family)
mouse_seqlet_annotations = mouse_seqlet_annotations[mouse_seqlet_annotations['attribution'] > 0.05]

# predict dog attr
dog_EPAS1_seq = oCRE_seq_dataframe[oCRE_seq_dataframe['CRE'] == 'Epas1_chr17_10063']
dog_EPAS1_seq = dog_EPAS1_seq[dog_EPAS1_seq['species'] == 'Lycaon_pictus']['seq'].tolist()[0]
X = substitute(X_random, dog_EPAS1_seq)

# Define sequence parameters
seq_length = 2114
middle = seq_length // 2
dog_string_length = len(dog_EPAS1_seq)
half_length = dog_string_length // 2
dog_start = middle - half_length - 2
dog_end = middle + half_length + 2

# DeepLIFT for Fine-tuned Model
dog_attr = deep_lift_shap(wrapper, X, device='cuda', random_state=0)
dog_attr_avg = dog_attr.mean(dim=0, keepdim=True)
seqlets = recursive_seqlets(dog_attr_avg.sum(dim = 1), threshold = 0.1, max_seqlet_len = 10)
motif_idxs, motif_pvalues = annotate_seqlets(X[0:1], seqlets, endo_motifs)
dog_seqlet_annotations = seqlets.copy()
dog_seqlet_annotations['example_idx'] = [motif_names[idx] for idx in motif_idxs]
dog_seqlet_annotations['example_idx'] = dog_seqlet_annotations['example_idx'].map(motif_to_family)
dog_seqlet_annotations = dog_seqlet_annotations[dog_seqlet_annotations['attribution'] > 0.05]

### get functional TFBS
mouse_epas1_tfbs = mouse_func_tfbs[mouse_func_tfbs['CRE'] == 'Epas1_chr17_10063']
epas1_oCRE_func_tfbs = pd.read_csv(FIG_DATA / "tfbs_tables" / "Epas1_oCRE_func_TFBS_mapping.txt", sep='\t')
dog_epas1_tfbs = epas1_oCRE_func_tfbs[epas1_oCRE_func_tfbs['species'] == 'Lycaon_pictus']

# Convert mm to inches (1 inch = 25.4 mm)
width_in = 210 / 25.4
height_in = 70 / 25.4

plt.clf()
plt.figure(figsize=(width_in, height_in))

# Subplot 1: mouse Epas1
ax1 = plt.subplot(2, 1, 1)  # 2 rows, 1 column, 1st plot
plot_logo(mouse_attr_avg[0], ax=ax1, start=mouse_start, end=mouse_end,) #annotations=mouse_seqlet_annotations, score_key='attribution')
ax1.set_title("Mus musculus Epas1:chr17_10063")
ax1.set_xlabel("")
ax1.set_ylabel("Attributions")
ax1.set_ylim(-0.1, 0.2)
ax1.grid(False)
ax1.set_xticks(range(0, mouse_string_length, 50))

# --- Add TFBS boxes ---
for _, row in mouse_epas1_tfbs.iterrows():
    start = row["TFBS_start"]
    end = row["TFBS_end"]
    tf = row["TF_name"]
    color = tf_colors.get(tf, "gray")
    rect = Rectangle(
        (start, -0.09),        # (x, y) position
        width=end - start,     # motif span
        height=0.02,           # bar height
        linewidth=0,
        facecolor=color,
        alpha=0.7
    )
    ax1.add_patch(rect)


# Subplot 2: african model
ax2 = plt.subplot(2, 1, 2)  # 3 rows, 1 column, 3rd plot
plot_logo(dog_attr_avg[0], ax=ax2, start=dog_start, end=dog_end, )#annotations=dog_seqlet_annotations, score_key='attribution')
ax2.set_title("Lycaon pictus (African dog) Epas1 CRE ortholog")
ax2.set_xlabel("bp position")
ax2.set_ylabel("Attributions")
ax2.set_ylim(-0.1, 0.2)
ax2.grid(False)
ax2.set_xticks(range(0, dog_string_length, 50))

# --- Add TFBS boxes ---
for _, row in dog_epas1_tfbs.iterrows():
    start = row["start_pos"]
    end = row["end_pos"]
    tf = row["TF_name"]
    color = tf_colors.get(tf, "gray")
    rect = Rectangle(
        (start, -0.09),        # (x, y) position
        width=end - start,     # motif span
        height=0.02,           # bar height
        linewidth=0,
        facecolor=color,
        alpha=0.7
    )
    ax2.add_patch(rect)

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig(OUT_DIR / "ext_data_fig2d.pdf", bbox_inches='tight')
 
#### now Gata4
### on mouse attr
mouse_GATA4_seq = mouse_seq_dataframe[mouse_seq_dataframe['CRE'] == 'Gata4_chr14_5729']['seq'].tolist()[0]
X = substitute(X_random, mouse_GATA4_seq)

# Define sequence parameters
seq_length = 2114
middle = seq_length // 2
mouse_string_length = len(mouse_GATA4_seq)
half_length = mouse_string_length // 2
mouse_start = middle - half_length - 2
mouse_end = middle + half_length + 2

# DeepLIFT for Fine-tuned Model
mouse_attr = deep_lift_shap(wrapper, X, device='cuda', random_state=0)
mouse_attr_avg = mouse_attr.mean(dim=0, keepdim=True)
seqlets = recursive_seqlets(mouse_attr_avg.sum(dim = 1), threshold = 0.1, max_seqlet_len = 10)
motif_idxs, motif_pvalues = annotate_seqlets(X[0:1], seqlets, endo_motifs)
mouse_seqlet_annotations = seqlets.copy()
mouse_seqlet_annotations['example_idx'] = [motif_names[idx] for idx in motif_idxs]
mouse_seqlet_annotations['example_idx'] = mouse_seqlet_annotations['example_idx'].map(motif_to_family)
mouse_seqlet_annotations = mouse_seqlet_annotations[mouse_seqlet_annotations['attribution'] > 0.05]

### on rat attr
rat_GATA4_seq = oCRE_seq_dataframe[oCRE_seq_dataframe['CRE'] == 'Gata4_chr14_5729']
rat_GATA4_seq = rat_GATA4_seq[rat_GATA4_seq['species'] == 'Rattus_norvegicus']['seq'].tolist()[0]
X = substitute(X_random, rat_GATA4_seq)

# Define sequence parameters
seq_length = 2114
middle = seq_length // 2
rat_string_length = len(rat_GATA4_seq)
half_length = rat_string_length // 2
rat_start = middle - half_length - 2
rat_end = middle + half_length + 2

# DeepLIFT for Fine-tuned Model
rat_attr = deep_lift_shap(wrapper, X, device='cuda', random_state=0)
rat_attr_avg = rat_attr.mean(dim=0, keepdim=True)
seqlets = recursive_seqlets(rat_attr_avg.sum(dim = 1), threshold = 0.1, max_seqlet_len = 10)
motif_idxs, motif_pvalues = annotate_seqlets(X[0:1], seqlets, endo_motifs)
rat_seqlet_annotations = seqlets.copy()
rat_seqlet_annotations['example_idx'] = [motif_names[idx] for idx in motif_idxs]
rat_seqlet_annotations['example_idx'] = rat_seqlet_annotations['example_idx'].map(motif_to_family)
rat_seqlet_annotations = rat_seqlet_annotations[rat_seqlet_annotations['attribution'] > 0.05]

### get functional TFBS
mouse_gata4_tfbs = mouse_func_tfbs[mouse_func_tfbs['CRE'] == 'Gata4_chr14_5729']
rat_gata4_tfbs = pd.read_csv(FIG_DATA / "tfbs_tables" / "Rattus_norvegicus_Gata4_CRE_func_TFBS_mapping.txt", sep='\t')

# Convert mm to inches (1 inch = 25.4 mm)
width_in = 210 / 25.4
height_in = 70 / 25.4

plt.clf()
plt.figure(figsize=(width_in, height_in))

# Subplot 1: mouse Epas1
ax1 = plt.subplot(2, 1, 1)  # 2 rows, 1 column, 1st plot
plot_logo(mouse_attr_avg[0], ax=ax1, start=mouse_start, end=mouse_end,) #annotations=mouse_seqlet_annotations, score_key='attribution')
ax1.set_title("Mus musculus Gata4:chr14_5729")
ax1.set_xlabel("")
ax1.set_ylabel("Attributions")
ax1.grid(False)
ax1.set_xticks(range(0, mouse_string_length, 50))
ax1.set_ylim(-0.1, 0.2)

# --- Add TFBS boxes ---
for _, row in mouse_gata4_tfbs.iterrows():
    start = row["TFBS_start"]
    end = row["TFBS_end"]
    tf = row["TF_name"]
    color = tf_colors.get(tf, "gray")
    rect = Rectangle(
        (start, -0.09),        # (x, y) position
        width=end - start,     # motif span
        height=0.02,           # bar height
        linewidth=0,
        facecolor=color,
        alpha=0.7
    )
    ax1.add_patch(rect)


# ax1.set_ylim(-0.005, 0.1)

# Subplot 2: sand rat Epas1
ax2 = plt.subplot(2, 1, 2)  # 2 rows, 1 column, 2nd plot
plot_logo(rat_attr_avg[0], ax=ax2, start=rat_start, end=rat_end,) #annotations=rat_seqlet_annotations, score_key='attribution')
ax2.set_title("Rattus norvegicus (Rat) Gata4 CRE ortholog")
ax2.set_xlabel("bp position")
ax2.set_ylabel("Attributions")
ax2.grid(False)
ax2.set_xticks(range(0, rat_string_length, 50))
ax2.set_ylim(-0.1, 0.2)

# --- Add TFBS boxes ---
for _, row in rat_gata4_tfbs.iterrows():
    start = row["start_pos"]
    end = row["end_pos"]
    tf = row["TF_name"]
    color = tf_colors.get(tf, "gray")
    rect = Rectangle(
        (start, -0.09),        # (x, y) position
        width=end - start,     # motif span
        height=0.02,           # bar height
        linewidth=0,
        facecolor=color,
        alpha=0.7
    )
    ax2.add_patch(rect)
# ax2.set_ylim(-0.005, 0.1)

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig(OUT_DIR / "ext_data_fig2c.pdf", bbox_inches='tight')

# Clear session to free memory
del model
del wrapper
tf.keras.backend.clear_session()