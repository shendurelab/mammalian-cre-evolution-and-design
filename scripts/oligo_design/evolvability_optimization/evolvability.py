import numpy as np
##################################################################
# Computational setup
##################################################################
#Modules
import os
import json
import sys
import pickle
import random
import pandas as pd
import numpy as np
import pyfaidx
import tensorflow as tf
from tqdm import tqdm
from itertools import permutations
import chrombpnet.training.utils.losses as losses
import chrombpnet.training.utils.one_hot as one_hot
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
from collections import ChainMap
import deepdish as dd

#Custom functions
sys.path.insert(0, f'/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/util_scripts/Brennan_Zelda_2023/analysis/scripts/py')
from chrombpnet_predict_functions import predict_chrombpnet
pd.options.mode.chained_assignment = None  # Disable the warning

#Independent Variables
BASE_DIR='/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet'
os.chdir(BASE_DIR)
INPUT_SEQLEN=2114
OUTPUT_SEQLEN=1000
trials = 1000 # 10
primary_position = 300
measurement_window_around_motif = 200
step  = 1

genome_fasta = pyfaidx.Fasta('/net/gs/vol1/home/tli824/bin/references/mouse/UCSC/mm10/mm10.fa')
#Dependent variables
MODEL_DIR=BASE_DIR + '/chrombpnet_models/endo/models/'


np.random.seed(42)

def get_seq(peaks_df, genome, width):
    """
    Same as get_cts, but fetches sequence from a given genome.
    """
    vals = dict()

    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
        vals[r['name']] = sequence.upper()
        
    return vals

def load_model_wrapper(model_path):
    # read .h5 model
    custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_path, compile=False)
    print("got the model")
    model.summary()
    return model

def generate_in_silico_mutatagenesis(sequence, start = 0, end = 2114):
    """
    Generate all possible single mutations in a given sequence at a certain location (default all)
    """
    mutations = dict()
    pos = 0
    
    for i in range(start, end):
        # Extract the base at the current position
        base = sequence[i]

        # Define the three alternative bases
        alternatives = ['A', 'C', 'G', 'T']
        alternatives.remove(base)  # Remove the original base from alternatives

        # Generate mutations by substituting the base at the current position
        for alt_base in alternatives:
            mutation = sequence[:i] + alt_base + sequence[i+1:]
            mutations[alt_base + '_' + str(i)] = mutation
        pos += 1
        
    mutations = pd.DataFrame(mutations.items(), columns=['pos', 'seq'])
    return mutations

#########################################################################################################
#########################################################################################################
### Find the single bp mutation that has highest expression 
def maximize_next_generation(population_current, population_current_fitness, model):
    population_one_hot = one_hot.dna_to_one_hot(population_current)
    
    pop_pf_fitness, pop_ct_fitness = predict_chrombpnet(model = model,
                                           seqs = population_one_hot)
    try:
        pop_ct_fitness = [np.log2(np.sum(x) + 1) for x in pop_ct_fitness]
        max_pop_fitness_idx = np.argmax(pop_ct_fitness)

        if pop_ct_fitness[max_pop_fitness_idx] > population_current_fitness:
            return population_current[max_pop_fitness_idx], max_pop_fitness_idx, pop_ct_fitness[max_pop_fitness_idx]
        else:
            return None, None, population_current_fitness
    except Exception as e:
        print(e)
        return None, None, population_current_fitness
    
#########################################################################################################
#########################################################################################################
#########################################################################################################


#########################################################################################################
#########################################################################################################
### INTRA FUNCTION NAMES ARE WRONG, VARIABLE IS RIGHT THOUGH , Find the single bp mutation that has least expression 
def minimize_next_generation(population_current, population_current_fitness, model):
    population_one_hot = one_hot.dna_to_one_hot(population_current)
    
    pop_pf_fitness, pop_ct_fitness = predict_chrombpnet(model = model,
                                           seqs = population_one_hot)
    try:
        pop_ct_fitness = [np.log2(np.sum(x) + 1) for x in pop_ct_fitness]
        min_pop_fitness_idx = np.argmin(pop_ct_fitness)

        if pop_ct_fitness[min_pop_fitness_idx] < population_current_fitness:
            return population_current[min_pop_fitness_idx], min_pop_fitness_idx, pop_ct_fitness[min_pop_fitness_idx]
        else:
            return None, None, population_current_fitness
    except Exception as e:
        print(e)
        return None, None, population_current_fitness

#########################################################################################################
#########################################################################################################
#########################################################################################################



#########################################################################################################
#########################################################################################################
### randomly introduce n mutations in a given sequence space S
### Simulate random drift
def neutral_next_generation(sequence, start = 1, end = 2114, n = 1): 
    """
    given a sequence space (S) with given start and end (default entire sequence space)
    introduce (n) number of mutations in the sequence space
    """
    seq_space = [x for x in range(start - 1, end)]
    if n > len(seq_space):
        print('Error number of mutations greater than sequence space')
        return None
    else:
        random_mutated_seq = sequence
        # Generate random indices to select n elements
        random_bases = np.random.choice(seq_space, n, replace=False)

        for i in random_bases:
            # Extract the base at the current position
            base = random_mutated_seq[i]

            # Define the three alternative bases
            alternatives = ['A', 'C', 'G', 'T']
            alternatives.remove(base)  # Remove the original base from alternatives
            alt_base = np.random.choice(alternatives)

            random_mutated_seq = random_mutated_seq[:i] + alt_base + random_mutated_seq[i+1:]
    
        return random_mutated_seq
    
#########################################################################################################
#########################################################################################################
#########################################################################################################


