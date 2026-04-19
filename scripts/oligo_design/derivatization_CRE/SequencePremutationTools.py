import sys,string,random
import pandas as pd
import numpy as np
from Bio.Seq import Seq

random.seed(42) ## set random seed

# altschulEriksonDinuclShuffle.py
# P. Clote, Oct 2003
# NOTE: One cannot use function "count(s,word)" to count the number
# of occurrences of dinucleotide word in string s, since the built-in
# function counts only nonoverlapping words, presumably in a left to
# right fashion.
def computeCountAndLists(s):
    #WARNING: Use of function count(s,'UU') returns 1 on word UUU
    #since it apparently counts only nonoverlapping words UU
    #For this reason, we work with the indices.

    #Initialize lists and mono- and dinucleotide dictionaries
    List = {} #List is a dictionary of lists
    List['A'] = []; List['C'] = [];
    List['G'] = []; List['T'] = [];
    nuclList = ["A","C","G","T"]
    s = s.upper()
    s = s.replace("T","T")
    nuclCnt = {}  #empty dictionary
    dinuclCnt = {}  #empty dictionary
    for x in nuclList:
        nuclCnt[x]=0
        dinuclCnt[x]={}
        for y in nuclList:
            dinuclCnt[x][y]=0

    #Compute count and lists
    nuclCnt[s[0]] = 1
    nuclTotal = 1
    dinuclTotal = 0
    for i in range(len(s)-1):
        x = s[i]; y = s[i+1]
        List[x].append( y )
        nuclCnt[y] += 1; nuclTotal  += 1
        dinuclCnt[x][y] += 1; dinuclTotal += 1
    assert (nuclTotal==len(s))
    assert (dinuclTotal==len(s)-1)
    return (nuclCnt,dinuclCnt,List)
 
def chooseEdge(x,dinuclCnt):
    numInList = 0
    for y in ['A','C','G','T']:
        numInList += dinuclCnt[x][y]
    z = random.random()
    denom=dinuclCnt[x]['A']+dinuclCnt[x]['C']+dinuclCnt[x]['G']+dinuclCnt[x]['T']
    numerator = dinuclCnt[x]['A']
    if z < float(numerator)/float(denom):
        dinuclCnt[x]['A'] -= 1
        return 'A'
    numerator += dinuclCnt[x]['C']
    if z < float(numerator)/float(denom):
        dinuclCnt[x]['C'] -= 1
        return 'C'
    numerator += dinuclCnt[x]['G']
    if z < float(numerator)/float(denom):
        dinuclCnt[x]['G'] -= 1
        return 'G'
    dinuclCnt[x]['T'] -= 1
    return 'T'


def connectedToLast(edgeList,nuclList,lastCh):
    D = {}
    for x in nuclList: D[x]=0
    for edge in edgeList:
        a = edge[0]; b = edge[1]
        if b==lastCh: D[a]=1
    for i in range(2):
        for edge in edgeList:
            a = edge[0]; b = edge[1]
            if D[b]==1: D[a]=1
    ok = 0
    for x in nuclList:
        if x!=lastCh and D[x]==0: return 0
    return 1
 

def eulerian(s):
    nuclCnt,dinuclCnt,List = computeCountAndLists(s)
    #compute nucleotides appearing in s
    nuclList = []
    for x in ["A","C","G","T"]:
        if x in s: nuclList.append(x)
    #compute numInList[x] = number of dinucleotides beginning with x
    numInList = {}
    for x in nuclList:
        numInList[x]=0
        for y in nuclList:
            numInList[x] += dinuclCnt[x][y]
    #create dinucleotide shuffle L 
    firstCh = s[0]  #start with first letter of s
    lastCh  = s[-1]
    edgeList = []
    for x in nuclList:
        if x!= lastCh: edgeList.append( [x,chooseEdge(x,dinuclCnt)] )
    ok = connectedToLast(edgeList,nuclList,lastCh)
    return (ok,edgeList,nuclList,lastCh)


def shuffleEdgeList(L):
    n = len(L); barrier = n
    for i in range(n-1):
        z = int(random.random() * barrier)
        tmp = L[z]
        L[z]= L[barrier-1]
        L[barrier-1] = tmp
        barrier -= 1
    return L


def dinuclShuffle(s):
    ok = 0
    while not ok:
        ok,edgeList,nuclList,lastCh = eulerian(s)
    nuclCnt,dinuclCnt,List = computeCountAndLists(s)

    #remove last edges from each vertex list, shuffle, then add back
    #the removed edges at end of vertex lists.
    for [x,y] in edgeList: List[x].remove(y)
    for x in nuclList: shuffleEdgeList(List[x])
    for [x,y] in edgeList: List[x].append(y)

    #construct the eulerian path
    L = [s[0]]; prevCh = s[0]
    for i in range(len(s)-2):
        ch = List[prevCh][0] 
        L.append( ch )
        del List[prevCh][0]
        prevCh = ch
    L.append(s[-1])
    t = "".join(L)
    return t


class SequencePremutationTools:
    def __init__(self, seq, tfbs_hits, max_length):
        self.seq = seq
        self.max_length = max_length
        self.tfbs_hits = tfbs_hits
        self.tfbs_coords = None
        self.gap_coords = None
        self.total_intervals = None
        self.tfbs_dict = dict()
    
    def dinuclShuffle_seq(self):
        return dinuclShuffle(self.seq)
    
    def process_intervals(self, add_bp = None):
        if add_bp:
            self.addTFBS_flanksites(add_bp)
        self.mergeIntervals()
        self.fill_gaps()
        self.total_intervals = sorted(self.tfbs_coords + self.gap_coords, key=lambda l:l[0])
    
    def fill_gaps(self):
        ### get gaps from coordinates
        last_x, last_y = 0, 0
        gaps = []

        for x, y in self.tfbs_coords:
            if x == 0:
                last_x = x
                last_y = y
                continue
            if len(gaps) == 0:
                if last_y > 0:
                    gaps.append([last_y, x])
                else:
                    gaps.append([last_x, x])
                last_x = x
                last_y = y
            else:
                if x == last_y:
                    last_x = x
                    last_y = y
                else:
                    if last_y > x:
                        last_x = x
                    else:
                        gaps.append([last_y, x])
                        last_x = x
                        last_y = y
        last_coord = self.tfbs_coords[-1]
        if len(gaps) > 0:
            last_gap = gaps[-1]
            if last_coord[-1] != self.max_length:
                if last_gap[-1] != self.max_length:
                    if last_coord[-1] > last_gap[-1]:
                        gaps.append([last_coord[-1], self.max_length])
                    else:
                        gaps.append([last_gap[-1], self.max_length])
        self.gap_coords = gaps

    def mergeIntervals(self):
        ### merge overlapping intervals (specifically TFBS regions)
        self.tfbs_dict = dict()
        intervals = [[x,y] for x, y in zip(self.tfbs_hits['start_pos'],self.tfbs_hits['end_pos'])]
        # Sort the array on the basis of start values of intervals.
        intervals.sort()
        tfbs_coords = []
        # insert first interval into stack
        tfbs_coords.append(intervals[0])
        for i in intervals[1:]:
            # Check for overlapping interval,
            # if interval overlap
            if tfbs_coords[-1][0] <= i[0] <= tfbs_coords[-1][-1]:
                tfbs_coords[-1][-1] = max(tfbs_coords[-1][-1], i[-1])
            else:
                tfbs_coords.append(i)
        self.tfbs_coords = tfbs_coords
        for i in self.tfbs_coords:
            self.tfbs_dict[tuple(i)] = self.seq[i[0]:i[1]]
        
    def shuffleFlankingTFBS(self):
        seq_parts = []
        ### shuffle regions next to TFBS 
        for i in self.total_intervals:
            segment = self.seq[i[0]:i[1]]
            if len(segment) == 1:
                seq_parts.append(segment)
            else:
                if i in self.tfbs_coords:
                    seq_parts.append(segment)
                else:
                    seq_parts.append(dinuclShuffle(segment))
        shuffled_seq = ''.join(seq_parts).strip()
        return(shuffled_seq)

    def shuffleTFBS(self):
        ### shuffle TFBS regions
        seq_parts = []
        ### shuffle regions next to TFBS 
        for i in self.total_intervals:
            segment = self.seq[i[0]:i[1]]
            if len(segment) == 1:
                seq_parts.append(segment)
            else:
                if i in self.tfbs_coords:
                    seq_parts.append(dinuclShuffle(segment))
                else:
                    seq_parts.append(segment)
        shuffled_seq = ''.join(seq_parts).strip()
        return(shuffled_seq)
    
    def shuffle_RandomSegments(self, num_cuts = 10, reverse_num = None, random_seed = 0):
        np.random.seed(random_seed)
        ### figure out delta difference where intervals have to be a certain bp size
        all_possible_selections = [] ### any breakpoints in the seqeunce outside the TFBS
        for i in self.gap_coords:
            all_possible_selections += [x for x in range(i[0],i[-1] + 1)]
        all_possible_selections = np.unique(all_possible_selections)
        cut_sites = np.sort(np.random.choice(all_possible_selections, num_cuts))
        intervals = []
        last_x = 0
        for i in range(0, len(cut_sites)):
            if len(intervals) == 0:
                intervals.append([last_x, cut_sites[i]])
                last_x = cut_sites[i]
            else:
                intervals.append([last_x, cut_sites[i]])
                last_x = cut_sites[i]
        if cut_sites[-1] != 500:
            intervals.append([last_x, 500])
        reverse_parts = []
        if reverse_num:
            reverse_parts = [x for x in range(0,num_cuts)]
            reverse_parts = np.sort(np.random.choice(reverse_parts, reverse_num))
        
        seq_parts = []
        ### shuffle regions next to TFBS 
        for i in range(0, len(intervals)):
            rg = intervals[i]
            segment = self.seq[rg[0]:rg[1]]
            if len(reverse_parts) > 0: ### reverse comp parts of the segments
                if i in reverse_parts:
                    segment = Seq(segment)
                    segment = str(segment.reverse_complement())
            seq_parts.append(segment)
        
        random.seed(random_seed)
        random.shuffle(seq_parts) ### shuffle dna parts
        shuffled_segments = ''.join(seq_parts).strip()
        return shuffled_segments
        
    def addTFBS_flanksites(self, add_bp = 1):
        ### add flanking regions to TFBS 
        intervals = [[x,y] for x, y in zip(self.tfbs_hits['start_pos'],self.tfbs_hits['end_pos'])]
        new_intervals = []
        for i in intervals:
            start = i[0]
            end = i[-1]
            if start - add_bp < 0:
                start = 0 ## cap at minimum sequence length
            else:
                start = start - add_bp
            if end + add_bp > self.max_length:
                end = self.max_length ## cap at maximum sequence length
            else:
                end = end + add_bp
            new_intervals.append([start,end])
        ### update tfbs_hits with extended flanking position   
        self.tfbs_hits[['start_pos','end_pos']] = new_intervals
    
    def TFBS_implantation(self, seq):
        seq = seq[:self.max_length]
        for key in self.tfbs_dict.keys():
            seq = seq[:key[0]] + self.tfbs_dict[key] + seq[key[1]:]
        return seq
    
    def shuffle_TFBS_implantation(self, seq, random_seed = 0):
        random.seed(random_seed)
        seq = seq[:self.max_length]
        # initialize the new dictionary
        new_tfbs_dict = {}
        used_parts = set()
        # iterate until we have created the desired number of intervals
        for intervals in self.tfbs_dict.keys():
            substr = self.tfbs_dict[intervals]
            overlap = False
            start = None
            end = None
            while not overlap:
                random_start = random.randint(1, self.max_length - len(substr))
                random_end = random_start + len(substr)
                if not any(random_start <= e <= random_end for e in used_parts):
                    for i in range(random_start, random_end + 1):
                        used_parts.add(i)
                    start = random_start
                    end = random_end
                    overlap = True
            new_tfbs_dict[(start,end)] = substr
#         print(self.tfbs_dict)
#         print(new_tfbs_dict)
        for key in new_tfbs_dict.keys():
            seq = seq[:key[0]] + new_tfbs_dict[key] + seq[key[1]:]
        return seq