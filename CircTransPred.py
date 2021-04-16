# -*- coding: utf-8 -*-

import joblib
import argparse
import numpy as np

# =============================================================================
# calculate features from sequence

# ORF
# =============================================================================
import re

start_coden = 'ATG'
stop_coden = 'TAA|TAG|TGA'
ORF = '(ATG)([A-Z][A-Z][A-Z])*?(TAA|TAG|TGA)'

def check_orf(pro, len_seq):
    res = True
    
    #匹配到的ORF长度小于环状RNA长度
    if len(pro.group()) >= len_seq:
        res = False
    
    #匹配到的ORF起点应该在第一个环状RNA序列中
    if pro.start() >= len_seq:
        res = False
    
    #匹配到的ORF能翻译成的肽链长度应该大等于20个氨基酸
    if len(pro.group()) < 63:
        res = False
    
    return res

def call_ORF(seq):
    len_seq = len(seq)
    double = seq + seq
    p = re.compile(ORF, re.I)
    orf_len = 0
    for pro in p.finditer(double):
        if check_orf(pro, len_seq):
            orf_seq = pro.group()
            orf_len = max(orf_len, len(orf_seq))
    
    return orf_len, orf_len / len_seq

# =============================================================================

BASE = ('A', 'T', 'C', 'G')

# Fickett
# =============================================================================
def call_base_pos(base, pos):
    return max(pos[base])/(min(pos[base])+1)

def call_base_fre(base, seq):
    return seq.count(base) / len(seq)

def call_fickett(seq):
    pos = {'A': [0, 0, 0],
           'T': [0, 0, 0],
           'C': [0, 0, 0],
           'G': [0, 0, 0]}
    
    for i in range(len(seq)):
        if seq[i] == 'N':
            for base in pos:
                pos[base][i%3] += 1
        else:
            pos[seq[i]][i%3] += 1
    
    fickett = []
    for base in BASE:
        fickett.extend([call_base_pos(base, pos), call_base_fre(base, seq)])
    
    return fickett

# =============================================================================

# Hexamer
# =============================================================================

bimer = []
for b1 in BASE:
    for b2 in BASE:
        bimer.append(b1 + b2)

trimer = []
for b1 in BASE:
    for b2 in BASE:
        for b3 in BASE:
            trimer.append(b1 + b2 + b3)

hexamer = []
for aa1 in trimer:
    for aa2 in trimer:
        hexamer.append(aa1+aa2)

def _hexamer(seq):
    hex_dict = dict(zip(hexamer, [[0, 0, 0] for i in hexamer]))
    if len(seq) > 6:
        for i in range(len(seq) - 6):
            hex_dict[seq[i:i+6]][i%3] += 1
    return hex_dict

def call_hexamer(seq):
    hex_dict = _hexamer(seq)
    
#    hex_scores = dict(zip(hexamer, [0 for i in hexamer]))
#    for h in hexamer:
#        frame_nums = hex_dict.get(h)
#        f1, f2, f3 = frame_nums
#        #hex_scores[h] = (f1 + f2 + f3) / (3 * len(seq))
#        hex_scores[h] = (f1 + f2 + f3) / len(seq)
    
    hex_scores = []
    for h in hexamer:
        frame_nums = hex_dict.get(h)
        f1, f2, f3 = frame_nums
        hex_score = (f1 + f2 + f3) / len(seq)
        hex_scores.append(hex_score)
    
    return hex_scores

# =============================================================================

# =============================================================================
# load data and extract features from sequences

def get_fa(file):
    info, seq = [], {}
    
    with open(file, 'r') as f:
        for lines in f.readlines():
            if '>' in lines:
                name = lines[1:].strip()
                info.append(name)
                seq[name] = ""
            else:
                seq[name] += lines.strip()
    
    return info, seq

def get_fea(lst, fa_dict):
    
    tmp = []
    
    for i, name in enumerate(lst):
        seq = fa_dict[name]
        
        fea = []
        fea.extend(call_ORF(seq))
        fea.extend(call_fickett(seq))
        fea.extend(call_hexamer(seq))
        
        tmp.append(fea)
    
    data_fea = np.array(tmp)
    
    return data_fea

def get_features(fa_file):
    info, seq = get_fa(fa_file)
    fea = get_fea(info, seq)
    
    return info, fea

# =============================================================================

# =============================================================================

# =============================================================================
# predict the protein-coding potential of circular RNAs by using trained model

class Circ_Tri:
    
    def __init__(self, model_files):
        self.classifiers = [joblib.load(f) for f in model_files]
    
    def predict(self, x):
        pred = np.asarray([self.classifiers[i].predict(x) for i in range(3)])
        pred[0][pred[1] == pred[2]] = pred[1][pred[1] == pred[2]]
        return pred[0]
    
    def predict_proba(self, x):
        score1 = self.classifiers[0].predict_proba(x)[:, 1]
        score2 = self.classifiers[0].predict_proba(x)[:, 1]
        score3 = self.classifiers[0].predict_proba(x)[:, 1]
        score = [(score1[i] + score2[i] + score3[i]) / 3 for i in range(x.shape[0])]
        return score
        
def get_prediction(input_file, out_file):
    model_dir = "model/model{}.pkl"
    model_files = [model_dir.format(i) for i in range(1, 4)]
    
    tri = Circ_Tri(model_files)
    
    info, fea = get_features(input_file)
    
    pred_result = tri.predict(fea)
    pred_prob = tri.predict_proba(fea)
    
    with open(out_file, 'w') as fo:
        for i in range(len(info)):
            fo.write("{}\t{}\t{}\n".format(info[i], pred_result[i], pred_prob[i]))

# =============================================================================

def run_CircTransPred(parser):
    input_file, out_file = parser.input_file, parser.output_file
    
    get_prediction(input_file, out_file)

def parse_arguments(parser):
    
    parser.add_argument('--input_file', type = str, default='data/test.fa',
                        help = 'input FASTA file for predicting data')
    
    parser.add_argument('--output_file', type = str, default = 'data/pred_result.txt',
                        help = 'The output file used to store the prediction label and probability of input data')
    
    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'the identification of the protein-coding potential of circular RNAs')
    
    args = parse_arguments(parser)
    
    run_CircTransPred(args)