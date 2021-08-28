from collections import Counter, OrderedDict
import numpy as np 
import scipy.stats
import gseapy as gp 

def enrichment_score(subset, gene_to_phenotype_correlation, p = 1, gene_to_rank=None):
    if gene_to_rank is None:
        gene_to_phenotype_correlation = collections.OrderedDict(gene_to_phenotype_correlation)
        gene_to_rank = {}
        rank = 1
        for k, v in gene_to_phenotype_correlation.items():
            gene_to_rank[k] = rank
            rank += 1
    scores = []
    N_R = 0
    for s in subset:
        N_R += abs(gene_to_phenotype_correlation[s])**p
    P_hit_i = []
    P_miss_i = []
    for i in range(1, len(gene_to_rank) + 1):
        P_hit = 0
        P_miss = 0
        for gene, j in gene_to_rank.items():
            r_j = gene_to_phenotype_correlation[gene]
            if j <= i:
                if gene in subset:
                    val = (abs(r_j)**p) / N_R
                    P_hit += val
                else:
                    val = 1 / (len(gene_to_rank) - len(subset))
                    P_miss += val
        P_hit_i.append(P_hit)
        P_miss_i.append(P_miss)
    
    ES = np.asarray(P_hit_i) - np.asarray(P_miss_i)
    argmax = np.argmax(list(map(lambda x : abs(x), ES)))
    score = ES[argmax]
    return score, ES
