from collections import Counter, OrderedDict
import numpy as np 
import scipy.stats
import gseapy as gp 
def calculate_p_value(subset, gene_to_phenotype_correlation, iterations = 1000, p = 1, gene_to_rank = None):
    real_score, _ = enrichment_score(subset, gene_to_phenotype_correlation, p = p, gene_to_rank = gene_to_rank)
    scores = []
    
    if gene_to_rank is None:
        dict_to_use = gene_to_phenotype_correlation
    else:
        dict_to_use = gene_to_rank
        
    keys = list(dict_to_use.keys())
    
    for i in range(iterations):
        values = list(dict_to_use.values())
        random.shuffle(values)
        shuffled_dict = {}
        for k,v in zip(keys, values):
            shuffled_dict[k] = v
        score, _ = enrichment_score(subset, gene_to_phenotype_correlation, p=p, gene_to_rank = shuffled_dict)
        scores.append(score)
        
    if real_score > 0:
        scores = list(filter(lambda x : x > 0, scores))
        p = 1 - scipy.stats.percentileofscore(scores, real_score)
        print(p)
    else:
        scores = list(filter(lambda x : x < 0, scores))
        p = scipy.stats.percentileofscore(scores, real_score)
#         ks = scipy.stats.kstest(scores, real_score)
        
        c = Counter(scores)
        something = float(c[real_score]) / len(x)
        print(ks)
        print(p)
    
    return p 
