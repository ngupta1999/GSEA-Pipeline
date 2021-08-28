import pandas as pd
def get_ordinal_gene_to_phenotype_correlation(ordered_genes):
    to_return = {}
    rank = 1
    step_size = 1 / len(ordered_genes)
    for gene in ordered_genes:
        corr = 1 - (step_size*rank)
        to_return[gene] = corr
        rank += 1
    return to_return 
  
 def get_gene_to_rank(ordered_genes):
    to_return = {}
    rank = 1
    for gene in ordered_genes:
        to_return[gene] = rank
        rank += 1
    return to_return 
  
def get_constant_gene_rank(genes):
    to_return = {k : 1 for k in genes}
    return to_return 
