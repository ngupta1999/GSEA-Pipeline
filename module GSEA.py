from collections import Counter, OrderedDict
import numpy as np 
import scipy.stats
import gseapy as gp 
def module_gsea_with_corrections(
    modules,
    rnk, 
    out_path = "out_path/",
    corrections = None, 
    use_correction = False,
    min_module_size = 1,
    permutation_num = 100
):
    rnk = rnk.dropna()
    gene_members = set(rnk[0].values)
    for index, row in modules.iterrows():
        if use_correction and corrections is not None:
            module_values = corrections.loc[corrections["module"] == index]
            to_correct = module_values.loc[module_values["res"] == True]
            gene_sets = {}
            for index2, row2 in to_correct.iterrows():
                if row2["gene"] in gene_members:
                    gene_loc = rnk.loc[rnk[0] == row2["gene"]].iloc[0].name
                    rnk.at[gene_loc, 1] = -1 * rnk.loc[gene_loc][1]

        rnk = rnk.sort_values(by = 1)
        members = row["Members"].split(" ")
        members = list(filter(lambda x : x in gene_members, members))
        gene_sets[f"module_{index}"] = members
        if len(members) >= min_module_size:
            pre_res = gp.prerank(
                rnk=rnk,
                gene_sets = gene_sets, 
                processes=4,
                permutation_num = permutation_num,
                outdir= out_path + f'{index}/', 
                format='png', 
                min_size = min_module_size,
                seed=6
            )
        

modules = pd.read_csv("new_modules_d_0.2.csv", index_col = 0)
corrections = pd.read_csv("pairwise_gene_corelations_d=0.2.csv", index_col =0 )
corrections = add_filter(corrections)
rnk = pd.read_excel('0621_Hepg2_THP1__U937_lipidscreen_summary_Nilesh.xlsx')
#print (rnk.columns)
cols_keeping = ["id", "U937 TF uptake LFC MAGECK"]
cols = list(rnk.columns)
for c in cols_keeping:
    cols.remove(c)
rnk = rnk.drop(columns = cols)
rnk = rnk.rename(columns = {"id" : 0, "U937 TF uptake LFC MAGECK" : 1})
rnk = rnk.sort_values(by = 1)
module_gsea_with_corrections(modules, rnk, out_path = "U937 TF uptake LFC/", corrections = corrections, use_correction = True)
