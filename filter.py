import pandas as pd
file = pd.read_csv("pairwise_gene_corelations_d=0.2.csv")
def add_filter(df):
    reverse = []
    for index, row in df.iterrows():
        num_negatives=row["num negative"]
        num_positives=row["num positive"]
        if num_negatives / (num_negatives +  num_positives) > 0.60:
            reverse.append(True)
        else:
            reverse.append(False)
            
    df["res"] = reverse
    return df
