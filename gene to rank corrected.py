import pandas as pd
def correct_gene_to_rank(df):
    reversed_sign = []
    for index, row in df.iterrows():
        res = row["res"]
        avg_positive = row["avg positive"]
        avg_negative = row["avg negative"]
        if res== True:
            avg_positive = -1*(avg_negative)
            reversed_sign.append(avg_positive)
        else:
            avg_positive = avg_positive
            reversed_sign.append(avg_positive)
    df["Corrected gene co-rel"] = reversed_sign 
    return df
