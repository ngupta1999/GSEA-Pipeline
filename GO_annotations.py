import pandas as pd
data1 = pd.read_excel("GSEA_GO_Annotated_Modules_d_0.2.xlsx", index_col=0)
data2 = pd.read_csv("Fatostatin.csv", index_col=0)
data2["Cluster"] = list(map(lambda x : int(x.split("_")[1]), data2["Term"].values))
annotated = pd.merge(data1, data2, on='Cluster', how="inner")
annotated.to_csv("Fatostatin_annotated.csv")
