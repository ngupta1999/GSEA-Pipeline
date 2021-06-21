import pandas as pd
import scipy.stats
import numpy as np

achilles = pd.read_csv("Achilles_proximity_corrected.csv", index_col = 0)
modules= pd.read_csv("new_modules_d_0.9.csv", index_col = 0)

rows = []

counter = 0 

for index, row in modules.iterrows():
    counter += 1
    print(counter)
    if counter % 100 == 0:
        print("counter is ", counter)
    
    members = row["Members"]
    genes = members.split(" ")
     
    for g1 in genes:
       
        num_positives = 0 
        num_negatives = 0 

        total_positive = 0 
        total_negative = 0
        for g2 in genes:
            if g1 != g2:
                z = scipy.stats.pearsonr(achilles.loc[g1].values, achilles.loc[g2].values)[0]
                if z >= 0 : 
                    num_positives += 1
                    total_positive += z
                else:
                    num_negatives += 1
                    total_negative += z
                    
                    
        
        avg_positive = 0
        avg_negative = 0
        if num_negatives > 0:
            avg_negative = total_negative / num_negatives
            
        if num_positives > 0:
            avg_positive = total_positive / num_positives
        
        r = [index, g1, num_positives, num_negatives, avg_positive, avg_negative]
        rows.append(r)
headers = ["module", "gene", "num positive", "num negative", "avg positive", "avg negative" ]
df = pd.DataFrame(rows, columns = headers)
df.to_csv("nile")
