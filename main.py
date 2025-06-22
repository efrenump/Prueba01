import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os
import time

from Bio import Entrez
from Bio import SeqIO

LOCAL_DATA_FILE = 'human_proteins_refseq.csv'

Entrez.email = "efrenump@gmail.com"
Entrez.api_key = "cfdb87856942d4a8dafe892ee2613c24c208"

if os.path.exists(LOCAL_DATA_FILE):
    print(f"Local data file found. Loading DataFrame from '{LOCAL_DATA_FILE}'...")
    df = pd.read_pickle(LOCAL_DATA_FILE)
    print("DataFrame loaded successfully.")

else:
    print("Searching for all human proteins...")
    handle = Entrez.esearch(db="protein", term="txid9606[Organism] AND refseq[filter]", usehistory="y")
    search_results = Entrez.read(handle)
    handle.close()
    
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    
    print(f"Found {count} human protein records.")
    
    batch_size = 500
    all_protein_records = []
    
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print(f"Downloading records for {start+1} to {end}")
        
        try:
         fetch_handle = Entrez.efetch(db="protein",
                                      rettype="gp",
                                      retmode="text",
                                      retstart=start,
                                      retmax=batch_size,
                                      webenv=webenv,
                                      query_key=query_key)
    
         for record in SeqIO.parse(fetch_handle, "genbank"):
             all_protein_records.append(record)
         
         fetch_handle.close()
         
        except Exception as e:
            print(f"An error occurred during fetch or parse: {e}")
            print("Continuing with the next batch...")
    
        continue
    
        time.sleep(0.11)
        
    print(f"\nFinished fetching. Total records in data structure: {len(all_protein_records)}")
    
    print("\nCreating Pandas DataFrame...")
    try:
        protein_data = {
            "id": [rec.id for rec in all_protein_records],
            "name": [rec.name for rec in all_protein_records],
            "description": [rec.description for rec in all_protein_records],
            "sequence": [str(rec.seq) for rec in all_protein_records],
            "length": [len(rec.seq) for rec in all_protein_records]
        }
        df = pd.DataFrame(protein_data)
        print("Successfully created a Pandas DataFrame.")
    
        print(f"Saving DataFrame to '{LOCAL_DATA_FILE}' for future use...")
        df.to_pickle(LOCAL_DATA_FILE)
        print("Save complete.")
    
    except Exception as e:
        print(f"Could not create Pandas DataFrame. Error: {e}")
        
normal_df_length = np.log(df["length"])

mean = normal_df_length.mean()
std = normal_df_length.std(ddof=0)
x = np.linspace(0, 10, 1000)
y = np.exp(-np.pow((x - mean)/std,2)/2)/(std * np.sqrt(2 * np.pi))

sns.set_theme(rc={'figure.figsize':(140,6)})
sns.kdeplot(normal_df_length, gridsize=2000, color="red", label='Dataset')
plt.plot(x,y, label='Prediction')
plt.legend()

normal_df_length.describe()





 