import pandas as pd
import os

# === Setting ===
rs_number = 'rs777630688'  # rsnumber of the variant
case_file = 'rs777630688_case_list.txt'  # The 'case' list

prefix = 'prefix_of_the_array_data' # prefix of the array data

fam_files = [f"{prefix}.fam"]

with open(case_file, 'r') as f:
    case_ids = set(line.strip() for line in f if line.strip())

for fam_path in fam_files:
    fam = pd.read_csv(fam_path, sep=' ', header=None)
    fam.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE']
    
    fam['SAMPLE_ID'] = fam['IID'].apply(lambda x: x.split('_')[0])
    
    fam['PHENOTYPE'] = fam['SAMPLE_ID'].apply(lambda x: 2 if x in case_ids else 1)
    fam.drop(columns=['SAMPLE_ID'], inplace=True)
    
    output_name = f"{prefix}_{rs_number}_phenotyped.fam"
    fam.to_csv(output_name, sep=' ', index=False, header=False)
    
    print(f"Saved: {output_name}")
