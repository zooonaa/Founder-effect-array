import pandas as pd
import os

# 輸入參數
rs_number = 'rs777630688'  # 用來標記輸出檔案
case_file = 'rs777630688_case_list.txt'  # 你的 case ID 清單檔案
fam_files = [
    'merged_nofamily_nonosex.qc.noindel.fam',  #一個rs number 產出一個phenotyped就可以了
#    'merged_all.qc.noindel_'+rs_number+'_100kb.fam'
]

# 讀取 case list
with open(case_file, 'r') as f:
    case_ids = set(line.strip() for line in f if line.strip())

# 處理每一個 fam 檔案
for fam_path in fam_files:
    fam = pd.read_csv(fam_path, sep=' ', header=None)
    fam.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE']
    
    # 從 IID 擷取 sample prefix
    fam['SAMPLE_ID'] = fam['IID'].apply(lambda x: x.split('_')[0])
    
    # 標註 phenotype：case = 2, control = 1
    fam['PHENOTYPE'] = fam['SAMPLE_ID'].apply(lambda x: 2 if x in case_ids else 1)
    fam.drop(columns=['SAMPLE_ID'], inplace=True)
    
    # 輸出檔案名稱
    base = os.path.basename(fam_path).replace('.fam', '')
    output_name = f"merged_nofamily_nonosex.qc.noindel_"+rs_number+"_phenotyped.fam"
    fam.to_csv(output_name, sep=' ', index=False, header=False)
    
    print(f"✅ Saved: {output_name}")
