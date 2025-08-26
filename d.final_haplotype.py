import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

# === 設定輸入檔案 ===
rs_number = 'rs116802390'  # 用來標記輸出檔案
output_prefix = 'rs116802390'  # prefix for result output
condition = '_20snp'
vcf_path = rs_number+condition+'_5.phased.vcf.gz'  # phased VCF
fam_path = 'merged_all.qc.noindel_'+rs_number+'_phenotyped.fam'  # .fam with phenotype

# === 讀取 .fam 檔案，建立 Phenotype 對應表 ===
fam_df = pd.read_csv(fam_path, sep='\s+', header=None,
                     names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])
pheno_map = dict(zip(fam_df['IID'], fam_df['PHENO']))

# === 擷取 VCF 所有 SNP ===

target_rows = []
import gzip
with gzip.open(vcf_path, 'rt') as f:  # 'rt' 表示 text 模式解壓讀取
    for line in f:
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            sample_names = header[9:]
        elif not line.startswith('#'):
            row = line.strip().split('\t')
            target_rows.append(row)

vcf_df = pd.DataFrame(target_rows, columns=header[:9] + sample_names)

# === 建立每個樣本的兩條 haplotype 組合 ===
hap_data = {s: {'H1': [], 'H2': []} for s in sample_names}
for _, row in vcf_df.iterrows():
    for s in sample_names:
        gt = row[s].split(':')[0]
        h1, h2 = gt.split('|') if '|' in gt else ('.', '.')
        hap_data[s]['H1'].append(h1)
        hap_data[s]['H2'].append(h2)

# === 整理成 haplotype dataframe ===
haplo_df = pd.DataFrame({
    'Sample': sample_names,
    'Haplotype_1': ['-'.join(hap_data[s]['H1']) for s in sample_names],
    'Haplotype_2': ['-'.join(hap_data[s]['H2']) for s in sample_names],
    'Cleaned_Sample': [s.split('_')[0] + '_' + s.split('_')[1] for s in sample_names]
})
haplo_df['Phenotype'] = haplo_df['Cleaned_Sample'].map(pheno_map)

# === 統計每個 haplotype 的 case/control 出現次數 ===
hap_count = {}
for _, row in haplo_df.iterrows():
    pheno = row['Phenotype']
    if pheno not in [1, 2]: continue
    for h in [row['Haplotype_1'], row['Haplotype_2']]:
        if h not in hap_count:
            hap_count[h] = {'case': 0, 'control': 0}
        hap_count[h]['case' if pheno == 2 else 'control'] += 1

total_case = sum([v['case'] for v in hap_count.values()])
total_control = sum([v['control'] for v in hap_count.values()])

# === Fisher's exact test ===
results = []
for hap, count in hap_count.items():
    case = count['case']
    control = count['control']
    table = [[case, control], [total_case - case, total_control - control]]
    _, p = fisher_exact(table)
    results.append({
        'Haplotype': hap,
        'Case_Count': case,
        'Case_Percent': case / total_case * 100,
        'Control_Count': control,
        'Control_Percent': control / total_control * 100,
        'Fisher_P': p
    })
df_result = pd.DataFrame(results)

# === 準備 permutation test ===
hap_expanded = []
for _, row in haplo_df.iterrows():
    for h in [row['Haplotype_1'], row['Haplotype_2']]:
        hap_expanded.append({'Haplotype': h, 'Phenotype': row['Phenotype']})
hap_exp = pd.DataFrame(hap_expanded).dropna()

hap_list = df_result['Haplotype'].tolist()
obs_cases = {h: hap_count[h]['case'] for h in hap_list}
perm_counts = {h: 0 for h in hap_list}

# === permutation test ===
np.random.seed(42)
true_pheno = hap_exp['Phenotype'].values
hap_seq = hap_exp['Haplotype'].values
for _ in range(100000):
    shuffled = np.random.permutation(true_pheno)
    temp = {h: 0 for h in hap_list}
    for h, p in zip(hap_seq, shuffled):
        if h in temp and p == 2:
            temp[h] += 1
    for h in hap_list:
        if temp[h] >= obs_cases[h]:
            perm_counts[h] += 1

perm_pvals = {h: (perm_counts[h] + 1) / 100001 for h in hap_list}
df_result['Permutation_P'] = df_result['Haplotype'].map(perm_pvals)

# === 輸出結果 ===
df_result = df_result.sort_values('Fisher_P')
df_result.to_csv(output_prefix + condition + '_haplotype_results.csv', index=False)
print(f"✅ 結果已儲存至 {output_prefix}_haplotype_results.csv")
