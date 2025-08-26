# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import gzip

# === Setting ===
prefix = 'prefix_of_the_array_data' # prefix of the array data
rs_number = 'rs761167763'  # rsnumber of the variant
range = '_1000kb' # range or snp
vcf_path = prefix+'_'+rs_number+range+'_5.phased.vcf.gz'  # phased VCF

fam_files = f"{prefix}_{rs_number}_phenotyped.fam"

fam_df = pd.read_csv(fam_files, sep=r'\s+', header=None,
                     names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])
pheno_map = dict(zip(fam_df['IID'], fam_df['PHENO']))

# === All SNP in VCF ===
target_rows = []
with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            sample_names = header[9:]
        elif not line.startswith('#'):
            row = line.strip().split('\t')
            target_rows.append(row)

vcf_df = pd.DataFrame(target_rows, columns=header[:9] + sample_names)

# === generate 2 haplotype combinations ===
hap_data = {s: {'H1': [], 'H2': []} for s in sample_names}
for _, row in vcf_df.iterrows():
    for s in sample_names:
        gt = row[s].split(':')[0]
        h1, h2 = gt.split('|') if '|' in gt else ('.', '.')
        hap_data[s]['H1'].append(h1)
        hap_data[s]['H2'].append(h2)

# ===  haplotype dataframe ===
haplo_df = pd.DataFrame({
    'Sample': sample_names,
    'Haplotype_1': ['-'.join(hap_data[s]['H1']) for s in sample_names],
    'Haplotype_2': ['-'.join(hap_data[s]['H2']) for s in sample_names],
    'Cleaned_Sample': [
        '_'.join(s.split('_')[:4]) if s.startswith('TV') else '_'.join(s.split('_')[:2])
        for s in sample_names
    ]
})
haplo_df['Phenotype'] = haplo_df['Cleaned_Sample'].map(pheno_map)

# === count haplotype - case/control number ===
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

print(f" Total case samples: {total_case}")
print(f" Total control samples: {total_control}")
if total_case == 0 or total_control == 0:
    print("No enough case/control data, skipping analysis.")
    exit()

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

# === Ready for permutation test ===
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

# === Output ===
df_result = df_result.sort_values('Fisher_P')
df_result.to_csv(rs_number + range + '_haplotype_results.csv', index=False)
print(f"DONE! {rs_number}{range}_haplotype_results.csv")
