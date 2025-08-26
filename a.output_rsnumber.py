import pandas as pd
## setting
bim_file = "merged_shared_exact.bim"
target_chr = "13"
target_pos = 49552182
output_prefix = "rs777630688"

kb_windows = [10000, 1000, 100, 50]
snp_window = 20
##
cols = ['chr', 'rsid', 'cm', 'pos', 'a1', 'a2']
bim_df = pd.read_csv(bim_file, sep='\t', header=None, names=cols, dtype={'chr': str})
chr_df = bim_df[bim_df['chr'] == str(target_chr)].reset_index(drop=True)

target_index = (chr_df['pos'] - target_pos).abs().idxmin()

start_idx = max(0, target_index - snp_window)
end_idx = min(len(chr_df), target_index + snp_window + 1)
snp_df = chr_df.loc[start_idx:end_idx]
snp_df.to_csv(f"{output_prefix}_{snp_window}snp_full.txt", sep='\t', index=False)
snp_df['rsid'].to_csv(f"{output_prefix}_{snp_window}snp_list.txt", index=False, header=False)

for kb in kb_windows:
    start_pos = target_pos - kb * 1000
    end_pos = target_pos + kb * 1000
    kb_df = chr_df[(chr_df['pos'] >= start_pos) & (chr_df['pos'] <= end_pos)]
    kb_df.to_csv(f"{output_prefix}_{kb}kb_full.txt", sep='\t', index=False)
    kb_df['rsid'].to_csv(f"{output_prefix}_{kb}kb_list.txt", index=False, header=False)

print("Step1 DONE!")

