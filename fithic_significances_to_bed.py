import pandas as pd

input_file = 'FXR/KO_GW_significances.txt'
output_file = 'FXR/KO_GW_significant.bed'
read_length = 50000

df = pd.read_csv(input_file, sep='\t')
print('num interactions', len(df))

# filter by p-value
df_filtered = df[df['p-value'] < 0.05]

print('num significant interactions:', len(df_filtered))

# sort by p-value and take top 5000
df_top = df_filtered.nsmallest(5000, 'p-value')

def midpoint_to_bed(chrom, mid, length=50000):
    start = int(mid) - length // 2
    end = int(mid) + length // 2
    return chrom, max(0, start), end  # make sure start is non-negative

bed_entries = []
for _, row in df_top.iterrows():
    bed_entries.append(midpoint_to_bed(row['chr1'], row['fragmentMid1']))
    bed_entries.append(midpoint_to_bed(row['chr2'], row['fragmentMid2']))

with open(output_file, 'w') as f:
    for chrom, start, end in bed_entries:
        f.write(f"{chrom}\t{start}\t{end}\n")

print(f"BED file written to {output_file} with {len(bed_entries)} entries.")
