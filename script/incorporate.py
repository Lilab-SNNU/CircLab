import pandas as pd

def load_and_label(path, source, cols):
    """Load dataframe from path and label it with a source column."""
    df = pd.read_csv(path, sep="\t", header=None, usecols=cols)
    df['Source'] = source
    return df

# Corrected the snippet to use dynamic file input using snakemake's input attribute
circ_bwa = load_and_label(snakemake.input[0], 'CIRCexplorer2_bwa', [0,1,2])
circ_star = load_and_label(snakemake.input[1], 'CIRCexplorer2_star', [0,1,2])
ciri = load_and_label(snakemake.input[2], 'CIRI2', [0,1,2])
findcirc = load_and_label(snakemake.input[3], 'findcirc', [0,1,2])

# Concatenating records from all sources
result0 = pd.concat([circ_bwa, circ_star, ciri, findcirc], ignore_index=True)
result0.columns = ['chr', 'start', 'end', 'Source']

# Grouping records and joining sources
grouped = result0.groupby(['chr', 'start', 'end'])
result = grouped['Source'].agg(','.join).reset_index()

# Generating circRNA_name column
result['circRNA_name'] = result.apply(lambda x: f"{x['chr']}:{x['start']}|{x['end']}", axis=1)

# Loading detailed information for merging
info1_cols = [0, 1, 2, 5, 9, 10, 12, 13, 14, 15]
info2_cols = [0, 1, 2, 5, 9, 10, 13, 14, 15]
info_ciri_cols = [0, 1, 2, 7, 10, 11, 12, 13]
info_findcirc_cols = [0, 1, 2, 4, 5, 15]

info1 = pd.read_csv(snakemake.input[0], sep="\t", header=None, usecols=info1_cols)
info_ciri = pd.read_csv(snakemake.input[2], sep="\t", header=None, usecols=info_ciri_cols)
info_findcirc = pd.read_csv(snakemake.input[3], sep="\t", header=None, usecols=info_findcirc_cols)

info_ciri[12] = info_ciri[12].astype(str).apply(lambda x: x.split(',')[0])

# Preparing for merging detailed info - Assuming info1 has all required details and structure as the template
merged_info = pd.merge(result, info1, left_on=['chr', 'start', 'end'], right_on=[0, 1, 2], how='left').drop(columns=[0, 1, 2])

# Renaming columns for clarity
merged_info.columns = [
    'chr', 'start', 'end', 'Source', 'circRNA_name', 'strand',
    'exon_count', 'exon_size', 'junction_reads', 'circRNA_type',
    'geneName', 'isoformName'
]

# Update the merged_info with info from CIRI and findcirc - Merging strategies need revising for actual logic
info_ciri_renamed = info_ciri.rename(columns={7: 'junction_reads', 10: 'junction_ratio', 11: 'circRNA_type', 12: 'geneName', 13: 'strand'})
info_findcirc_renamed = info_findcirc.rename(columns={5: 'strand', 15: 'signal'})

# It is crucial to first update the dataframe with one set of extra data and then with the other
merged_info.update(info_ciri_renamed)
merged_info.update(info_findcirc_renamed)

# Saving the processed data to the file specified by snakemake output
merged_info.to_csv(snakemake.output[0], sep='\t', index=False)