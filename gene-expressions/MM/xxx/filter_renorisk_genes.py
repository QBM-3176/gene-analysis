import pandas as pd

counts_df = pd.read_csv("GSE226353_ex_count.csv", index_col=0, sep=';')
filter_df = counts_df[counts_df['Geneid'].isin(['ENSMUSG00000025075', 'ENSMUSG00000058952', 'ENSMUSG00000026365', 'ENSMUSG00000058806', 'ENSMUSG00000026405'])]

# Pivot: set Geneid as index and transpose so samples are rows and genes are columns
# filter_df = filter_df.set_index('Geneid').T

filter_df.to_csv("GSE226353_renorisk_gene_counts.csv")