import pandas as pd

fpath=snakemake.input[0]
anno=pd.read_csv(fpath,sep='\t',header=None)
#print(anno.head())
gene = anno[11]
anno.insert(0,'gene',gene)
anno.to_csv(snakemake.output[0],sep='\t',index=False,header=None)
