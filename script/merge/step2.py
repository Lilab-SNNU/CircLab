import pandas as pd

circ = pd.read_csv(snakemake.input[0], sep='\t', header=None)
circ[1]= circ[1]+1
circ[3] = circ[3].str.replace(' ', '')
circ[3] = circ[5].apply(lambda x: 'forward' if x == '+' else 'reverse')
circ_i = circ[[0,1,2,3,5,6]]
circ_i.to_csv(snakemake.output[0], sep='\t', header=False, index=False)

circ2 = pd.read_csv(snakemake.input[1], sep='\t', header=None)
circ2[1]= circ2[1]+1
circ2[3] = circ2[3].str.replace(' ', '')
circ2[3] = circ2[5].apply(lambda x: 'forward' if x == '+' else 'reverse')
circ2_i = circ2[[0,1,2,3,5,6]]
circ2_i.to_csv(snakemake.output[1], sep='\t', header=False, index=False)

ciri = pd.read_csv(snakemake.input[2], sep='\t', header=None)
ciri[3] = ciri[3].str.replace(' ', '')
ciri[3] = ciri[5].apply(lambda x: 'forward' if x == '+' else 'reverse')
ciri_i = ciri[[0,1,2,3,5,6]]
ciri_i.to_csv(snakemake.output[2], sep='\t', header=False, index=False)

findcirc = pd.read_csv(snakemake.input[3], sep='\t', header=None)
findcirc[1]= findcirc[1]+1
findcirc[3] = findcirc[3].str.replace(' ', '')
findcirc[3] = findcirc[5].apply(lambda x: 'forward' if x == '+' else 'reverse')
findcirc_i = findcirc[[0,1,2,3,5,6]]
findcirc_i.to_csv(snakemake.output[3], sep='\t', header=False, index=False)