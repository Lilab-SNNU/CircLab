import pandas as pd

def deal_circ(circout, circ_name1, recirc1):
    circ = pd.read_csv(circout,sep='\t',header=None)
    name = pd.read_csv(circ_name1,sep='\t',header=None)
    circ[1] = circ[1] + 1
    merged_df = pd.merge(name, circ, left_on=[0, 1, 2], right_on=[0, 1, 2])
    #print(merged_df)
    merged_df.to_csv(recirc1, sep='\t', header=None, index=None)

def deal_ciri(ciriout, circ_name2, recirc2):
    ciri = pd.read_csv(ciriout, sep='\t')
    name = pd.read_csv(circ_name2, sep='\t', header=None)
    name[1] = name[1].astype(int)
    ciri['circRNA_start'] = ciri['circRNA_start'].astype(int)
    # Merge the DataFrames
    merged_df = pd.merge(name, ciri, left_on=[0, 1, 2], right_on=['chr', 'circRNA_start', 'circRNA_end'])
    #result = merged_df.iloc[:, 6:]
    merged_df.to_csv(recirc2, sep='\t', header=None, index=None)

def deal_findcirc(findcircout, circ_name3, recirc3):
    findcirc = pd.read_csv(findcircout,sep='\t',header=None)
    name = pd.read_csv(circ_name3,sep='\t',header=None)
    findcirc[1] = findcirc[1] + 1
    df_find = findcirc.loc[findcirc[4] >= 2, :]
    merged_df = pd.merge(name, df_find, left_on=[0, 1, 2], right_on=[0, 1, 2])
    merged_df.to_csv(recirc3, sep='\t', header=None, index=None)

# circ_name1s = snakemake.input[0:1]
# circouts = snakemake.input[4:5]
# recirc1s = snakemake.output[0:1]
#
# for circout, circ_name1, recirc1 in zip(circouts, circ_name1s, recirc1s):
#     print(circout,circ_name1,recirc1)
#     deal_circ(circout, circ_name1, recirc1)

def main():
    # 直接从 snakemake 命名输入中获取文件路径
    # ----------------------------------------------
    # 1. 处理 CIRCexplorer2 (BWA)
    deal_circ(
        circout=snakemake.input.bwa_raw,
        circ_name1=snakemake.input.bwa_filtered,
        recirc1=snakemake.output[0]
    )

    # 2. 处理 CIRCexplorer2 (STAR)
    deal_circ(
        circout=snakemake.input.star_raw,
        circ_name1=snakemake.input.star_filtered,
        recirc1=snakemake.output[1]
    )

    # 3. 处理 CIRI2
    deal_ciri(
        ciriout=snakemake.input.ciri2_raw,
        circ_name2=snakemake.input.ciri2_filtered,
        recirc2=snakemake.output[2]
    )

    # 4. 处理 FindCirc
    deal_findcirc(
        findcircout=snakemake.input.findcirc_raw,
        circ_name3=snakemake.input.findcirc_filtered,
        recirc3=snakemake.output[3]
    )


if __name__ == "__main__":
    main()