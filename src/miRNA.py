import sys
import os
import argparse
import re
import time
import subprocess
from os import path

def circRNA_mi(miRNA_file, circRNA_file, circRNA_miRNA):
    GSTAr_cmd=f'perl script/GSTAr.pl {miRNA_file} {circRNA_file} > {circRNA_miRNA}'
    print(GSTAr_cmd)
    return_code=os.system(GSTAr_cmd) >> 8
    if return_code:
        sys.exit('Error: cannot work! Please check your file')
    return circRNA_miRNA

def m_miRNA(miRNA_file, mRNA, m_miRNA):
    GSTAr_cmd=f'perl script/GSTAr.pl {miRNA_file} {mRNA} > {m_miRNA}'
    print(GSTAr_cmd)
    return_code = os.system(GSTAr_cmd) >> 8
    if return_code:
        sys.exit('Error: cannot work! Please check your file')
    return m_miRNA

def circRNA_decoy(miRNA_file, circRNA_file, output_dir):
    '''

    '''
    print(get_time(), f"Predict circRNA-miRNA-decoy...")
    print(f"-miRNA file: {miRNA_file}") # unique miRNA
    print(f"-circRNA file: {circRNA_file}") # circRNA fasta file

    result_path = path.join(output_dir, "circRNA-mi-decoy.txt") # -o 参数 不包含文件名
    print(f"-out_file: {result_path}") # out_dir+file

    circRNA_miRNA=path.join(f'{output_dir}',"circRNA-miRNA.txt") #circRNA_mi中的形参
    circ_miRNA=circRNA_mi(miRNA_file, circRNA_file, circRNA_miRNA)
    #上一步定义的形参传入函数circRNA_mi中，返回circRNA-mi-decoy.txt保存在变量circ_miRNA中

    new = False
    result = ''
    result0 = ''
    with open(circ_miRNA, 'rt') as f_in,open(result_path, 'wt') as f_out:
        for line in f_in:
            if re.match(r'^5', line):  # re.match()匹配开头
                gene_inf = line.rstrip('\n')
                gene_infs = gene_inf.split(' ')
                gene_seq = gene_infs[1]
                new = True
            if new and re.match(r'^ ', line):
                align0 = line[3:].rstrip('\n')
            if re.match(r'^3', line):
                new = False
                miRNA_inf = line.rstrip('\n')
                miRNA_infs = miRNA_inf.split(' ')
                miRNA_seq = miRNA_infs[1][::-1]  # 提取miRNA序列并且反转
                nucleotides = list(miRNA_seq)

            mfeperfect_match = re.search(r'MFE of perfect match: (\S+)', line)
            if mfeperfect_match:
                mfeperfect = mfeperfect_match.group(1)

            mfesite_match = re.search(r'MFE of this site: (\S+)', line)
            if mfesite_match:
                mfesite = mfesite_match.group(1)

            mferatio_match = re.search(r'MFEratio: (\S+)', line)
            if mferatio_match and 'Sorted by' not in line:
                mferatio = mferatio_match.group(1)
                condition1 = condition2 = condition3 = 0
                '''
                condition1&2
                '''
                start_end9_12 = []
                start_end2_8_0 = []
                m = 0  # m为有效碱基，n为有效碱基的位置
                for n, char in enumerate(nucleotides):
                    if char in ['A', 'U', 'C', 'G']:
                        m += 1
                    if char in ['A', 'U', 'C', 'G'] and (m == 2 or m == 8):
                        start_end2_8_0.append(n)  # 有效碱基的位置
                        if len(start_end2_8_0) == 2:
                            start_end2_8 = start_end2_8_0
                    if char in ['A', 'U', 'C', 'G'] and (m == 9 or m == 12):
                        start_end9_12.append(n)
                        if len(start_end9_12) == 2:
                            break
                m = 0

                align = align0[::-1]
                aligns = list(align)
                indel9_12 = 0
                for i in range(start_end9_12[0], start_end9_12[1] + 1):
                    if aligns[i] == ' ':
                        indel9_12 += 1

                match2_8 = 0
                for i in range(start_end2_8[0], start_end2_8[1] + 1):
                    if aligns[i] == '|':
                        match2_8 += 1

                condition1 = 1 if indel9_12 > 1 and indel9_12 < 6 else 0
                condition2 = 1 if match2_8 == 7 and (start_end2_8[0] == 1 and start_end2_8[1] == 7) else 0
                '''
                condition3
                '''
                indel_all = aligns.count(' ')
                condition3 = 1 if indel_all - indel9_12 <= 4 else 0

                if condition1 and condition2 and condition3:
                    miRNA_name = miRNA_infs[4]
                    gene_name_position = gene_infs[4]
                    gene_name = ':'.join(gene_name_position.split(":")[-3:-1])  # 左开右闭
                    position = gene_name_position.split(":")[1]
                    miRNA_decoy = f'miRNA_target: 5 {gene_infs[1]} 3'
                    miRNA = f'miRNA:        3 {miRNA_infs[1]} 5'
                    q = ' ' * 16
                    result += f'{miRNA_name}\t{gene_name}\t{position}\t{mfeperfect}\t{mfesite}\t{mferatio}\n{miRNA_decoy}\n{q}{align0}\n{miRNA}\n'
                    result0 += f'{miRNA_name}\t{gene_name}\t{position}\t{mfeperfect}\t{mfesite}\t{mferatio}\t{miRNA_seq}\t{gene_seq}\n'
                indel9_12 = match2_8 = indel_all = condition1 = condition2 = condition3 = 0
                start_end9_12 = start_end2_8 = []
        result0 = f'miRNA_name\tTranscript\tstart_end\tMFEperfect\tMFEsite\tMFEratio\n{result0}'
        f_out.write(result0 + '\n')

def miRNA_target_prediction(miRNA_file, mRNA, output_dir):
    '''
    select appropriate sequence as miRNA target
    '''
    print(get_time(), f"Predict miRNA target...")
    print(f"-miRNA file: {miRNA_file}")
    print(f"-mRNA file: {mRNA}")
    result_path = path.join(output_dir, 'm-mi-target.txt')
    print(f"-out_file: {result_path}")

    mRNA_miRNA_decoy = path.join(output_dir, "mi-mRNA.txt")
    m_miRNA_result = m_miRNA(miRNA_file, mRNA, mRNA_miRNA_decoy)

    new = False
    m = 0
    nucleotides = []
    aligns = []
    result = ''
    result0 = '' # different type
    indel9_12 =  condition1 = condition2 = condition3 = 0
    with open(m_miRNA_result,'r') as f_in,open(result_path,'wt') as fo:
        for line in f_in:
            if line.startswith('5'):
                gene_inf = line.rstrip('\n')
                gene_infs = gene_inf.split(' ')
                gene_seq = gene_infs[1]
                new = True
            if new and line.startswith(' '):
                align0 = line[3:].rstrip('\n')
                align = align0[::-1]
                aligns = list(align)
            if line.startswith('3'):
                new = False
                miRNA_inf = line.rstrip('\n')
                miRNA_infs = miRNA_inf.split(' ')
                miRNA_seq = miRNA_infs[1][::-1]
                nucleotides = list(miRNA_seq)
            if 'MFE of perfect match' in line:
                mfeperfect = line.split(' ')[4].rstrip('\n')
            if 'MFE of this site' in line:
                mfesite = line.split(' ')[4].rstrip('\n')
            if 'MFEratio' in line and 'Sorted by' not in line:
                mferatio = line.split(' ')[1].rstrip('\n')
                '''
                condition1
                '''
                start_end9_12 = []
                for n, char in enumerate(nucleotides):
                    if char in ['A', 'U', 'C', 'G']:
                        m += 1
                    if char in ['A', 'U', 'C', 'G'] and (m == 9 or m == 12):
                        start_end9_12.append(n)
                        if len(start_end9_12) == 2:
                            break
                m = 0
                indel9_12_l = []
                indel9_12 = 1
                if len(start_end9_12) >= 2:
                    for i in range(start_end9_12[0],start_end9_12[1]+1):
                        if aligns[i] == ' ':
                            indel9_12 += 1

                    indel9_12_l.append(indel9_12)
                condition1 = 1 if indel9_12 <= 1 else 0
                '''
                condition2
                '''
                mismatch = 0
                for n in range(len(aligns)-1):
                    bi = aligns[n] + aligns[n+1]
                    if '  ' in bi:
                        mismatch += 1
                condition2 = 1 if not mismatch else 0
                '''
                condition3
                '''
                indel_all = aligns.count(' ')
                condition3 = 1 if indel_all - indel9_12 <= 4 else 0
                if condition1 and condition2 and condition3:
                    # f_out.write('1')
                    miRNA_name = miRNA_infs[4]
                    gene_name_position = gene_infs[4]
                    gene_name = ':'.join(gene_name_position.split(":")[-3:-1])
                    position = gene_name_position.split(":")[1]  # target in mRNA
                    # position = gene_name_position.split(":")[2] # miRNA_target in circ
                    miRNA_target = f'miRNA_target: 5 {gene_infs[1]} 3'
                    miRNA = f'miRNA:        3 {miRNA_infs[1]} 5'
                    q = ' '*16
                    result += f'{miRNA_name}\t{gene_name}\t{position}\t{mfeperfect}\t{mfesite}\t{mferatio}\n{miRNA_target}\n{q}{align0}\n{miRNA}\n'
                    result0 += f'{miRNA_name}\t{gene_name}\t{position}\t{mfeperfect}\t{mfesite}\t{mferatio}\t{miRNA_seq}\t{gene_seq}\n'
                indel9_12 = indel_all = condition1 = condition2 = condition3 = 0

        result0 = f'miRNA_name\tTranscript\tstart_end\tMFEperfect\tMFEsite\tMFEratio\n{result0}'
        result = f'miRNA_name\tTranscript\tstart_end\tMFEperfect\tMFEsite\tMFEratio\n{result}'
        # f_out.write(result + '\n')
        fo.write(result0 + '\n')
        fo.close()

def get_time():
    nowtime = time.strftime('[%Y-%m-%d %H:%M:%S]', time.localtime())
    return nowtime

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='circRNA function detection module',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Add command line arguments
    parser.add_argument('-mi', '--miRNA', help='Input file(s) for miRNA sequence, FASTA/FASTQ file', required=True)
    parser.add_argument('-m', '--mRNA', help='Input file for mRNA sequence, FASTA/FASTQ file', required=True)
    parser.add_argument('-i', '--circRNA', dest='input_file', help='Input file(s) for circRNA sequence, FASTA/FASTQ file', required=True)
    parser.add_argument('-o', '--output', dest='output_dir', help='Output directory', required=True)

    args = parser.parse_args(argv)

    circRNA_decoy(args.miRNA, args.input_file, args.output_dir)
    print(get_time(), f"circRNA-miRNA-decoy predict finished!!,Please check your file in:{args.output_dir}")

    miRNA_target_prediction(args.miRNA, args.mRNA, args.output_dir)
    print(get_time(), f"Predict miRNA target finished!!,Please check your file in:{args.output_dir}")

if __name__ == "__main__":
    # detection()
    main(sys.argv[1:])