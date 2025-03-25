import os
from os import path
script_path = path.abspath(__file__)
install_dir = path.dirname(os.path.dirname(script_path))

rule bowtie2_index:
    input:
        ref=path.join(config['ref_path'], config['ref'])
    output:
        touch(path.join(config['ref_path'], config['ref'] + ".bowtie2.indexed"))
    conda:
        "../envs/findcirc.yaml"
    log:
        path.join(config['out_path'],config['out_dir'],"logs/bowtie2/index.log")
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/bowtie2_index.txt")
    shell:
        "bowtie2-build {input.ref} {input.ref} 2> {log}"

rule bowtie2_1stmap:
    input:
        fa=path.join(config['ref_path'], config['ref']),  # 基因组文件路径
        indexed=path.join(config['ref_path'], config['ref'] + ".bowtie2.indexed"),
        R1=path.join(config['out_path'], config['out_dir'], "fastp/{sample}_1.fq"),
        R2=path.join(config['out_path'], config['out_dir'], "fastp/{sample}_2.fq")
    output:
        path.join(config['out_path'], config['out_dir'], "map/bowtie2/{sample}.bam")
    conda:
        "../envs/bowtie2_samtools.yaml"
    log:
        path.join(config['out_path'], config['out_dir'], "logs/bowtie2/{sample}.log")
    threads:
        8
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_bowtie2_1stmap.txt")
    shell:
        """
        bowtie2 -p {threads} --very-sensitive --score-min=C,-15,0 --mm -x {input.fa} -1 {input.R1} -2 {input.R2} \
        | samtools view -buS - \
        | samtools sort -o {output}
        """

rule unmapped_reads:
    input:
        path.join(config['out_path'],config['out_dir'],"map/bowtie2/{sample}.bam")
    output:
        path.join(config['out_path'],config['out_dir'],"map/bowtie2/{sample}.unmapped.bam")
    conda:
        "../envs/bowtie2_samtools.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_unmapped_reads.txt")
    shell:
        "samtools view -hf 4 {input[0]} | samtools view -Sb - > {output[0]}"

rule anchors:
    input:
        path.join(config['out_path'],config['out_dir'],"map/bowtie2/{sample}.unmapped.bam")
    output:
        path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.anchors.fastq.gz")
    conda:
        "../envs/findcirc.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_anchors.txt")
    shell:
        "python2 script/unmapped2anchors.py {input[0]} | gzip > {output[0]}"

rule bowtie2_2ndmap:
    input:
        anchors_fq=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.anchors.fastq.gz"),
        fa=path.join(config['ref_path'],config['ref'])
    output:
        sites_log=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.sites.log"),
        sites_bed=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.fincirc.sites.bed"),
    log:
        path.join(config['out_path'],config['out_dir'],"logs/findcirc/{sample}.log")
    threads:
        8
    conda:
        "../envs/findcirc.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_bowtie2_2ndmap.txt")
    shell:
        """
        bowtie2 -p {threads} --reorder --mm --score-min=C,-15,0 -q -U {input.anchors_fq} -x {input.fa} | python2 script/find_circ.py -G {input.fa} -p prefix -s {output.sites_log} > {output.sites_bed}
        """

rule findcirc:
    input:
        sites_log=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.sites.log"),
        sites_bed=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.fincirc.sites.bed")
    output:
        txt=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.txt")
    conda:
        "../envs/findcirc.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_findcirc.txt")
    shell:
        '''
         grep CIRCULAR {input.sites_bed} | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | python2 script/maxlength.py 100000 > {output.txt}
        '''
