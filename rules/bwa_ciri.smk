import os
from os import path

rule bwa_index:
    input:
        ref=path.join(config['ref_path'],config['ref'])
    output:
        touch(path.join(config['ref_path'], config['ref'] + ".bwa.indexed"))
    conda:
        "../envs/CIRI2.yaml"
    shell:
        "bwa index {input[0]} {output[0]}"

rule bwa_mem:
    input:
        fa=path.join(config['ref_path'],config['ref']),# 基因组文件路径
        indexed=path.join(config['ref_path'],config['ref'] + ".bwa.indexed"),
        R1=path.join(config['out_path'],config['out_dir'],"fastp/{sample}_1.fq"),
        R2=path.join(config['out_path'],config['out_dir'],"fastp/{sample}_2.fq")
    output:
        path.join(config['out_path'],config['out_dir'],"map/bwa/{sample}.bam")
    conda:
        "../envs/CIRI2.yaml"
    log:
        path.join(config['out_path'],config['out_dir'],"logs/bwa/{sample}.log")
    threads:
        8
    shell:
        "bwa mem -T 19 -t {threads} {input.fa} {input[2]} {input[3]} -o {output} 2> {log}"

rule CIRI2:
    input:
        rules.bwa_mem.output,
        gtf=config['gtf'],
        fa=path.join(config['ref_path'],config['ref'])
    output:
        path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}")
    log:
        path.join(config['out_path'],config['out_dir'],"logs/CIRI2/{sample}.log")
    threads:
        4
    conda:
        "../envs/CIRI2.yaml"
    shell:
        "perl {config[tools_path]}/script/CIRI2.pl -F {input[2]} -I {input[0]} -O {output[0]} -A {input[1]} -T {threads} 2> {log}"

