import os
from os import path

script_path = path.abspath(__file__)
install_dir = path.dirname(os.path.dirname(script_path))
config_path = path.join(install_dir, 'config', 'config.yaml')
configfile: config_path

module bwa:
    snakefile: "bwa_ciri.smk"
    config: config

module gtf2anno:
    snakefile: "gtf2anno.smk"
    config: config

use rule bwa_mem from bwa as bwa
use rule bwa_index from bwa as bwa_index

use rule gtf2ano1 from gtf2anno as gtf2ano1
use rule gtf2ano2 from gtf2anno as gtf2ano2

rule circ_parse_bwa:
    input:
        rules.bwa.output
    output:
        bsj=path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_bsj.bed")
    log:
        parse=path.join(config['out_path'],config['out_dir'],"logs/CIRCexplorer2/bwa/{sample}_parse.log")
    conda:
        "../envs/CIRC2.yaml"
    threads:
        6
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_circ_parse_bwa.txt")
    shell:
        "CIRCexplorer2 parse -t BWA {input[0]} -b {output[bsj]} > {log[parse]}"

rule circ_anno_bwa:
    input:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/annotation.txt"),
        fa=path.join(config['ref_path'],config['ref']),
        bsj=path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_bsj.bed")
    output:
        txt=path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}")
    log:
        anno = path.join(config['out_path'],config['out_dir'],"logs/CIRCexplorer2/bwa/{sample}_annotate.log")
    conda:
        "../envs/CIRC2.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_circ_anno_bwa.txt")
    shell:
        "CIRCexplorer2 annotate -r {input[0]} -g {input[1]} -b {input[bsj]} -o {output[txt]} > {log[anno]}"