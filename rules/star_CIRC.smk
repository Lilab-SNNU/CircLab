import os
from os import path
script_path = path.abspath(__file__)
install_dir = path.dirname(os.path.dirname(script_path))
config_path = path.join(install_dir, 'config', 'config.yaml')
configfile: config_path

module gtf2anno:
    snakefile: "gtf2anno.smk"
    config: config

use rule gtf2ano1 from gtf2anno as gtf2ano1
use rule gtf2ano2 from gtf2anno as gtf2ano2

rule star_index:
    input:
        ref=path.join(config['ref_path'],config['ref']),
        gtf=config['gtf']
    output:
        directory(path.join(config['ref_path'],"STAR_index"))
    conda:
        "../envs/CIRC2.yaml"
    threads:
        8
    log:
        path.join(config['out_path'],config['out_dir'],"logs/star/index.log")
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/star_index.txt")
    shell:
        '''
        STAR --runThreadN {threads} --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.ref} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 150 \
        --genomeSAindexNbases 12 > {log}
        '''

rule star:
    input:
        path.join(config['ref_path'],"STAR_index"),
        R1 = path.join(config['out_path'],config['out_dir'],"fastp/{sample}_1.fq"),
        R2 = path.join(config['out_path'],config['out_dir'],"fastp/{sample}_2.fq")
    output:
        outdir=path.join(config['out_path'],config['out_dir'],"map/star/{sample}"),
        bam=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Aligned.out.bam"),
        # log_final=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Log.final.out"),
        # log=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Log.out"),
        # log_progress=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Log.progress.out"),
        # sj_out=path.join(config['out_path'],config['out_dir'],"map/star/{sample}SJ.out.tab"),
        junction=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Chimeric.out.junction")
        # ReadsPerGene=path.join(config['out_path'],config['out_dir'],"map/star/{sample}ReadsPerGene.out.tab")
    #Ara_sampleLog.progress.out,Ara_sampleLog.out,Ara_sampleChimeric.out.junction,Ara_sampleLog.final.out,
    #Ara_sampleReadsPerGene.out.tab,Ara_sampleAligned.out.bam,Ara_sampleSJ.out.tab
    conda:
        "../envs/CIRC2.yaml"
    log:
        path.join(config['out_path'],config['out_dir'],"logs/star/{sample}.log")
    threads:
        8
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_star.txt")
    shell:
        "touch {output[0]} "
        '''
        STAR \
        --genomeDir {input[0]} \
        --runThreadN {threads} \
        --readFilesIn {input.R1} {input.R2} \
        --outFileNamePrefix {output.outdir} \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --outSAMtype BAM Unsorted \
        --outSAMattributes All \
        --quantMode GeneCounts > {log}
        '''
#  --readFilesCommand zcat \
#  --outFilterIntronMotifs RemoveNoncanonical > {log} 2>&1
#bam=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Aligned.out.bam"),
#log_final=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Log.final.out"),
#log=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Log.out"),
#log_progress=path.join(config['out_path'],config['out_dir'],"map/star/{sample}Log.progress.out"),
#sj_out=path.join(config['out_path'],config['out_dir'],"map/star/{sample}SJ.out.tab")
rule star_bam_index:
    input:
        bam=path.join(config['out_path'], config['out_dir'], "map/star/{sample}Aligned.out.bam")
    output:
        bam_index=path.join(config['out_path'], config['out_dir'], "map/star/{sample}Aligned.out.bam.bai")
    conda:
        "../envs/CIRC2.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_star_bam_index.txt")
    shell:
        """
        samtools sort {input.bam} -o {input.bam}.sorted.bam
        samtools index {input.bam}.sorted.bam {output.bam_index}
        """

rule circ_parse_star:
    input:
        junction=path.join(config['out_path'], config['out_dir'], "map/star/{sample}Chimeric.out.junction")
    output:
        bsj=path.join(config['out_path'], config['out_dir'], "CIRCexplorer2/star/{sample}_bsj.bed")
    log:
        parse = path.join(config['out_path'],config['out_dir'],"logs/CIRCexplorer2/star/{sample}_parse.log")
    conda:
        "../envs/CIRC2.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_circ_parse_star.txt")
    shell:
        '''
        CIRCexplorer2 parse -t STAR {input.junction} -b {output.bsj} > {log[parse]}
        '''

# 5. CIRCexplorer2 annotation step
rule circ_anno_star:
    input:
        gtf2anno_output=rules.gtf2ano2.output,
        fa=path.join(config['ref_path'], config['ref']),
        bsj=path.join(config['out_path'], config['out_dir'], "CIRCexplorer2/star/{sample}_bsj.bed")
    output:
        txt=path.join(config['out_path'], config['out_dir'], "CIRCexplorer2/star/{sample}")
    log:
        anno=path.join(config['out_path'], config['out_dir'], "logs/CIRCexplorer2/star/{sample}_annotate.log")
    conda:
        "../envs/CIRC2.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_circ_anno_star.txt")
    shell:
        '''
        CIRCexplorer2 annotate -r {input.gtf2anno_output} -g {input.fa} -b {input.bsj} -o {output.txt}  > {log.anno}
        '''
