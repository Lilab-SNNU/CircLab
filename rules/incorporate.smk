from os import path
configfile: "config/config.yaml"
# 将四个软件的circRNA去除冗余
rule merge1:
    input:
        CIRC2_bwa=path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}"),
        CIRC2_star=path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}"),
        CIRI2=path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}"),
        findcirc=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.txt")
    output:
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_circname.bed")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_circname.bed")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_circname.bed")),
        temp(path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_circname.bed"))
    log:
        merge = (path.join(config['out_path'],config['out_dir'],"logs/{sample}_merge.txt"))
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_merge1.txt")
    shell:
        '''
        script/merge/step1.sh {input[CIRC2_bwa]} {input[CIRC2_star]} {input[CIRI2]} {input[findcirc]} \
         {output[0]} {output[1]} {output[2]} {output[3]} > {log[merge]}
        '''

rule merge2:
    input:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_circname.bed"),
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_circname.bed"),
        path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_circname.bed"),
        path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_circname.bed")
    output:
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_circname_filtered.bed")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_circname_filtered.bed")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_circname_filtered.bed")),
        temp(path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_circname_filtered.bed"))
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_merge2.txt")
    script:
        "../script/merge/step2.py"

rule merge3:
    input:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_circname_filtered.bed"),
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_circname_filtered.bed"),
        path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_circname_filtered.bed"),
        path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_circname_filtered.bed"),
        ref=path.join(config['ref_path'],config['ref'])
    output:
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_circseq.fasta")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_circseq.fasta")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_circseq.fasta")),
        temp(path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_circseq.fasta"))
        #"results1/findcirc/Ara_sample_circseq.fasta"
    conda:
        "../envs/merge.yaml"
    log:
        merge = (path.join(config['out_path'],config['out_dir'],"logs/{sample}_merge2.txt"))
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_merge3.txt")
    shell:
        #"script/merge/step3.sh {input.ref} {input[0]} {input[1]} {input[2]} {output[0]} {output[1]} {output[2]}"
        '''
        script/merge/step3.sh {input.ref} {input[0]} {input[1]} {input[2]} {input[3]} \
        {output[0]} {output[1]} {output[2]} {output[3]} > {log[merge]}
        '''

rule merge4:
    input:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_circseq.fasta"),
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_circseq.fasta"),
        path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_circseq.fasta"),
        path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_circseq.fasta")
    output:
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_all_circ_name.txt")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_all_circ_name.txt")),
        temp(path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_all_circ_name.txt")),
        temp(path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_all_circ_name.txt"))
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_merge4.txt")
    script:
        "../script/merge/step4.py"

rule merge5:
    input:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_all_circ_name.txt"),
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_all_circ_name.txt"),
        path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_all_circ_name.txt"),
        path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_all_circ_name.txt"),
        CIRC2_bwa=path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}"),
        CIRC2_star=path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}"),
        CIRI2=path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}"),
        findcirc=path.join(config['out_path'],config['out_dir'],"findcirc/{sample}.txt")
    output:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_all_recirc.txt"),
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_all_recirc.txt"),
        path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_all_recirc.txt"),
        path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_all_recirc.txt")
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_merge5.txt")
    script:
        "../script/merge/step5.py"

rule incorporate:
    input:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/bwa/{sample}_all_recirc.txt"),
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/star/{sample}_all_recirc.txt"),
        path.join(config['out_path'],config['out_dir'],"CIRI2/{sample}_all_recirc.txt"),
        path.join(config['out_path'],config['out_dir'],"findcirc/{sample}_all_recirc.txt")
    output:
        path.join(config['out_path'],config['out_dir'],"{sample}_recirc.txt")
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/{sample}_incorporate.txt")
    script:
        "../script/incorporate.py"

# rule all:
#     input:
#         expand("results1/all_recirc")
#
# rule merge1:
#     input:
#         CIRC2_bwa = "results1/CIRCexplorer2/bwa/Ara_sample",
#         CIRC2_star = "results1/CIRCexplorer2/star/Ara_sample",
#         CIRI2 = "results1/CIRI2/Ara_sample",
#         findcirc = "results1/findcirc/Ara_sample.txt"
#     output:
#         temp("results1/CIRCexplorer2/bwa/Ara_sample_circname.bed"),
#         temp("results1/CIRCexplorer2/star/Ara_sample_circname.bed"),
#         temp("results1/CIRI2/Ara_sample_circname.bed"),
#         temp("results1/findcirc/Ara_sample_circname.bed")
#     log:
#         merge = "results1/logs/merge/merge.log"
#     shell:
#         "script/merge/step1.sh {input[CIRC2_bwa]} {input[CIRC2_star]} {input[CIRI2]} {input[findcirc]} \
#          {output[0]} {output[1]} {output[2]} {output[3]} > {log[merge]}"
#
# rule merge2:
#     input:
#         "results1/CIRCexplorer2/bwa/Ara_sample_circname.bed",
#         "results1/CIRCexplorer2/star/Ara_sample_circname.bed",
#         "results1/CIRI2/Ara_sample_circname.bed",
#         "results1/findcirc/Ara_sample_circname.bed"
#     output:
#         temp("results1/CIRCexplorer2/bwa/Ara_sample_circname_filtered.bed"),
#         temp("results1/CIRCexplorer2/star/Ara_sample_circname_filtered.bed"),
#         temp("results1/CIRI2/Ara_sample_circname_filtered.bed"),
#         temp("results1/findcirc/Ara_sample_circname_filtered.bed")
#     script:
#         "../script/merge/step2.py"
#
# rule merge3:
#     input:
#         "results1/CIRCexplorer2/bwa/Ara_sample_circname_filtered.bed",
#         "results1/CIRCexplorer2/star/Ara_sample_circname_filtered.bed",
#         "results1/CIRI2/Ara_sample_circname_filtered.bed",
#         "results1/findcirc/Ara_sample_circname_filtered.bed",
#         ref=path.join(config['ref_path'],config['ref'])
#     output:
#         temp("results1/CIRCexplorer2/bwa/Ara_sample_circseq.fasta"),
#         temp("results1/CIRCexplorer2/star/Ara_sample_circseq.fasta"),
#         temp("results1/CIRI2/Ara_sample_circseq.fasta"),
#         temp("results1/findcirc/Ara_sample_circseq.fasta"),
#         temp("results1/CIRCexplorer2/bwa/Ara_sample_circseq.fasta.clstr"),
#         temp("results1/CIRCexplorer2/star/Ara_sample_circseq.fasta.clstr"),
#         temp("results1/CIRI2/Ara_sample_circseq.fasta.clstr"),
#         temp("results1/findcirc/Ara_sample_circseq.fasta.clstr")
#     conda:
#         "../envs/merge.yaml"
#     log:
#         merge = "results1/logs/merge.log"
#     shell:
#         '''
#         script/merge/step3.sh {input.ref} {input[0]} {input[1]} {input[2]} {input[3]} \
#         {output[0]} {output[1]} {output[2]} {output[3]} >> {log[merge]}
#         '''
#
# rule merge4:
#     input:
#         "results1/CIRCexplorer2/bwa/Ara_sample_circseq.fasta",
#         "results1/CIRCexplorer2/star/Ara_sample_circseq.fasta",
#         "results1/CIRI2/Ara_sample_circseq.fasta",
#         "results1/findcirc/Ara_sample_circseq.fasta"
#     output:
#         temp("results1/CIRCexplorer2/bwa/Ara_sample_all_circ_name.txt"),
#         temp("results1/CIRCexplorer2/star/Ara_sample_all_circ_name.txt"),
#         temp("results1/CIRI2/Ara_sample_all_circ_name.txt"),
#         temp("results1/findcirc/Ara_sample_all_circ_name.txt")
#     script:
#         "../script/merge/step4.py"
#
# rule merge5:
#     input:
#         "results1/CIRCexplorer2/bwa/Ara_sample_all_circ_name.txt",
#         "results1/CIRCexplorer2/star/Ara_sample_all_circ_name.txt",
#         "results1/CIRI2/Ara_sample_all_circ_name.txt",
#         "results1/findcirc/Ara_sample_all_circ_name.txt",
#         CIRC2_bwa= "results1/CIRCexplorer2/bwa/Ara_sample",
#         CIRC2_star="results1/CIRCexplorer2/star/Ara_sample",
#         CIRI2="results1/CIRI2/Ara_sample",
#         findcirc="results1/findcirc/Ara_sample.txt"
#     output:
#         "results1/CIRCexplorer2/bwa/Ara_sample_all_recirc.txt",
#         "results1/CIRCexplorer2/star/Ara_sample_all_recirc.txt",
#         "results1/CIRI2/Ara_sample_all_recirc.txt",
#         "results1/findcirc/Ara_sample_all_recirc.txt"
#     script:
#         "../script/merge/step5.py"
#
# rule incorporate:
#     input:
#         "results1/CIRCexplorer2/bwa/Ara_sample_all_recirc.txt",
#         "results1/CIRCexplorer2/star/Ara_sample_all_recirc.txt",
#         "results1/CIRI2/Ara_sample_all_recirc.txt",
#         "results1/findcirc/Ara_sample_all_recirc.txt"
#     output:
#         "results1/all_recirc"
#     script:
#         "../script/incorporate.py"