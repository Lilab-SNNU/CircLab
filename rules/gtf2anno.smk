from os import path
configfile: "config/config.yaml"

rule gtf2ano1:
    input:
        gtf=config['gtf']
    output:
        temp(path.join(path.dirname(config['gtf']),"annotation0.txt"))
    conda:
        "../envs/CIRC2.yaml"
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/gtf2ano1.txt")
    shell:
        "gtfToGenePred -genePredExt {input[0]} {output[0]}" #envs中加入了ucsc-gtfToGenePred

rule gtf2ano2:
    input:
        rules.gtf2ano1.output
    output:
        path.join(config['out_path'],config['out_dir'],"CIRCexplorer2/annotation.txt")
    benchmark:
        path.join(config['out_path'],config['out_dir'],"benchmark/gtf2ano2.txt")
    script:
        "../script/gtf2ano.py"