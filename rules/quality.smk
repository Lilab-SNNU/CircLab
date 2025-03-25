"""
quality
"""
from os import path
#if config["sample_format"] == "gz":
#print(path.join(config['out_path'],config['out_dir'],"results/fastp/{sample}_1.fq"))
rule fastp:
    input:
        input_fastq
    output:
        path.join(config['out_path'],config['out_dir'],"fastp/{sample}_1.fq"),
        path.join(config['out_path'],config['out_dir'],"fastp/{sample}_2.fq"),
        json=path.join(config['out_path'],config['out_dir'],"logs/fastp/{sample}.json")
    log:
        path.join(config['out_path'],config['out_dir'],"logs/fastp/{sample}.log")
    threads:
        3
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} "
        "-j {output.json} -h /dev/null 2> {log}"

# elif config["sample_format"] == "not_gz":
#     rule fastp:
#         input:
#             input_fastq
#         output:
#             "data/deal.fastq/A_1.fq",
#             "data/deal.fastq/A_2.fq",
#         shell:
#             "trim_galore --paired --quality 20  --length 20  -o data/deal.fastq {input[0]} {input[1]}"
