import subprocess
import argparse
import os
import sys
from os import path
import yaml
import time
from collections import defaultdict


def run_snakemake(genome, anno, input_file, output_dir, cores):
    start_time = time.time()
    snakemake_cmd = [
        'snakemake',
        '--snakefile', 'pipes/Snakefile',
        '--cores', str(cores),
        '-p',
        '--use-conda',
        '--keep-going'
        # '--unlock'
        # '--forceall',
    ]
    end_time = time.time()
    print('=' * 100 + f'\nsnakemake pipeline : {" ".join(snakemake_cmd)}\n')
    subprocess.run(snakemake_cmd)
    end_time = time.time()
    print(f"Detection module runtime: {end_time - start_time:.2f} seconds")

def snakemake_dag():
    snakemake_cmd = [
        'snakemake',
        '--snakefile', 'pipes/Snakefile',
        '--dag | dot -Tsvg > dag.svg',
    ]

def load_config(config_path):
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def generate_config(genome, anno, input_file, output_file):
    print(get_time(), f"Running detection module...")
    print(f"--Genome: {genome}")
    print(f"--Annotation File: {anno}")
    print(f"--Output Path: {output_file}")

    # Create config dictionary
    if ',' in input_file:
        samples = input_file.split(',')
        print(f'Paired-End reads: {samples[0]},{samples[1]}')
    else:
        samples = [input_file]
        print(f'Single-read {samples[0]}')

    # Initialize the dictionary that saves the sample name and file list
    samples_dict = defaultdict(list)
    for file_path in input_file.split(','):
        filename = os.path.basename(file_path)
        sample_name = filename.rsplit('_', 1)[0]
        samples_dict[sample_name].append(file_path)

    samples_path = path.dirname(samples[0])
    samples_dict = dict(samples_dict)
    sample_format = ['gz', 'not_gz']

    ref_path = path.dirname(genome)
    ref = path.basename(genome)
    output_path = path.dirname(output_file)
    out_dir = path.basename(output_file)

    config_dict = {
        'samples': samples_dict,
        'sample_format': sample_format,
        'samples_path': samples_path,
        'suffix': 'fq',
        'f_suffix': '_1.fq',
        'r_suffix': '_2.fq',
        'ref': ref,
        'ref_path': ref_path,
        'gtf': anno if anno else '',
        'out_path': output_path,
        'out_dir': out_dir,
    }

    # Ensure output directory exists
    if not os.path.exists(output_file):
        os.makedirs(output_file)

    # Write YAML config file to the specified location
    with open('config/config.yaml', 'w') as f:
        yaml.dump(config_dict, f)
    print(f"Configuration file has been saved to config/config.yaml")


def get_time():
    nowtime = time.strftime('[%Y-%m-%d %H:%M:%S]', time.localtime())
    return nowtime


def detection(config_output_path=None, genome=None, anno=None, input_file=None, output_file=None, cores=10):
    # If the config file exists, load it, otherwise generate one
    if config_output_path:
        if os.path.exists(config_output_path):
            print(f"Loading configuration from {config_output_path}...")
            config_dict = load_config(config_output_path)
        else:
            print(f"Config file {config_output_path} not found, generating a new one...")
            generate_config(genome, anno, input_file, output_file)
    else:
        print(f"Generating config file...")
        generate_config(genome, anno, input_file, output_file)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Detection circRNA detect module",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Add command line arguments
    parser.add_argument('-g', '--genome', help="Genome reference file")
    parser.add_argument('-a', '--anno', help="Annotation file in GTF/GFF format", required=False)
    parser.add_argument('-i', '--input', dest='input_file', help='Input file(s) name, fastq/fq file')
    parser.add_argument('-o', '--output', dest='output_dir', help="Output directory")
    parser.add_argument('--config', dest='config_output', help="Path to an existing config.yaml file")
    parser.add_argument('-c', '--cores', type=int, default=12, help="Number of cores to use for Snakemake (default=12)")

    args = parser.parse_args(argv)

    # If --config is provided, use the config file, otherwise generate one based on input
    if args.config_output:
        detection(config_output_path=args.config_output)
    else:
        if not args.genome or not args.input_file or not args.output_dir:
            print("Error: Missing required arguments for auto-config generation.")
            sys.exit(1)
        detection(genome=args.genome, anno=args.anno, input_file=args.input_file, output_file=args.output_dir)

    # After detecting, run Snakemake with the user-defined number of cores
    run_snakemake(args.genome, args.anno, args.input_file, args.output_dir, args.cores)
    #snakemake_dag()
    print(get_time(), f"Finished detection, please see the output file “{args.output_dir}” for details.")

    output_path = os.path.abspath(args.output_dir)
    if os.path.exists(output_path) and os.listdir(output_path):
        print(get_time(), "Detection output exists. Do you want to proceed with annotation?(y/n):",end="")
        user_input = input().strip().lower()
        if user_input == 'y':
            run_annotation(output_path)
    else:
        print(get_time(), "No valid output detected.Annotation will not proceed.")


if __name__ == "__main__":
    main(sys.argv[1:])
