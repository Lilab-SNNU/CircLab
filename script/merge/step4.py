def circ_name(fasta_file, output_file):
    with open(fasta_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip()[1:].split('(')[0]
                chrom, position = seq_id.split(':')
                start, end = position.split('-')
                out.write(f"{chrom}\t{start}\t{end}\n")

fasta_files = snakemake.input[:4]
output_files = snakemake.output[:4]

for fasta_file, output_file in zip(fasta_files, output_files):
    circ_name(fasta_file, output_file)