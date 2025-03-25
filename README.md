# CircLab
A comprehensive software tool for circRNA detection,  function prediction, functional enrichment, and visualization.

Install
-----------
#### - Open dist to download the .whl file and install it by pip
```bash
  pip install circlab-0.1.0-py3-none-any.whl
```
#### - Download whole project by git
```bash
  git clone https://github.com/LiLab-SNNU/CircLab.git
```
Use **Circlab -h** to check.

Use
----

1. ​**CircRNA Detection**
- **We construct a circRNA detection pipeline by snakemake, there are 20 rules in this pipeline which contain fastp, four algorithms to detect, merge.**
   - Integrates multiple detection algorithms (e.g., STAR+CIRCexplorer2, BWA+CIRI2, Bowtie2+find_circ).  
   - Reduces false positives by combining results from multiple tools.  
   - Supports both single-end and paired-end sequencing data.
  ```bsh
     CircLab detection -h 
     CircLab detection -i -g data/genome/aragenome.fa -i data/raw.fastq/Ath/Ath_simu_1.fq,data/raw.fastq/Ath/Ath_simu_2.fq -a data/gene.gtf -o outfile 
  ```

2. ​**Functional Prediction**  
   - Predicts circRNA-miRNA-mRNA interaction networks.
   
    ```bsh
       CircLab function miRNA
    ```
   - Identifies conserved circRNAs across species using sequence alignment and conservation matrix analysis.
   
    ```bsh
       CircLab function conservation 
    ```

3. ​**Visualization**  
   - Provides six visualization functions, including circRNA length distribution, junction reads, tissue-specific expression, and codon composition.  

4. ​**Gene Annotation**  
   - Performs GO and KEGG enrichment analysis for circRNA parent genes.  
   - Outputs structured tables and visualizations for downstream analysis.  
