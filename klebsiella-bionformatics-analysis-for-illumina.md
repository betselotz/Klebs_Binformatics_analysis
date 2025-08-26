# Bioinformatics Analysis of Klebsiella pneumoniae of Illumina# 

---  

###### **_Group_members_**: [Betselot Z Ayano ](https://github.com/), [], [Gilbert Kibet] & [Ouso Daniel]

---

## Table of Contents
  - [Overview](#overview)
  - [Learning Objectives](#learning-objectives)
  - [Prerequisites](#background)
  - [prerequisites](#prerequisites)
  - [Scope of the Tutorial](#scope-of-the-tutorial)
  - [Set-up](#set-up)
      - [Workshop Environment](#workshop-environment)
      - [Logging into the HPC](#logging-into-the-hpc)
      - [Compute Node](compute-node)
      - [Project organisation](#project-organisation)
  - [Bioinformatics Analysis](#bioinformatics-analysis)
      - [About the Sample](#about-the-sample)
      - [Step 1: Load required modules](#step-1-load-required-modules)
      - [Step 2: Retrieve reference genome in GenBank and FASTA format](#step-2-retrieve-reference-genome-in-genBank-and-fasta-format)
      - [Step 3: Data Quality Assessment](#step-3-data-quality-assessment)
      - [Step 4: Genome Assembly](#step-4-genome-assembly)
      - [Step 5: Genome Annotation](#step-5-genome-annotation)
      - [Step 6: Pathogen Typing](#step-6-pathogen-typing)
        - [MLST](#mlst)
        - [MLST Output Format](#mlst-output-format)
        - [MLST Results Interpretation](#mlst-results-interpretation)
        - [Visualising MLST Results](#visualising-mlst-results)
        - [MLST Interpretation Limitations](#mslt-interpretation-limitations)
        - [Serotyping using Kaptive](#serotyping-using-kaptive)
        - [K antigen and locus](#k-antigen-and-locus)
        - [O antigen and locus](#0-antigen-and-locus)
      - [Step 7: AMR genes detection](#step-7-amr-genes-detection)
          - [AMR genes detection using
            ResFinder](#amr-genes-detection-using-resfinder)
          - [Batch AMR detection](#batch-amr-detection)
          - [AMR Detection using CARD and RGI web resources](#amr-detection-using-card-and-rgi-web-resources)
          - [AMR genes detection using
            AMRFinder](#amr-genes-detection-using-amrfinder)
          - [Output Format](#output-format)
      - [Step 8: Variant Calling and Consensus Assemblies](#step-8-variant-calling-and-consensus-assemblies)
          - [Fast Bacterial Variant Calling with Contigs](#fast-bacterial-variant-calling-with-contigs)
          - [Snippy Outputs](#snippy-outputs)
          - [Visualising Snippy Variants](#visualising-snippy-variants)
          - [Build Core and Whole Genome Alignments from Snippy Output](#build-core-and-whole-genome-alignments-from-snippy-output)
          - [Run Snippy-core](#run-snippy-core)
          - [Snippy Core Outputs](#snippy-core-outputs)
          - [Cleanup the Snippy SNP Alignment Intermediates](#cleanup-the-snippy-snp-alignment-intermediates)
          - [Compute Pairwise SNP Distances](#compute-pairwise-snp-distances)
    - [Step 9: Phylogenetic Analysis](#step-9-phylogenetic-analysis)
        - [Phylogenetic Analysis of Gubbins
          Output](#phylogenetic-analysis-of-gubbins-output)
        - [Maximum likelihood phylogenetic inference](#maximum-likelihood-phylogenetic-inference)
    - [Step 10: Visualization](#step-10-visualization)
        - [Visualize the phylogeny alongside typing, antibiotic resistance or
          epidemiological
          data](#visualize-the-phylogeny-alongside-typing-antibiotic-resistance-or-epidemiological-data)
        - [Visualization using other tools](#visualization-using-other-tools)
    


## Introduction – Klebsiella pneumoniae

    **Opportunistic Pathogen**

        Gram-negative bacterium

        Causes hospital-acquired infections: pneumonia, bloodstream, urinary tract

    **Antimicrobial Resistance (AMR)** 

        Produces ESBLs and carbapenemases (KPC, NDM, etc.)

        Limits treatment options, including last-resort carbapenems

    **Drug-Resistant Strains**

        MDR (multidrug-resistant) & XDR (extensively drug-resistant)

        Spread via plasmid-mediated genes and successful clonal lineages

    **Genomic Insights**

        Whole-genome sequencing (WGS): tracks resistance genes & plasmids

        Phylogenomics: reveals clonal expansion and transmission routes

        Pangenome analysis: highlights core vs. accessory genes driving adaptation

        Mobile genetic elements (MGEs): key role in resistance gene dissemination

    **Public Health Impact**

        Rapid global dissemination in healthcare settings

        Major challenge for infection prevention and treatment
---


## Setting up the Bioinformatics analysis environment

### Set up directories  

 - Before starting the analysis, ensure that you are logged into the HPC, create an interactive session on the assigned compute node, and change directory to the project folder which is `ACDC_AMR2025`.  
 - First log in to the HPC.   

```
ssh <user_name>@hpc.ilri.cgiar.org
```
  
>Compute05  
```
interactive -w compute05 -c 2 -J amr-surveillance -p batch
```
>Compute06  
```
interactive -w compute06 -c 2 -J amr-surveillance -p batch
```
 - Setting up the project directory structure.  

```bash
mkdir -p /var/scratch/$USER/illumina_downloads
cd /var/scratch/$USER/illumina_downloads
```

 - Within our project directory let us create a structured project directories:

```bash
mkdir -p \
results/illumina/ecoli/{fastqc,fastp,fastq-scan,shovill,prokka,resfinder,amrfinder,mlst,tmp/{shovill,prokka,resfinder,amrfinder}}
```

 - Then let us link some data.  
 
Our analysis data is in the group leader diretory `ln -sf /var/scratch/user3/illumina_downloads/`. We will create a symbolic link of the folders with data to current directory.

```bash
ln -sf /var/scratch/user3/illumina_downloads/[dps]* .
```

### Load modules

```
module load fastqc/0.11.9
module load fastp/0.22.0
module load bbmap/38.95
module load shovill/1.1.0
module load prokka/1.14.6
module load resfinder/4.6.0
module load mlst/2.23.0
module load amrfinder/4.0.22
```

## Bioinformatics Analysis
we make directory called data inside illumina downloads to download our raw data in it. 
```bash
mkdir data
cd data
```
### Download data from NCBI using sra explorer  
### Option 1: SRA
 - Go to SRA (https://sra-explorer.info/#)  
 - search: `SRR28370701`  
 - checkboxes to check in SRA:  
`   - `WGS of Klebsiella pneumoniae isolate`;
    - `add 1 to collection`;
    - `1 data saved`;
    - `Bash script for downloading FastQ files` 

```bash
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR283/001/SRR28370701/SRR28370701_1.fastq.gz -o SRR28370701_WGS_of_Klebsiella_pneumoniae_isolate_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR283/001/SRR28370701/SRR28370701_2.fastq.gz -o SRR28370701_WGS_of_Klebsiella_pneumoniae_isolate_2.fastq.gz
```


### Assessing Read Quality using fastqc before quality trimming  

The **raw `FASTQ`** files have sequences as generated by the sequencer. They include **poor quality reads**, **adapters**, **PhiX** and some reads may be **duplicates**. We need to check the quality of these suquences and clean up if we need to.  
We use our first module - **`FASTQC`**, a bioinformatics software used to analyse quality of raw FASTQ format reads and generate visual plots. The report generated by this step will have information on **Number of reads**,**Sequence Quality**, **Sequence length**, **GC content**, **Adapter content**, **duplication rate** and others. Run the command below to execute fastqc on our reads:

```bash
fastqc \
    -o ./results/illumina/klebs/fastqc \
    --noextract \
    -t 2 \
    --nogroup \
    ./data/klebs/SRR28370701_1.fastq.gz \
    ./data/klebs/SRR28370701_2.fastq.gz
```
Takes less than 1 minute

 - Now download the results of the fastqc command to your local laptop for evaluation. The results are in a `HTML` file.  
 - First copy to `~/`(home).

```bash
cp ./results/illumina/klebs/fastqc/*html ~/
```
 - To copy the report to a directory on **`local pc`**

```bash
rsync -avP --partial <USERXX>@hpc.ilri.cgiar.org:~/SRR28370701*.html ~/
```
 - Opened the HTML file and explore.


### Quality Trimming fastq files with fastp and Trims adapter sequences

After assessing the quality, we will proceed and do some Quality Control (QC). With **`fastp`** module we will trim reads with qualities less than 20 phred score, remove adapters, remove duplicates and remove PhiX if any.

```bash
fastp \
    --in1 ./data/klebs/SRR28370701_1.fastq.gz \
    --in2 ./data/klebs/SRR28370701_2.fastq.gz \
    --out1 ./results/illumina/klebs/fastp/SRR28370701_trim_R1.fastq.gz \
    --out2 ./results/illumina/klebs/fastp/SRR28370701_trim_R2.fastq.gz \
    --detect_adapter_for_pe \
    --thread 2 \
    --json ./results/illumina/klebs/fastp/SRR28370701.fastp.json \
    --html ./results/illumina/klebs/fastp/SRR28370701.fastp.html \
    --cut_mean_quality 20 \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --qualified_quality_phred 25 \
    --unqualified_percent_limit 40 \
    --length_required 20 \
    2>&1 | tee ./results/illumina/klebs/fastp/SRR28370701.fastp.log
```


 - `fastp` also generates a report in `HTML` format. Let us download it and explore. First copy it to home (`~/`)

```bash
cp ./results/illumina/klebs/fastp/SRR28370701.fastp.html ~/
```
 - Then download:

```bash
rsync -avP --partial <USERXX>@hpc.ilri.cgiar.org:~/SRR28370701.fastp.html ~/
```
 - Examine the report.
 - Links available online: [SRR25008769.fastp.html](https://hpc.ilri.cgiar.org/~gkibet/AMR-Genomic-Surveillance/SRR25008769.fastp.html)

### De novo assembly pipeline for Illumina paired reads
Assemble bacterial isolate genomes from Illumina paired-end reads

Shovill is a pipeline which uses SPAdes at its core, but alters the steps before
and after the primary assembly step to get similar results in less time. Shovill
also supports other assemblers like SKESA, Velvet and Megahit, so you can take
advantage of the pre- and post-processing the Shovill provides with those too.


>**Note**
Shovill is for isolate data only, primarily small haploid organisms. 
It will NOT work on metagenomes or larger genomes. 
Please use Megahit directly instead.


#### Main steps

1. Estimate genome size and read length from reads (unless --gsize provided) - mash
2. Reduce FASTQ files to a sensible depth (default --depth 100) - seqtk
3. Trim adapters from reads (with --trim only) - trimmomatic
4. Conservatively correct sequencing errors in reads - lighter
5. Pre-overlap ("stitch") paired-end reads - flash
6. Assemble with SPAdes/SKESA/Megahit with modified kmer range and PE + long SE reads
7. Correct minor assembly errors by mapping reads back to contigs (bwa-mem + pilon)
8. Remove contigs that are too short, too low coverage, or pure homopolymers
9. Produce final FASTA with nicer names and parseable annotations


```
export TMPDIR="./results/illumina/ecoli/tmp/shovill"
```

```
shovill \
  --R1 ./results/illumina/klebs/fastp/SRR28370701_trim_R1.fastq.gz \
  --R1 ./results/illumina/klebs/fastp/SRR28370701_trim_R2.fastq.gz \
  --gsize 5249449 \
  --outdir ./results/illumina/klebs/shovill/SRR28370701 \
  --assembler skesa \
  --minlen 500 \
  --mincov 2 \
  --force \
  --keepfiles \
  --depth 0 \
  --noreadcorr \
  --namefmt "SRR28370701_%05d" \
  --cpus 4 \
  --ram 16 \
  --tmpdir $TMPDIR
```

```
mv ./results/illumina/klebs/shovill/SRR28370701/contigs.fa ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa
```

### Assembly evaluation
Assembly evaluation

We need to verify the quality of the resulting assembly before any downstream processes. We will use the utility script below:

``` bash
stats.sh in=./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa
```
>**Note**
Unfortunately, the N50 and L50 values generated by stats.sh are switched. N50
should be a length and L50 should be a count. The results table below shows the
corrected values based on stats.sh outputs.


``` bash
module avail seqkit
```
``` bash
module load seqkit/0.11.0 
```
``` bash
seqkit -h
```
``` bash
seqkit fx2tab -nl ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa
```
We can use another tool assembly-scan to generate summary statistics of the assembly.

``` bash
module avail seqkit
```
``` bash
module load assembly-scan/1.0.0
```
``` bash
assembly-scan /var/scratch/user3/illumina_downloads/results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa \
  --transpose \
  | tee /var/scratch/user3/illumina_downloads/results/illumina/klebs/shovill/SRR28370701/SRR28370701-assembly-scan.tsv
```

Compute the GC content of the assembly using the output from assembly-scan output
``` bash
grep \
 'contig_percent_[cg]' \
 ./results/illumina/klebs/shovill/SRR28370701/SRR28370701-assembly-scan.tsv | \
 awk -F '\t' '{sum+=$3} END {print "GC%=",sum}'
 awk -F '\t' '{sum+=$3} END {print "GC%=",sum}'
```
Transfer assembly to local computer

Copy the assembly to your home directory on hpc
``` bash
rsync -avP \
    --partial \
    ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa \
    ~/
```
``` bash
rsync -avP --partial user3@hpc.ilri.cgiar.org:~/SRR28370682.fa ~/AMR_training/group1/
```


# Step 5: Genome Annotation

After assembling our genome, we need to characterise it through annotation to enable inferencing. Prokka (Prokaryotic annotation) is a pipeline tool designed for the rapid annotation of prokaryotic genomes, including bacteria and archaea. After Prokka annotation, tools like ABRicate or RGI can be run on the annotated genome to identify resistance genes in their proper genomic context.

Prokka respects the default Linux temporary directory `--TMPDIR`. Therefore we
assign one if none exists.  

However, because we are all running the same pipeline and we don't want to
overwrite each other's files, we will use our project-specific temporary
directories.

```
export TMPDIR=./results/ont/klebsiella/tmp/prokka/
```

The tool also comes with its own database which is stored in the variable
`PROKKA_DBDIR`. To point it to a custom DB, the variable can be tailored and
exported as follows `export PROKKA_DBDIR=<path/to/custom/db>`


```
prokka \
    --evalue 1e-09 \
    --coverage 80 \
    --centre ILRI \
    --cpus 2 \
    --prefix SRR28370701 \
    --locustag SRR28370701 \
    --proteins ./databases/prokka/proteins.faa \
    --force \
    --outdir ./results/illumina/klebs/prokka \
    ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa
```

Prokka generates multiple output files in standard bioinformatics formats:

|File Extension | Format | Description|
|--- | --- |---|
`.gff` | GFF3 | General feature format (coordinates of all features)
`.gbk` | GenBank | Standard sequence with annotations format
`.fna` | FASTA | Nucleotide FASTA file of the input contig sequences
`.faa` | FASTA | Protein FASTA file of the translated CDS sequences
`.ffn` | FASTA | Nucleotide FASTA file of all feature sequences
`.sqn` | Sequin | ASN.1 format for GenBank submission
`.fsa` | FASTA | Nucleotide FASTA file of the input contig sequences with annotations
`.tbl` | Feature Table | Feature table for GenBank submission
`.txt` | Text | Statistics of the annotation run

# Add additional genomes from pathogenwatch

Here we will use 11 Klebs isolates collected in Kenya between January 14 and January 31, 2019
https://pathogen.watch/genomes/all?country=ke&genusId=570&maxDate=2019-01-31T20%3A59%3A59.999Z&minDate=2018-12-31T21%3A00%3A00.000Z&sort=date&speciesId=573

```
mkdir -p ./pathogenwatch/klebs/assemblies-to-test
```


1. add our genome assembly to the directory for global assemblies

```
rsync -avP --partial \
    ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa \
    ./pathogenwatch/klebs/assemblies-to-test/SRR28370701.fasta
```

2. add reference genome assembly to the directory for global assemblies

```
rsync -avP --partial \
    ./genomes/klebs/GCF_000016305.1_ASM1630v1_genomic.fna \
    ./pathogenwatch/klebs/assemblies-to-test/Reference.fasta
```

3. select assemblies to test and add to the directory for global assemblies

```
rsync -avP --partial \
    ./genomes/klebs/*.fna \
    ./pathogenwatch/klebs/assemblies-to-test/
```

# Step 6: Pathogen Typing

Understanding the differences and relationships within circulating pathogens is an important aspect of genomics epidemiology. Whole-genome comparison has better resolution for pathogen characterisation than fragments per genome equivalent (FPGE) or gene-based or multi-locus sequence typing (MLST). Distances between genomes can be compared based a reference or _de novo_; _k-mer_-composition-based and _core-genome_ assembly-based.

#### MLST 

Multi-Locus Sequence Typing (MLST) is an essential tool in analyzing multiple antimicrobial resistant (MAR) bacterial isolates. It provides a standardized approach to characterizing bacterial strains based on the sequences of multiple housekeeping genes. It focuses on bacterial population structure and evolutionary relationships.

MLST is a molecular typing method that characterizes bacterial isolates by sequencing internal fragments (typically 450-500 bp) of multiple housekeeping genes (usually 7-8 genes). Each unique sequence for a gene is assigned an allele number, and the combination of allele numbers defines the sequence type (ST) of an isolate.

``` bash
MLST_DB=$(find /export/apps/mlst/2.23.0/db -name "mlst.fa" | sed 's=blast/mlst.fa==')
```

``` bash
mlst \
    --threads 2 \
    --blastdb $MLST_DB/blast/mlst.fa \
    --datadir $MLST_DB/pubmlst \
    --scheme klebsiella \
    --minid 100 \
    --mincov 10 \
    --minscore 50 \
    ./results/illumina/klebs/prokka/SRR28370701.fna \
    > ./results/illumina/klebs/mlst/SRR28370701.tsv

```


**Batch MLST**

``` bash
for fn in ./pathogenwatch/klebs/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running mlst on: $sample - $fn"

    mlst \
        --threads 2 \
        --blastdb $MLST_DB/blast/mlst.fa \
        --datadir $MLST_DB/pubmlst \
        --scheme klebsiella \
        --minid 100 \
        --mincov 10 \
        --minscore 50 \
        $fn \
        > ./results/illumina/klebs/mlst/${sample}.tsv
done

```


##### MLST Output Format

The standard output of MLST analysis is a tabular plain text file with the following columns.

Column | Description
--- | ---
FILE | Input sequence name
SCHEME | The specific bacterial species or genus
ST | Sequence type number
Allelic profile | Depends on how many gene are used in the scheme

> **Note:** When reporting MLST results the scheme used for the profiling must be provided for accurate results interpretation.

<details>
  <summary>
    Click to toggle <b style="color:blue">BIGSdb platform for assigning STs</b>
  </summary>
  <p>
    <a href="https://bigsdb.pasteur.fr/klebsiella" target="_blank">BIGSdb</a> is curated by the Institut Pasteur.
  </p>
</details>

##### MLST Result Interpretation
The most important information of the results in the ST:
- **Know STs:** which matches a database hit a number is assigned.
- **Novel allele combinations:** may be represented as `?` or `novel` if no known database profile is matched
- **Incomplete matches:** may be shown as `ST-like` or with an `*` (asterisk) if most but not all allele profiles match a known database ST
- **Coverage, identity and depth:** are importance to consider in interpratation
- Look out for **STs with specific characteristics**, eg., _K. pneumoniae_ ST258 is a major global clone carrying KPC carbapenemases
- **Clonal complexes:** STs sharing alleles at most loci, usually sharing identical alleles at 5 or more loci; reported as `CC` followed by the number of the central/founding ST
- MLST reuslts can indicate evolutionary relationships

##### Visualising MLST Results
A number of visualisation tools are available, examples:
- [**GrapeTree**](https://achtman-lab.github.io/GrapeTree/MSTree_holder.html): creates hierarchical clustering of MLST data
- [**goeBURST**](https://www.phyloviz.net/goeburst/#Software): a classic tool for visualizing MLST data that focuses on identifying clonal complexes

##### MLST Interpretation Limitations
When analyzing MLST results, be aware of certain limitations:

- Limited Resolution: MLST may not distinguish between closely related isolates
- Temporal Dynamics: Does not capture all evolutionary changes over time
- Geographic Bias: Some STs may be overrepresented in databases due to sampling bias
- Horizontal Gene Transfer: May complicate interpretation of evolutionary
  relationships


#### Merge MLST reports

``` bash
cat \
    ./results/illumina/klebs/mlst/*.tsv > \
    ./results/illumina/klebs/mlst/klebs-mlst.txt
```

#### Assign sequence types using web resources

BIGSdb platform curated by the Institute Pasteur
(https://bigsdb.pasteur.fr/klebsiella)


# Serotyping using Kaptive

Like other human pathogens, members of the KpSC produce immunogenic surface
polysaccharides which protect against host defences, phages and desiccation. In
recent years there has been a resurgence of interest in these polysaccharides as
targets for anti-KpSC vaccines 💉. In order to design effective vaccines, we
need to understand the underlying diversity and epidemiology of polysaccharide
variants in the population. We can now achieve this goal at scale by predicting
surface polysaccharide antigens from rapidly growing collections of whole genome
sequences.

Kaptive is a system for surface polysaccharide typing from bacterial genome
sequences.

For each input assembly, Kaptive runs the kaptive.assembly.typing_pipeline which does the following:

- Aligns locus gene nucleotide sequences to the assembly contig sequences using minimap2.
- Identifies the best matching locus type using the scoring algorithm.
- Extracts the locus gene sequences from the assembly contig sequences.
- Predicts the serotype/phenotype based on the gene content.


### K antigen and locus
The capsular polysaccharide (CPS), or K-antigen, is an extracellular layer of polysaccharide chains, anchored to the outer-membrane (OM).
``` bash
mkdir -p ./results/illumina/klebs/kaptive
```

``` bash
kaptive assembly \
  /export/apps/kaptive/3.1.0/lib/python3.10/site-packages/reference_database/Klebsiella_k_locus_primary_reference.gbk \
  ./pathogenwatch/klebs/assemblies-to-test/*.fasta \
  --min-cov 70 \
  -t 2 \
  -o ./results/illumina/klebs/kaptive/kaptive_k_locus.tsv
```

### O antigen and locus
The outer lipopolysaccharide (LPS), or O-antigen, locus is located next to the K-locus. Each KpSC genome contains a single O-locus, for which X variants have been defined to-date. 


``` bash 
kaptive assembly \
  /export/apps/kaptive/3.1.0/lib/python3.10/site-packages/reference_database/Klebsiella_o_locus_primary_reference.gbk \
  ./pathogenwatch/klebs/assemblies-to-test/*.fasta \
  --min-cov 70 \
  -t 2 \
  -o ./results/illumina/klebs/kaptive/kaptive_o_locus.tsv
```  

### Kaptive Web
You can also upload your assemblies to
[Kaptive-Web](https://kaptive-web.erc.monash.edu/). 


### Summary Table

Column | Description
--- | ---- |
`Assembly` | the name of the input assembly, taken from the assembly filename.
`Best match locus` | the locus type which most closely matches the assembly
`Best match type` | the predicted serotype/phenotype
`Match confidence` | a categorical measure of locus match quality.
`Problems` | characters indicating issues with the locus match. An absence of any such characters indicates a very good match. `?` = the match was not in a single piece, possible due to a poor match or discontiguous assembly.The number of pieces will be indicated with an integer. `-` = one or more genes expected in the locus were not found.`+` = one or more extra genes were found in the locus.`*` = one or more expected genes was found but with identity below the minimum threshold (default threshold for KpSC is 82.5%, default for A. baumannii is 85%). `!` = one or more locus genes is truncated
`Identity` | weighted percent identity of the best matching locus to the assembly.
`Coverage` | weighted percent coverage of the best matching locus in the assembly.
`Length discrepancy` | the difference in length between the locus match and the corresponding part of the assembly. Only available if the locus was found in a single piece (i.e. the ? problem character is not used).
`Expected genes in locus` | a fraction indicating how many of the genes in the best matching locus were found in the locus part of the assembly.
`Expected genes in locus, details` | gene names and percent identity and percent coverage for the expected genes found in the locus part of the assembly.
`Missing expected genes` | a string listing the gene names of expected genes that were not found.
`Other genes in locus` | the number of unexpected genes (genes from loci other than the best match) which were found in the locus part of the assembly.
`Other genes in locus, details` | gene names, percent identity and coverage for the other genes found in the locus part of the assembly.
`Expected genes outside locus` | the number of expected genes which were found in the assembly but not in the locus part of the assembly (usually zero)
`Expected genes outside locus, detail` | gene names, percent identity and coverage for the expected genes found outside the locus part of the assembly.
`Other genes outside locus` | the number of unexpected genes (genes from loci other than the best match) which were found outside the locus part of the assembly.
`Other genes outside locus, details` | gene names, percent identity and coverage for the other genes found outside the locus part of the assembly.
`Truncated genes, details` | gene names for the truncated genes found in the assembly.
`Extra genes, details` | Gene names for the extra genes found in the assembly.



# Step 7: AMR genes detection


## AMR genes detection using ResFinder

Now that we have an annotated genome, we can query it for antimicrobial
resistance. There a variety of tools for the task. We will use commonly used
tools whose AMR databases are regularly/frequently.

[**ResFinder**](https://genepi.dk/resfinder) identifies acquired genes and/or finds chromosomal mutations mediating antimicrobial resistance in total or partial DNA sequence of bacteria.

``` bash 
python -m resfinder \
    -ifa ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa \
    -o ./results/illumina/klebs/resfinder/SRR28370701 \
    -s klebsiella \
    --min_cov 0.6 \
    --threshold 0.9 \
    --min_cov_point 0.6 \
    --threshold_point 0.9 \
    --ignore_stop_codons \
    --ignore_indels \
    --acquired \
    --point

```

## Batch AMR detection 

``` bash 
for fn in ./pathogenwatch/klebs/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running ResFinder on: $sample - $fn"

    python -m resfinder \
        -ifa $fn \
        -o ./results/illumina/klebs/resfinder/${sample} \
        -s klebsiella \
        --min_cov 0.6 \
        --threshold 0.9 \
        --min_cov_point 0.6 \
        --threshold_point 0.9 \
        --ignore_stop_codons \
        --ignore_indels \
        --acquired \
        --point
done
```

## AMR Detection using CARD and RGI web resources
[**CARD/RGI**](https://card.mcmaster.ca/analyze/rgi) can be used to predict
resistomes from protein or nucleotide data based on homology and SNP models.


``` bash
module unload prokka/1.14.6
module unload amrfinder/4.0.22
module load amrfinder/4.0.22
```

## AMR genes detection using AMRFinder

NCBI--[***AMRFinderPlus***](https://github.com/ncbi/amr/wiki/Home). It
identifies acquired antimicrobial resistance genes in bacterial protein and/or
assembled nucleotide sequences as well as known resistance-associated point
mutations for several taxa.

The most critical part of any AMR tool is the underlying database. So we take time to understand the databse on which AMRFinderPlus is based.

"This database is derived from the [Pathogen Detection Reference Gene Catalog](https://www.ncbi.nlm.nih.gov/pathogens/isolates#/refgene/), [Pathogen Detection Reference Gene Hierarchy](https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy/), and [Reference HMM Catalog](https://www.ncbi.nlm.nih.gov/pathogens/hmm/) and is used by the [Pathogen Detection](https://ncbi.nlm.nih.gov/pathogens/) isolate analysis system to provide results to the [Isolates browser](https://www.ncbi.nlm.nih.gov/pathogens/isolates) and [MicroBIGG-E](https://www.ncbi.nlm.nih.gov/pathogens/microbigge) as well as the command-line version of AMRFinderPlus. The 'core' subset version focuses on acquired or intrinsic AMR gene products including point mutations in a limited set of taxa. As of version 4.0 AMRFinderPlus also includes [StxTyper](https://github.com/ncbi/stxtyper) which has a separate DNA-sequence database and algorithm for typing Stx operons.

The 'plus' subset include a less-selective set of genes of interest including genes involved in virulence, biocide, heat, metal, and acid resistance.

>**Note:** that AMRFinderPlus reports gene and point mutation presence/absence; it does not infer phenotypic resistance. Many of the resistance genes detected by AMRFinderPlus may not be relevant for clinical management or antimicrobial surveillance, AMRFinderPlus.

Let's get started by defining a temporary directory for the tool and define DB
path:


``` bash
export TMPDIR=./results/illumina/klebs/tmp/amrfinder/

```

<!-- ```
AMRFINDER_DB=$(find ./databases/amrfinderplus/2023-11-15.1/ -name "AMR.LIB" | sed 's=AMR.LIB==')
``` -->


You can update the database before the analysis to ensure the latest
information using the command `amrfinder --update`

``` bash
amrfinder \
    --nucleotide ./results/illumina/klebs/prokka/SRR28370701.fna \
    --protein   ./results/illumina/klebs/prokka/SRR28370701.faa \
    --gff       ./results/illumina/klebs/prokka/SRR28370701.gff \
    --annotation_format prokka \
    --organism Klebsiella_pneumoniae \
    --plus \
    --ident_min -1 \
    --coverage_min 0.5 \
    --translation_table 11 \
    --database $AMRFINDER_DB \
    --threads 2 \
    --name SRR28370701 \
    > ./results/illumina/klebs/amrfinder/SRR28370701.tsv
```


## Output Format
AMRFinderPlus produces a tab-separated output file with detailed information:
Column | Description
---|---
Protein identifier | Protein/contig identifier in the input sequence
Gene symbol | Gene symbol for the AMR gene
Sequence name | Matching sequence name in the AMR database
Scope | Core, Plus (e.g., virulence factors)
Element type | AMR gene class (e.g., AMR, STRESS, VIRULENCE)
Element subtype | Specific resistance mechanism
Class | Antibiotic class
Subclass | Specific antibiotic subclass
Method | Detection method used (PARTIALP, EXACTP, BLASTX, HMM)
Target length | Reference sequence length
Reference sequence length | Length of matching reference
% Coverage of reference | Percentage of reference covered by alignment
% Identity to reference | Percent identity to reference sequence
Alignment length | Length of the alignment
Accession of closest sequence | Accession number of the matching sequence
Name of closest sequence | Name of the matching sequence
HMM id | Identifier of HMM used for detection (if applicable)
HMM description | Description of HMM (if applicable)




# Step 8: Variant Calling and Consensus Assemblies


Prepare the working environment

``` bash 
module purge
module load snippy/4.6.0
module load gubbins/3.4
```

We will now focus on using an alignment-based comparison approach to identify
relationships within pathogen isolates. To model an epidemilogical situation involving different pathogens in circulation, 11 isolate assemblies collected in Kenya between January 14 and January 31, 2019 were retrieved from [PathogenWatch](https://pathogen.watch/genomes/all?country=ke&genusId=570&maxDate=2019-01-31T20%3A59%3A59.999Z&minDate=2018-12-31T21%3A00%3A00.000Z&sort=date&speciesId=573)  as context to ourassembled isolate. They are located in `./data/klebs/pathogenwatch/genomes`.


#### Fast Bacterial Variant Calling with Contigs
Snippy finds SNPs between a haploid reference genome and your NGS sequence
reads. It will find both substitutions (snps) and insertions/deletions (indels).
It will use as many CPUs as you can give it on a single computer (tested to 64
cores). It is designed with speed in mind, and produces a consistent set of
output files in a single folder. It can then take a set of Snippy results using
the same reference and generate a core SNP alignment (and ultimately a
phylogenomic tree).

Define the temporary directory for `snippy` analysis.

``` bash 
export TMPDIR=$(pwd)/results/illumina/klebs/tmp/snippy/

mkdir -p ./results/illumina/klebs/tmp/snippy
```


Contig-based variant calling with `snippy`:


``` bash
for fn in ./pathogenwatch/klebs/assemblies-to-test/*.fasta; do
    sample=$(basename $fn)
    sample="${sample%.*}"
    echo -e "-------------------------------\n"
    echo -e "running snippy on: $sample - $fn"
    echo -e "-------------------------------\n"
    
    snippy \
        --force \
        --prefix $sample \
        --mapqual 60 \
        --basequal 13 \
        --cpus 2 \
        --ram 8 \
        --tmpdir $TMPDIR \
        --outdir ./results/illumina/klebs/snippy/$sample \
        --reference ./genomes/klebs/GCF_000016305.1_ASM1630v1_genomic.gbff \
        --ctgs $fn
    echo -e "\n-------------------------------\n"
done

```

`--mapqual` is the minimum mapping quality to accept in variant calling. 
BWA MEM using `60` to mean a read is "uniquely mapped".

`--basequal` is minimum quality a nucleotide needs to be used in variant
calling. We use `13` which corresponds to error probability of ~5%. It is a
traditional SAMtools value.

The variant calling is done by Freebayes. The key parameters under user control are:

`--mincov` - the minimum number of reads covering a site to be considered
(default=10)

`--minfrac` - the minimum proportion of those reads which must differ from the
reference

`--minqual` - the minimum VCF variant call "quality" (default=100)



##### `Snippy` Outputs
Snippy can perform variant calling based on reads or contigs. We used the latter approach which is associated with the following output files (Some of them are universal to either approches except the amount of infomation may differ).

File | Description
--- | ---
`*.tab` | Lists all variants in a human-readable form; favourable for quick epidemilogical insight
`*.vcf` | Universal variant calling format parsible by intergrating bioinformatics tools
`*.gff` | Provides genomic feature context for each variant
`*.txt` | Lists numbers of different variant types
`*.consensus.fa` | FASTA file of consensus sequence; variants are incorporated into the reference
`*.aligned.fa` | The alignment between the reference and the consensus in FASTA format; useful for phylogenetics
`*.bam` | BAM file for **contigs** alignment agains the reference
`*.html` | HTML summary of variants
`*.txt` | A simple summary of the run; lists numbers of different variant types


##### Visualising the Snippy Variants

Copy the following files from the output directory to your `home` directory on
`hpc`:

- Reference genome file `ref.fa` and Reference genome index file `ref.fa.fai`

``` bash
rsync \
    -avP \
    --partial \
    ./results/illumina/klebs/snippy/SRR28370701/reference/ref.fa \
    ./results/illumina/klebs/snippy/SRR28370701/reference/ref.fa.fai \
    ~/
```

- BAM file `SRR28370682.bam` and Indexed BAM file `SRR28370682.bam.bai`

```
rsync \
    -avP \
    --partial \
    ./results/illumina/klebs/snippy/SRR28370701/SRR28370701.bam* \
    ~/
```


> Note. On a new terminal session outside the `hpc`. Replace **<user_name>** with your
> `user` name

Copy the files from your `home` directory on `hpc` to your `home` directory on
`local` machine:

```
rsync -avP --partial <user_name>@hpc.ilri.cgiar.org:~/ref.fa* ~/
```

```
rsync -avP --partial <user_name>@hpc.ilri.cgiar.org:~/SRR28370682.bam* ~/
```


1. Go to https://igv.org/app/
2. Load the Genome `ref.fa` in the `Genome` tab
3. Load the alignment and its index in the `Tracks` tab
4. For the `ref.fa` select `NC_009648` as chromosome and type `NC_009648:1-200`
   as the region of interest


<br>
<left><img src="img/igv-screenshot.png" alt="Screenshot of IGV" width="1532"/></left>
<br>


##### Build Core and Whole Genome Aligments from Snippy Output

Snippy-core combines multiple Snippy outputs to:

- Identify the conserved-sequence genome ("core genome") which are regions present in all isolates.
- Extract SNPs from these core regions
- Create a multiple sequence alignment of core SNPs (SNP absence/presence-based avoiding the full compilation of variation such as `INS`, `DELS`, variant types)
- Generate various outputs suitable for downstream phylogenetic analysis

#### Run snippy-core
`snippy-core` is a tool that allows you to compile a core **SNP alignment** from
multiple Snippy analyses. It's particularly useful for phylogenetic analysis of closely related bacterial isolates.


If you want to mask certain regions of the genome, you can provide a BED file 
with the `--mask` parameter. Any SNPs in those regions will be excluded. 
This is common for genomes like M.tuberculosis where pesky repetitive 
PE/PPE/PGRS genes cause false positives, or masking phage regions. 
A `--mask` bed file for M.tb is provided with Snippy in the 
`etc/Mtb_NC_000962.3_mask.bed` folder. 
It is derived from the XLSX file from https://gph.niid.go.jp/tgs-tb/

```
snippy-core \
    --maxhap 100 \
    --mask-char N \
    --ref ./genomes/klebs/GCF_000016305.1_ASM1630v1_genomic.gbff \
    --prefix ./results/illumina/klebs/snippy-core/core-snp \
    ./results/illumina/klebs/snippy/*
```

#### Snippy Core Outputs
File | Description
--- | ---
`*.aln` | Alignment of core SNPs in FASTA format
`*.full.aln` | Alignment of entire core genome (including invariant sites)
`*.vcf` | Multi-sample VCF file with all core SNP sites
`*.tab` | Tab-separated table of core SNPs with annotations
`*.text` | Summary statistics text file
`*.ref.fa` | The reference sequence used
`*.aligned.fa` | Contains one sequence per isolate including the reference
<!-- `*.nwk` | Fast neighbor-joining tree based on core SNP alignment -->

#### Cleanup the Snippy SNP Alignment
The resultant `core-snp.full.aln` is loaded with alphabetical encoding (see [Snippy](https://github.com/tseemann/snippy/tree/master) for details) which may not be compatible with downstream analyses like tree-building or recombination-removal. The numerous encodings may be simplified  as below:

``` bash 
snippy-clean_full_aln \
    /results/illumina/klebs/snippy-core/core-snp.full.aln > \
    ./results/illumina/klebs/snippy-core/core-snp-clean.full.aln
```

# Step 9: Phylogenetic analysis

## Detect Recombination and mask recombinant regions

Gubbins (Genealogies Unbiased By recomBinations In Nucleotide Sequences) is a toolfor detecting and accounting for homologous recombination in bacterial whole genome alignments. Since homologous recombination can obscure the true relatedness/vertical inheritance, it may be desirable to find and mask such regions before phylogenetic inference.

Gubbins:
- Takes a core genome alignment (e.g., from Snippy-core) as input.
- Detects recombination regions using elevated SNP density or maximum likelihood methods (RAxML/FastTree/PhyML).
- Iteratively masks/"remove" recombination and infers a clonal phylogeny.

Predicting recombination regions:

We will submit the job using the job scheduler SLURM. 

1. Download the `run_gubbins.sh` to a scripts directory in your home directory

``` bash 
mkdir -p ~/scripts
wget -c -N https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/scripts/run_gubbins.sh \
-P ~/scripts/
```

2. Edit the downloaded script

>Note Ensure that you edit the partition option `-w` in the script accordingly to
>correspond to the compute node that you were assigned.


3. Submit it as a job

```
sbatch ~/scripts/run_gubbins.sh
```

>**Note: DO NOT RUN THE SCRIPT INTERACTIVELY**

Interactively, we'd have run the job as follows:

```
run_gubbins.py \
    --threads 2 \
    --prefix ./results/ont/klebsiella/gubbins/core-snp \
    --iterations 5 \
    --min-snps 3 \
    --min-window-size 100 \
    --max-window-size 10000 \
    --filter-percentage 25.0 \
    ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln
```

#### Output files
File | Description
--- | ---
`.recombination_predictions.embl` | Recombination predictions in EMBL file format.
`.recombination_predictions.gff` | Recombination predictions in GFF format
`.branch_base_reconstruction.embl` | Base substitution reconstruction in EMBL format
`.summary_of_snp_distribution.vcf` | VCF file summarising the distribution of point mutations
`.per_branch_statistics.csv` | per branch reporting of the base substitutions inside and outside recombination events
`.filtered_polymorphic_sites.fasta` | FASTA format alignment of filtered polymorphic sites used to generate the phylogeny in the final iteration
`.filtered_polymorphic_sites.phylip` | Phylip format alignment of filtered polymorphic sites used to generate the phylogeny in the final iteration
`.final_tree.tree` | final phylogeny in Newick format; branch lengths are in point mutations
`..node_labelled.final_tree.tre` | final phylogenetic tree in Newick format but with internal node labels; branch lengths are in point mutations


To generate a recombination-masked alignment (i.e., with sequences predicted to
have been introduced by recombination removed, leaving just the clonal frame),
the post-processing script `mask_gubbins_aln.py` can be used:


```
mask_gubbins_aln.py \
    --aln ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln \
    --gff ./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.gff \
    --missing-char N \
    --out ./results/ont/klebsiella/gubbins/core-snp.masked.aln
```


<!-- We can then convert the Gubbins `.gff` masking feature output into a `BED`
format file usable with `snippy` as argument to the `--mask` option.

Convert Gubbins GFF to BED format: 
```
module load bedops/2.4.29
```

```
gff2bed < ./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.gff > ./results/ont/klebsiella/gubbins/gubbins_recomb.bed
```

<details>
    <summary>Click to toggle a <b style='color:blue'>Challenge</b>
    </summary>

Now, can you use the resultant `BED` file to re-run Snippy all the way to building a phylogenrtic tree. How do the trees compare with(out) masking of recombinant regions?
</details> -->


To produce publication-ready figures of Gubbins analyses and annotate with
metadata, clades or any other information of the samples.

```
awk -F"\t" 'BEGIN{printf "id,ST\n"} {printf "%s,%s\n", $1,"ST"$3}' ./results/ont/klebsiella/mlst/klebs-mlst.txt | sed s'/.fasta//g' | sed s'/.*\///g' > ./results/ont/klebsiella/gubbins/core-snp-metadata.csv
```

```
Rscript ./scripts/plot_gubbins.R \
  -t ./results/ont/klebsiella/gubbins/core-snp.final_tree.tre \
  -r  ./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.gff \
  --meta ./results/ont/klebsiella/gubbins/core-snp-metadata.csv \
  --legend-height 0.35  \
  --tree-axis-expansion 30 \
  --heatmap-x-nudge 0.05 \
  --heatmap-y-nudge -0.05 \
  -o ./results/ont/klebsiella/gubbins/core-snp-plots.png
```

<br>
<left><img src="img/core-snp-plots.png" alt="predicted regions of recombination" width="1532"/></left>
<br>


<!-- To create a PDF showing	the	**predicted regions of recombination** against a	
phylogenetic reconstruction	based on the final iteration of the	Gubbins
analysis:

```
gubbins_drawer.py \
  –o core-snp-tree-recombination-regions.pdf \
  -t ./results/ont/klebsiella/gubbins/core-snp.final_tree.tre \
  ./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.embl
```


For each isolate, blocks representing the regions identified as recombinations
by gubbins are indicated by coloured blocks. Blue blocks are unique to a single 
isolate	while red blocks are shared	by multiple	isolates. The horizontal
position of	the	blocks represents their	position in the alignment.


To create a PDF showing	the	**reconstructed SNPs** against the same tree:

```
gubbins_drawer.py \
  –o core-snp-tree-SNPs.pdf \
  -t ./results/ont/klebsiella/gubbins/core-snp.final_tree.tre \
  ./results/ont/klebsiella/gubbins/core-snp.branch_base_reconstruction.embl
``` -->


## Maximum likelihood phylogenetic inference

We can use another phylogenetic tool `iqtree` to visualise the relationships
between the recombination-masked isolates from Gubbins.

```
module purge
module load iqtree/1.6.12
```

```
iqtree \
    -m HKY \
    -bb 1000 \
    -alrt 1000 \
    -alninfo \
    -s ./results/ont/klebsiella/gubbins/core-snp.masked.aln \
    -nt 2 \
    -redo \
    -pre ./results/ont/klebsiella/iqtree/core-snp
```

> **Note:** It is important to use phylogenetic algorithms that take into
> account SNP alignments. These algorithms usually include some form of
> ascertainment bias correction that corrects for the 'missing' nucleotides in
> the alignment that were masked/removed because they did not show polymorphism.


# Step 10: Visualization

## Visualize the phylogeny alongside typing, antibiotic resistance or epidemiological data


We will explore [**Microreact**](https://microreact.org/) to visualize genomic
data alongside AMR, typing (MLST, serotyping) and epidemiological data.


Microreact is software developed by the Centre for Genomic Pathogen Surveillance
(CGPS) that allows you to upload, visualise and explore any combination of
**clustering (trees)**, **geographic (map)** and **temporal (timeline)** data.
Other metadata variables are displayed in a table. You can specify colours
and/or shapes to display on the map, tree and/or timeline. A permanent URL is
produced for you to share your Microreact, or a .microreact file can be
downloaded for sharing with collaborators.


[**Microreact**](https://microreact.org/) is a web-based application for the
provision of interactive data visualisations. It enables the rapid generation
and linkage of trees, maps, networks, charts and timelines, enabling
epidemiologists and key decision makers to react faster and with greater
accuracy.


1. Microreact is freely available (at http://microreact.org).
2. Microreact – Hierarchical and Geographical Analysis Tool.
3. Microreact allows you to upload, visualize and explore dendrograms (trees)
   linked to adata containing geographical locations.


1. Download the data to your local computer through either the command line or
   the **Download button**

    ```
    wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/microreact/data.csv
    ```

    ```
    wget https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/microreact/tree.nwk
    ```

2. On your preferred web browser (Firefox, Google Chrome, Safari, Microsoft
   Edge), open a new window and type http://microreact.org on the address bar.

3. Login here: https://microreact.org/api/auth/signin?callbackUrl=/my-account



## Visualization using other tools
We can also use open-source tools to visualize the complex genomic data
integrated with other data types - AMR profiles, Serotying and MLST information.

1. Login to the Rstudio server https://hpc.ilri.cgiar.org/rstudio/
2. Start a new `Terminal` session on the `Terminal Tab`
3. Download the `tree.nwk` file into `home` directory

    ```
    wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/microreact/tree.nwk -P ~/
    ```
4. Download the `data.csv` file into `home` directory

    ```
    wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/microreact/data.csv -P ~/
    ```
5. Download the `plot_tree.R` script into home directory

    ```
    wget -c https://raw.githubusercontent.com/ILRI-Genomics-Platform/AMR-Genomic-Surveillance/refs/heads/main/scripts/plot_tree.R -P ~/
    ```
6. Navigate to the `Output` pane/panel (Output pane, containing the Files,
   Plots, Packages, Help, Viewer, and Presentation tabs) on Rstudio and click on
   the `Refresh icon`

7. Click on the `plot_tree.R` script to open it on the `Source` pane/panel.


# Assignment

## 1. Visualize the generated phylogenetic tree of the 11 samples alongside the
   metadata in Microreact. 

   - Create a new project and name it accordingly.
   - Create a metadata with the relevant information. You can be guided on how
     to prepare your metadata using the guidelines outlined here:
     https://docs.microreact.org/instructions/creating-a-microreact-project/metadata-column-types
     


## 2. Visualize the data using R

We will copy the `mlst`, `core-snp.treefile` and `resfinder` results from the 
compute nodes (`compute05` or `compute06`) to the head node (`hpc`) to our
`home` directory interactively as follows:

```
rsync -avP --partial \
    /var/scratch/$USER/ACDC_AMR2025/results/ont/klebsiella/mlst ~/
```

```
rsync -avP --partial /var/scratch/$USER/ACDC_AMR2025/results/ont/klebsiella/resfinder --exclude="*_blast*" ~/
```


```
rsync -avP --partial \
    /var/scratch/$USER/ACDC_AMR2025/results/ont/klebsiella/iqtree ~/
```



Symbolically link the script to the `scripts` directory in the `home` directory


```
ln -sf /var/scratch/global/jjuma/ACDC_AMR2025/scripts/visualizeAMR* ~/scripts/
```


>**Note: Run the script interactively on the `rstudio` server `Terminal`** 

```
module load R/4.4
```

```
Rscript ~/scripts/visualizeAMR.R \
    --tree ~/iqtree/core-snp.treefile \
    --mlst ~/mlst/*.tsv \
    --pointfinder ~/resfinder/*/PointFinder_results.txt \
    --resfinder ~/resfinder/*/ResFinder_results_tab.txt \
    --prefix phylogeny-amr \
    --outdir ~/plots
```