# Bioinformatics Analysis of Klebsiella pneumoniae of Illumina
# Group 1

---  

###### **Group_members **: [Betselot Zerihun Ayano ](https://github.com/betselotz),  [Vera Morangi Onwonga](https://github.com/Vmorangi), Kiconco Benadine, Nkiruka Lynda Uzoebo & Nambozo Eunice Jennifer

---
##  Introduction – *Klebsiella pneumoniae*

#### Opportunistic Pathogen  
- Gram-negative bacterium  
- Causes hospital-acquired infections: **pneumonia, bloodstream infections, urinary tract infections**  

#### Antimicrobial Resistance (AMR)  
- Produces **extended-spectrum β-lactamases (ESBLs)** and **carbapenemases** (e.g., KPC, NDM)  
- Severely limits treatment options, including last-resort **carbapenems**  

#### Drug-Resistant Strains  
- **MDR (multidrug-resistant)** & **XDR (extensively drug-resistant)**  
- Spread via **plasmid-mediated genes** and **successful clonal lineages**  

#### Genomic Insights  
- **Whole-genome sequencing (WGS):** tracks resistance genes & plasmids  
- **Phylogenomics:** reveals clonal expansion and transmission routes  
- **Pangenome analysis:** highlights core vs. accessory genes driving adaptation  
- **Mobile genetic elements (MGEs):** key role in resistance gene dissemination  

### Public Health Impact  
- Rapid global dissemination in **healthcare settings**  
- Major challenge for **infection prevention and treatment**  

##  Setting up the Bioinformatics Analysis Environment  

###  Log in to the HPC  
First, log in to the High-Performance Computing (HPC) cluster using your assigned credentials:  
```bash
ssh <user_name>@hpc.ilri.cgiar.org
```
We all worked on 

>Compute05  
``` bash
interactive -w compute05 -c 2 -J amr-surveillance -p batch
```

###  Setting Up the Project Directory Structure  

Create a working directory on the HPC scratch space and move into it:  

```bash
mkdir -p /var/scratch/$USER/illumina_downloads
```
We changed directory into 

```bash
cd /var/scratch/$USER/illumina_downloads
```


##  Bioinformatics Analysis  
Group leader created a `data` directory inside `illumina_downloads` to download our raw sequencing data:  

```bash
mkdir data
cd data
```
### 📥 Download Data from NCBI  (Group leader)

---

#### Using SRA Explorer  

1. Go to **[SRA Explorer](https://sra-explorer.info/#)**  
2. Search for:  SRR28370701

3. In the search results, check the following boxes:  
- ✅ *WGS of Klebsiella pneumoniae isolate*  
- ✅ *Add 1 to collection*  
- ✅ *1 data saved*  
- ✅ *Bash script for downloading FastQ files*  

4. Download the generated Bash script (e.g., `sra_download.sh`)  
5. Run the script to download the FASTQ file(s):  
```bash
#!/usr/bin/env bash     # tells system to run this with bash

# download read 1
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR283/001/SRR28370701/SRR28370701_1.fastq.gz \
  -o SRR28370701_WGS_of_Klebsiella_pneumoniae_isolate_1.fastq.gz

# download read 2
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR283/001/SRR28370701/SRR28370701_2.fastq.gz \
  -o SRR28370701_WGS_of_Klebsiella_pneumoniae_isolate_2.fastq.gz
```
### 📝 Renaming FASTQ Files  

When downloading from SRA Explorer, the FASTQ files may have long descriptive names such as:  


### Loop through every FASTQ file ending with .fastq.gz in the current directory (Group leader)
``` bash
for f in *.fastq.gz; do

    # Extract the accession ID (SRR followed by digits) from the filename
    # Example: "SRR28370701_WGS_of_Klebsiella_pneumoniae_isolate_2.fastq.gz" → "SRR28370701"
    base=$(echo "$f" | grep -oP 'SRR[0-9]+')

    # Extract the read number (_1 or _2) right before ".fastq.gz"
    # Example: "_2" from "..._isolate_2.fastq.gz"
    readnum=$(echo "$f" | grep -oP '_[12](?=\.fastq\.gz)')

    # Rename the file to "<SRR_accession><readnum>.fastq.gz"
    # Example: SRR28370701_WGS_of_Klebsiella_pneumoniae_isolate_2.fastq.gz → SRR28370701_2.fastq.gz
    mv "$f" "${base}${readnum}.fastq.gz"
done
```

### Create a Structured Project Directory  

Within our project directory, we can create well-organized subfolders for the *Klebsiella* analysis:  

```bash
mkdir -p \
results/illumina/klebs/{fastqc,fastp,fastq-scan,shovill,prokka,amrfinder,mlst,kaptive,resfinder,iqtree,snippy,snippy-core,gubbins,tmp/{shovill,prokka,resfinder,amrfinder,snippy}}
```

Group members created symbolic links to the relevant `data` folders in our current working directory.  

```bash
ln -sf /var/scratch/user3/illumina_downloads/[dps]* .
```

###  Load Required Modules  

Before running the analyses, load the necessary bioinformatics tools on the HPC:  

```bash
module load fastqc/0.11.9       # Quality control of raw reads
module load fastp/0.22.0        # Read trimming and filtering
module load bbmap/38.95         # Read mapping and QC utilities
module load shovill/1.1.0       # Genome assembly (Illumina)
module load prokka/1.14.6       # Genome annotation
module load resfinder/4.6.0     # Antimicrobial resistance gene detection
module load mlst/2.23.0         # MLST typing
module load amrfinder/4.0.22    # AMR gene detection and annotation
```


### Assessing Read Quality using FastQC (before quality trimming)

We used **FastQC** to check the quality of raw sequencing reads before trimming.  

```bash
fastqc \
    -o ./results/illumina/klebs/fastqc \
    --noextract \
    -t 2 \
    --nogroup \
    ./data/klebs/SRR28370701_1.fastq.gz \
    ./data/klebs/SRR28370701_2.fastq.gz
```

## Download and Explore FastQC Results

After running `FastQC`, we downloaded the resulting HTML reports to our HPC home directory and our local computer for visualization.  

1. **Copy the HTML results to your home directory on the cluster**  

```bash
cp ./results/illumina/klebs/fastqc/*html ~/
```
2. **Transfer the Report to Your Local Computer**  

We use `rsync` to securely copy the file:  

```bash
rsync -avP --partial <USERXX>@hpc.ilri.cgiar.org:~/SRR28370701*.html ~/AMR_training/group1/     
```

### Quality Trimming FASTQ Files with fastp (Adapter and Quality Trimming)

We used **fastp** to trim adapters, filter low-quality bases, and clean up paired-end reads before assembly.  

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
After running `fastp`, we downloaded the resulting HTML reports to our HPC home directory and our local computer for visualization.  

1. **Copy the HTML results to your home directory on the cluster**

```bash
cp ./results/illumina/klebs/fastp/SRR28370701.fastp.html ~/
```
2. **Transfer the Report to Your Local Computer** 

```bash
rsync -avP --partial <USERXX>@hpc.ilri.cgiar.org:~/SRR28370701.fastp.html ~/AMR_training/group1/ 
```
 - Examine the report.
 
### De novo Genome Assembly: Shovill

**Shovill** is a pipeline designed for **fast and automated assembly of bacterial genomes** from Illumina paired-end reads.  

**Key points:**  
- Uses **SPAdes, SKESA, or other assemblers** internally to build contigs.  
- Performs **adapter trimming, error correction, and read filtering** before assembly.  
- Produces a **high-quality draft genome** suitable for downstream analyses such as annotation, MLST, and AMR gene detection.  
- Optimized for **speed and low memory usage** on HPC environments.  

**Why we use it:**  
- Converts raw sequencing reads into **contiguous genome sequences** (contigs).  
- Provides a foundation for **genome annotation, typing, and comparative genomics**.

```bash
# Set temporary directory for Shovill intermediate files
export TMPDIR="./results/illumina/klebs/tmp/shovill"
```
###  Run Shovill assembly
```bash
shovill \
  --R1 ./results/illumina/klebs/fastp/SRR28370701_trim_R1.fastq.gz \
  --R2 ./results/illumina/klebs/fastp/SRR28370701_trim_R2.fastq.gz \
  --gsize 5700000 \
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

### Rename Assembly Output

After running Shovill, the main contigs file is named `contigs.fa` by default.  
We rename it to include the **sample accession** for clarity:

```bash
mv ./results/illumina/klebs/shovill/SRR28370701/contigs.fa \
   ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa
```

###  Assembly Evaluation

Before downstream analyses, it is important to verify the quality of the assembled genome.

#### 1. Quick assembly stats using `stats.sh`

```bash
stats.sh in=./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa
```

#### 2. Using seqkit to explore assembly
##### List available seqkit modules
``` bash
module avail seqkit
```
###### Load seqkit
``` bash
module load seqkit/0.11.0 
```
###### Display help
``` bash
seqkit -h
```
###### Convert FASTA to tab-delimited table (sequence length and name)
``` bash
seqkit fx2tab -nl ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa
```
We can use another tool assembly-scan to generate summary statistics of the assembly.

#### 3. Assembly summary with assembly-scan
######  List available modules
``` bash
module avail assembly-scan
```
######  Load assembly-scan
``` bash
module load assembly-scan/1.0.0
```
######  Generate summary statistics
``` bash
assembly-scan ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa \
  --transpose \
  | tee ./results/illumina/klebs/shovill/SRR28370701/SRR28370701-assembly-scan.tsv
```
##### 4. Compute GC content from assembly-scan output
``` bash
grep 'contig_percent_[cg]' \
  ./results/illumina/klebs/shovill/SRR28370701/SRR28370701-assembly-scan.tsv \
  | awk -F '\t' '{sum+=$3} END {print "GC%=",sum}'
```


### 💻 Transfer Assembly to Local Computer

#### 1. Copy assembly to HPC home directory
```bash
rsync -avP --partial \
    ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa \
    ~/
```
```bash                             
#### 2. Download assembly to local machine

rsync -avP --partial <USERXX>@hpc.ilri.cgiar.org:~/SRR28370682.fa \
    ~/AMR_training/group1/
```

## Step 5: Genome Annotation with Prokka

**Prokka** is a widely used tool for rapid **annotation of bacterial genomes**. It predicts coding sequences (CDS), rRNAs, tRNAs, and other features, and generates standardized output files like GFF, GBK, FAA, and FNA. Using Prokka ensures annotations are reproducible and compatible with downstream bioinformatics tools.

---

### 1. Set Temporary Directory for Prokka
```bash
export TMPDIR=./results/illumina/klebs/tmp/prokka/
```

### 2.  Run Prokka for Genome Annotation
```bash
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
##  Add Additional Genomes from Pathogenwatch

To better understand the genomic diversity and epidemiology of Klebsiella in Kenya, we will include **11 additional Klebsiella isolates** collected between January 14 and January 31, 2019.  
This allows us to:  

- Compare our locally sequenced isolate (`SRR28370701`) to other regional isolates.  
- Place our genome in a broader **phylogenetic and genomic context**.  
- Evaluate the presence and spread of antimicrobial resistance genes across multiple isolates.  

Data source: [Pathogenwatch Klebsiella genomes](https://pathogen.watch/genomes/all?country=ke&genusId=570&maxDate=2019-01-31T20%3A59%3A59.999Z&minDate=2018-12-31T21%3A00%3A00.000Z&sort=date&speciesId=573)

---

#### 1. Create a directory for global assemblies
```bash
mkdir -p ./pathogenwatch/klebs/assemblies-to-test   # -p creates parent directories if they do not exist
```
#### 2.  Add our genome assembly to the directory
```bash
rsync -avP --partial \
    ./results/illumina/klebs/shovill/SRR28370701/SRR28370701.fa \
    ./pathogenwatch/klebs/assemblies-to-test/SRR28370701.fasta
```
#### 3. Add reference genome assembly
```bash
rsync -avP --partial \
    ./genomes/klebs/GCF_000016305.1_ASM1630v1_genomic.fna \
    ./pathogenwatch/klebs/assemblies-to-test/Reference.fasta
```
#### 4. Add selected global assemblies to test
```bash
rsync -avP --partial \
    ./genomes/klebs/*.fna \
    ./pathogenwatch/klebs/assemblies-to-test/
```

##  Step 6: Pathogen Typing (MLST)

### Introduction to Multi-Locus Sequence Typing (MLST)

**Multi-Locus Sequence Typing (MLST)** is a widely used method for **characterizing bacterial isolates at the strain level** based on the sequences of several (typically 5–7) conserved housekeeping genes. MLST provides a standardized, reproducible, and portable way to assign **sequence types (STs)** to bacterial strains, allowing comparisons across studies and global databases.

#### Key Features:
- **Purpose:** Identify clonal lineages, track outbreaks, and study population structure.
- **Resolution:** Intermediate – higher than species-level typing, lower than whole-genome SNP analysis.
- **Species:** MLST schemes are species-specific; for example, Klebsiella pneumoniae has its own curated MLST scheme in pubMLST.

#### Input:
- **Annotated genome sequences** (FASTA format) or assembled contigs.
- Can also use raw reads, but most pipelines prefer **assembled genomes** to ensure accuracy.

#### Output:
- **Allelic profile table**: Lists the alleles present at each housekeeping gene locus.
- **Sequence Type (ST)**: A unique identifier corresponding to the combination of alleles.
- **Optional details**: BLAST scores, percent identity, coverage, and unassigned loci.
- **Format**: TSV, CSV, or plain text files suitable for downstream analysis or integration into global MLST databases.

MLST is widely used in **epidemiology, surveillance, and population genomics**, providing a consistent method to **compare bacterial isolates worldwide**.


#### 1. Set the MLST database path
```bash
MLST_DB=$(find /export/apps/mlst/2.23.0/db -name "mlst.fa" | sed 's=blast/mlst.fa==')
```

#### 2. Run MLST typing
```bash
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

### Batch MLST Typing

We can perform MLST on multiple assemblies in a single step using a loop. This is useful when comparing our isolate to additional genomes, such as those from Pathogenwatch.
```bash
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
```

#### Merge MLST Reports

After running MLST on multiple assemblies, we combined all individual TSV files into a single consolidated file for easier comparison and downstream analysis.

```bash
cat \
    ./results/illumina/klebs/mlst/*.tsv > \
    ./results/illumina/klebs/mlst/klebs-mlst.txt
```

### Serotyping using Kaptive

Kaptive is a bioinformatics tool used for **in silico identification of Klebsiella capsule (K) and O antigen loci** from genome assemblies.  

**Why we do this:**  
- **Capsule (K) and O antigen typing** is critical because these surface polysaccharides are major virulence factors in Klebsiella.  
- Helps identify **high-risk clones** associated with outbreaks or severe infections.  
- Enables **comparative genomics** and tracking of **capsule and LPS diversity** across isolates.  
- Supports **vaccine research and epidemiology**, since K and O antigens are targets for vaccine development.  
- Complements MLST typing by providing **additional discriminatory power** for strain characterization.

---

#### K antigen and locus
-  The **K antigen** is the **capsular polysaccharide** surrounding Klebsiella cells.  
- It plays a major role in **immune evasion and virulence** by protecting bacteria from phagocytosis and complement-mediated killing.  
- The **K locus** in the genome encodes enzymes responsible for synthesizing the capsular polysaccharide.  
- Kaptive uses the K locus to **predict the capsule type** from genome assemblies. 

```bash
mkdir -p ./results/illumina/klebs/kaptive
kaptive assembly \
  /export/apps/kaptive/3.1.0/lib/python3.10/site-packages/reference_database/Klebsiella_k_locus_primary_reference.gbk \
  ./pathogenwatch/klebs/assemblies-to-test/*.fasta \
  --min-cov 70 \
  -t 2 \
  -o ./results/illumina/klebs/kaptive/kaptive_k_locus.tsv
```
#### O antigen and locus
- The **O antigen**  is part of the **lipopolysaccharide (LPS)**  on the outer membrane of Gram-negative bacteria.  
- It is important for **host immune recognition** and contributes to bacterial pathogenicity.  
- The **O locus**  contains genes responsible for synthesizing the O antigen sugar structures.  
- Kaptive identifies the **O antigen type** by comparing the O locus in genome assemblies to reference databases.  

```bash
kaptive assembly \
  /export/apps/kaptive/3.1.0/lib/python3.10/site-packages/reference_database/Klebsiella_o_locus_primary_reference.gbk \
  ./pathogenwatch/klebs/assemblies-to-test/*.fasta \
  --min-cov 70 \
  -t 2 \
  -o ./results/illumina/klebs/kaptive/kaptive_o_locus.tsv
```


### Kaptive Web
we uploaded our assemblies to [Kaptive-Web](https://kaptive-web.erc.monash.edu/). 

# Step 7: AMR genes detection


## AMR genes detection using ResFinder

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

``` bash
module unload prokka/1.14.6
module unload amrfinder/4.0.22
module load amrfinder/4.0.22
```

## AMR genes detection using AMRFinder

``` bash
export TMPDIR=./results/illumina/klebs/tmp/amrfinder/

```

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

# Step 8: Variant Calling and Consensus Assemblies


Prepare the working environment

``` bash 
module purge
module load snippy/4.6.0
module load gubbins/3.4
```


#### Fast Bacterial Variant Calling with Contigs

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



##### Build Core and Whole Genome Aligments from Snippy Output

#### Snippy Core Alignments

``` bash 
snippy-core \
    --maxhap 100 \
    --mask-char N \
    --ref ./genomes/klebs/GCF_000016305.1_ASM1630v1_genomic.gbff \
    --prefix ./results/illumina/klebs/snippy-core/core-snp \
    ./results/illumina/klebs/snippy/*
```

#### Cleanup the Snippy SNP Alignment

``` bash 
snippy-clean_full_aln \
    /results/illumina/klebs/snippy-core/core-snp.full.aln > \
    ./results/illumina/klebs/snippy-core/core-snp-clean.full.aln
```

# Step 9: Phylogenetic analysis

## Detect Recombination and mask recombinant regions

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

``` bash 
sbatch ~/scripts/run_gubbins.sh
```

>**Note: DO NOT RUN THE SCRIPT INTERACTIVELY**

Interactively, we'd have run the job as follows:

``` bash 
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


``` bash 
mask_gubbins_aln.py \
    --aln ./results/ont/klebsiella/snippy-core/core-snp-clean.full.aln \
    --gff ./results/ont/klebsiella/gubbins/core-snp.recombination_predictions.gff \
    --missing-char N \
    --out ./results/ont/klebsiella/gubbins/core-snp.masked.aln
```



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


## Maximum likelihood phylogenetic inference

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
