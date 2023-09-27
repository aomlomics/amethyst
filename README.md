# amethyst

AOML Metagenomics Workflow

## 0. Conda environments

A separate environment is used for each section of the workflow. Environments can be user-created or shared (e.g., `conda activate /home/cns.local/luke.thompson/miniconda3/envs/mg-qc`).

### 0.1. Sequence QC and trimming

* fastqc
* multiqc
* fastp
* seqkit

```
mamba create -n mg-qc -c bioconda fastqc multiqc fastp seqkit
```

### 0.2. Coverage and normalization

* nonpareil
* bbtools
* seqkit
```
mamba create -n mg-coverage -c bioconda -c agbiome nonpareil bbtools seqkit
```

### 0.3. Taxonomic and functional diversity

* mash
* sourmash
* krona
* humann
* metaphlan
* seqkit

```
mamba create -n mg-diversity -c bioconda -c biobakery sourmash krona humann metaphlan seqkit
```

### 0.4. Assembly and annotation

* spades
* prodigal
* prokka
* kofamscan
* seqkit

```
mamba create -n mg-assembly -c bioconda spades prodigal prokka kofamscan seqkit
```

### 0.5. Binning, mapping, QC, dereplication, and taxonomy

* bowtie2
* minimap2
* maxbin2
* metabat2
* checkm
* drep
* gtdbtk
* seqkit

```
mamba create -n mg-binning -c bioconda bowtie2 minimap2 maxbin2 metabat2 drep gtdbtk seqkit
mamba create -n mg-checkm -c auto checkm
```

## 1. Read quality control

Programs for inspecting FASTQ files before and after trimming and QC:

* [seqkit](https://bioinf.shenwei.me/seqkit/) -- general purpose fasta and fastq manipulations
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info/)

Example FastQC and MultiQC commands:

```
conda activate mg-qc
cd /path/to/fastq
mkdir fastqc-R1
mkdir fastqc-R2
fastqc *_R1_001.fastq.gz -o fastqc-R1
fastqc *_R2_001.fastq.gz -o fastqc-R2
cd fastqc-R1
multiqc --export .
cd ..
cd fastqc-R2
multiqc --export .
```

MultiQC is a great tool to visualize the quality and content of your reads. You can find a good tutorial on it [here](https://www.youtube.com/watch?v=qPbIlO_KWN0).

### Recommended program for adapter trimming and QC

* [fastp](https://github.com/OpenGene/fastp) 
    - adapter detection, quality score sliding window, very fast
    - occasionally misses adapters

```
mkdir QC
fastp \
-i 01.RAW/${n}_R1.fastq.gz \
-o QC/${n}_1.fq.gz \
-I 01.RAW/${n}_R2.fastq.gz \
-O QC/${n}_2.fq.gz \
--detect_adapter_for_pe \
-g -l 50 -W 4 -M 20 -w 16 \
--cut_front \
-R ${n}_report
```

TEST METAGENOME: You can explore the different outputs based on different settings for the following parameters: -W (sliding window size), -M (mean quality score in sliding window), and --cut_{right/front}. Three different outputs exist in the 01.fastp directory.

### Alternative programs for adapter trimming and QC

* FaQCs
    - adapter detection, per-base quality scores
* Cutadapt

```
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG
T -q 20 -o sample1_trimmed_F.fastq.gz -p sample1_trimmed_R.fastq.gz sample1_F.fastq.gz sample1_R.fastq.gz
```

* Trimmomatic 
* [Atropos](https://github.com/jdidion/atropos)

Example Atropos command:

```
atropos trim \
--front file:fwd_primer_permutations.fasta \
--front2 file:rev_primer_permutations.fasta \
--adapter file:fwd_primer_permutations_revcomp.fasta \
--adapter2 file: rev_primer_permutations_revcomp.fasta \
--times 2 \
--max-n 0 \
--input1 sample1_F.fastq.gz \
--input2 sample1_R.fastq.gz \
--output sample1_F_trim.fastq.gz \
--paired-output sample1_R_trim.fastq.gz \
--too-short-output sample1_F_tooshort.fastq.gz \
--too-short-paired-output sample1_R_tooshort.fastq.gz \
>/sample1.report
```

* Falco
* [Multitrim](https://github.com/KGerhardt/multitrim) (combo of FaQCs and fastp)

Example multitrim command:

```
multitrim -1 sample1_F.fastq.gz -2 sample1_R.fastq.gz -o samples/ --max
```

Multitrim command for multiple samples:

```
for f in 01.RAW/*R1.fastq.gz; do
	n=`echo ${f##*/} | cut -d _ -f 1-2` &&
	multitrim -1 $f -2 $n'_R2.fastq.gz' -o 02.QC/multitrim/ --max 
done
```

### Generate stats for pre- and post-trimming files using seqkit

```
seqkit stats *.fq.gz -T -j 20 > seqkit_fastq_stats.txt
```

Post-trimmed files should have higher average quality and contain no super-long/super-short reads.

## 2. Read-based analyses

Some questions are better answered using unassembled metagenomic reads rather than assembled contigs. Many reads will fail to assemble and that information is lost downstream.

### 2.1. Coverage analysis

Metagenome coverage using [Nonpareil](https://nonpareil.readthedocs.io/en/latest/).

Example Nonpareil command:
* Use only one of the read pairs (eg forward reads)

```
nonpareil -s reads_1.fa -b out_prefix -d 0.7 -t 20 -R 9500 -T alignment
```

Nonpareil output (.npo file) can then be loaded into R and a sequencing depth curve generated with the Nonpareil library.

```
install.packages('Nonpareil')
library(Nonpareil)
Nonpareil.curve('sample1.npo')
```

![](https://i.imgur.com/p6csMid.png)

The slope of the curve represents the sequence coverage (steeper slope = better coverage). In this example the blue samples on average have higher coverage with less sequencing effort than the red samples.

### 2.2. Normalization

Metagenome normalization (for better assembly)	

* BBNorm

Example BBnorm command:

* This command is for separate forward and reverse reads (output will be interleaved)
* Can also run using interleaved reads
```
bbnorm.sh in=sample1_F.fq.gz in2=sample1_R.fq.gz out=sample1-bbnorm.fq.gz \
target=100 min=5 interleaved=FALSE
 ```

### 2.3. Diversity analysis

Taxonomic composition:
- Sourmash, Kraken2, phyloFlash, singleM, kaiju
- BLAST against reference DB followed by MEGAN on BLAST results
- k-mer composition analyses using [Mash](https://mash.readthedocs.io/en/latest/) (MinHash method) or [Simka](https://github.com/GATB/simka)

### Recommended programs for diversity analysis

#### Mash

To generate a distance matrix of all your metagenomes based on k-mer composition, which correlates closely with taxonomic composition, use Mash. First create k-mer sketches of your unassembled reads ('mash sketch'), then combine all sketches ('mash paste') and generate the distance matrix (all-vs-all comparison with 'mash dist'). The example below uses k=25 (-k 25), reads instead of contigs (-r), a Bloom filter to filter out most unique k-mers with constant memory (-b 500M) and a threshold for stopping the calculation once a threshold of estimated average coverage is reached (-c 50). These settings may need to change based on the size of your metagenomes.

```
mash sketch -r -b 500M -c 50 -k 25 -o ${output_name} ${sample}

mash paste zr2760_all *.msh 

mash dist zr2760_all.msh zr2760_all.msh -t > zr2760_all_dist_matrix.txt
```

In a loop:

```
for f in *.fa; do
	n=`echo $f | cut -d . -f 1` &&
	mash sketch -r -b 100M -c 20 -k 25 -o $n $f
done
```

The resulting distance matrix can be used to generate ordination plots like NMDS to get a sense of metagenome similarities.

#### Sourmash

Similar to Mash, sourmash computes has sketches from DNA sequences and uses those sketches to compare genomes/metagenomes with each other or with query genomes/metagenomes. It also uses a k-mer based approach for genomic/metagenomic taxonomic composition with the NCBI and GTDB taxonomies. Documentation is [here](https://sourmash.readthedocs.io/en/latest/index.html).

Taxonomic assignment relies on pre-prepared databases. There are databases for GTDB (reference genomes only or all genomes) as well as GenBank genomes. You can download the pre-prepared databases [here](https://sourmash.readthedocs.io/en/latest/databases.html). We recommend k31-formatted .zip databases. For GTDB, you will also need to download the accompanying taxonomy sheet, which you can find at the species and strain levesls at the link above.

First interleave all paired reads with bbtools:

```
for f in *1.fastq.gz; do 
	n=`echo $f | rev | cut -d _ -f 2- | rev` && 
	reformat.sh in1=$f in2=$n'_2.fastq.gz' out=$n'_interleaved.fq.gz'
done
```

Then run in sourmash:

```
sourmash sketch dna {sample-interleaved}.fq.gz
```

in a loop
```
for f in *.fq.gz; do 
    sourmash sketch dna $f
done
```
run the signature against GTDB k31 pre-prepared database using 'gather'
```
sourmash gather {sample}.sig /phodnet/OCED/databases/sourmash/GTDB_k31.sbt -o {sample}_sourmash_gather_out.csv
```
run the gather output against the GTDB taxonomy sheet to generate (a) krona plot-formatted output and (b) taxonomic summary csv

```
sourmash tax metagenome -g {sample}_sourmash_gather_out.csv -t gtdb-rs202.taxonomy.v2.with-strain.csv -o {sample}_sourmash_tax --output-dir {sample}_sourmash_tax_out --output-format krona --rank species

sourmash tax metagenome -g {sample}_sourmash_gather_out.csv -t gtdb-rs202.taxonomy.v2.with-strain.csv -o {sample}_sourmash_tax --output-dir {sample}_sourmash_tax_out --output-format csv_summary
```

### Alternative programs for diversity analysis

#### Kraken2

Kraken2 matches each k-mer within a query sequence to the lowest common ancestor (LCA) of all genomes containing the given k-mer. The k-mer assignments inform the classification algorithm. Bracken allows users to estimate relative abundances within a specific sample from Kraken2 classification results. Bracken uses a Bayesian model to estimate abundance at any standard taxonomy level, including species/genus-level abundance.

Kraken2/Bracken databases can be downloaded from https://benlangmead.github.io/aws-indexes/k2.

Example Kraken2 and Bracken (run at species level) commands (where `${i}` is a variable representing the root of the sample name):

```
kraken2 \
--db /path/to/kraken2/database \
--threads 8 \
--paired \
--use-names \
--gzip-compressed \
--output /path/to/kraken2/output/${i}.kraken \
--report /path/to/kraken2/output/${i}.report \
/path/to/input/${i}_R1_trim.fastq.gz \
/path/to/input/${i}_R2_trim.fastq.gz

bracken \
-d /path/to/kraken2/database \
-l S \
-i /path/to/kraken2/output/${i}.report \
-o /path/to/kraken2/output/bracken-species/${i}.brackenreport
```

#### Visualization with Krona plots

In a Krona conda environment, run the following command to convert krona-formatted text files to an html visualization file:
```
ktimporttext {sample}-krona.tsv -o {sample}.krona.html
```
![](https://i.imgur.com/YtDL0ZM.jpg)

### 2.4. Functional analysis

[HUMAnN3](https://github.com/biobakery/humann) calculates community functional profiles stratified by known and unclassified organisms. It uses the UniRef database to provide gene family definitions, MetaCyc to provide pathway definitions by gene family, and MinPath to identify the set of minimum pathways. Accelerated mapping of reads to reference databases (including run-time generated databases tailored to the input) is performed by Bowtie2 for  nucleotide-level searches and Diamond for translated searches.

Instructions: https://huttenhower.sph.harvard.edu/humann

Generating a joint taxonomic profile: https://github.com/biobakery/humann#joint-taxonomic-profile
```
humann --input sample.interleaved.fastq --output humann_out --taxonomic-profile max.fixed.txt --threads 24
```

Guide to paired-end reads: https://github.com/biobakery/humann#markdown-header-humann-and-paired-end-sequencing-data
>HUMAnN3 does not take end-pairing relationships into account. As a result, the best way to use paired-end sequencing data with HUMAnN 3.0 is simply to concatenate all reads into a single FASTA or FASTQ file.

Joining and normalizing HUMAnN3's output tables: https://github.com/biobakery/humann

```
humann_join_tables --input genefamilies --output joined_genefamilies.tsv --file_name genefamilies

humann_renorm_table --input joined_genefamilies.tsv --output joined_genefamilies.relab.tsv --units relab
```

## 3. Assembly-based analyses

### 3.1. Assembly

Find overlapping short reads to make longer ones. Most tools use de Bruijn graphs to assemble contigs.

Popular assemblers:

* SPAdes (preferred; if slow see MEGAHIT below)

```
spades.py -pe -1 sample1_R1_trim.fq.gz -2 sample1_R2_trim.fq.gz -o output_dir 
```

* MEGAHIT (faster for larger datasets)

```
for f in 02.QC/*R1_trim.fq.gz; do
    n=`echo ${f##*/} | cut -d _ -f 1-2` &&
    megahit -1 $f -2 $n'_R2_trim.fq.gz --min-contig-len $MIN_CONTIG_SIZE -m 0.85 -o 03.ASSEMBLY/ -t $NUM_THREADS
done
```

* IDBA-UD (convert fastq to fasta before assembly)

```
idba_ud -i unassembled_input.fasta -o assembled_contigs.fasta --num_threads 24
```

* MetaVelvet
* OMEGA
* BIGMAC

### 3.2. Gene prediction

Identify likely protein-coding sequences (CDS) from assembled contigs. Usually done by extracting Open Reading Frames (ORFs), defined by start and stop codons.

ORF prediction programs:

* Prodigal (most common)

```
prodigal -i contigs.fa -o contig_cords.gbk -a contig_orfs.faa -d contig_orfs.fna
```

* FragGeneScan
* ORF Finder
* geneRFinder (machine-learning based)

At a minimum, you want nucleotide (.fna) and amino acid (.faa) outputs to search against different databases

### 3.3. Gene annotation

Functional annotation is the least standardized part of this work flow, with many different approaches and programs available for assigning biochemical function to genes. Here, we provide instructions for two options: (1) [Prokka](https://github.com/tseemann/prokka/blob/master/README.md#installation) and (2) KofamScan. By default, Prokka uses the manually curated SwissProt database, which is less comprehensive than databases like UniProt or GenBank but the annotations will have higher confidence. KofamScan uses the KEGG database, which is also curated. We recommend both tools for a comprehensive assessment of metagenomic functional potential.

In addition, there are tools for estimating metabolic pathways that take advantage of multiple databases and generate graphical outputs.  Three examples of these are: (1) [DRAM](https://github.com/WrightonLabCSU/DRAM) (Distilled and Refined Annotation of Metabolism, (2) [MEBS](https://github.com/valdeanda/mebs) (Multigenomic Entropy Based Score), and (3) [MicrobeAnnotator](https://github.com/cruizperez/MicrobeAnnotator). We will not discuss these pipelines but recommend checking them out once you're comfortable with basic annotation. 


#### 3.3.a. Prokka is a wrapper that incorporates several different tools, including RNAmmer, BLAST+, and Prodigal. 
```
prokka {contigs.fa} --outdir {prokka_annot_directory} --prefix {sample_name}
```

#### 3.3.b KofamScan is the downloadable software version of [KofamKOALA](https://www.genome.jp/tools/kofamkoala/), which uses Hidden Markov Models to search a KEGG Ortholog database of genes. The input file should be an amino acid sequence file of ORFs (e.g. from a Prodigal output). Both the program and databases are included in the conda installation (see section 0.4 for instructions).

```
ruby exec_annotation -o {outfile_name} {contig_orfs.faa}
```

### 3.4. Coverage of assembly (contigs), genes, and genomes

Calculate representation of sequences in the community.
Aligning (“mapping”) short reads onto assembled contigs and genes can provide relative abundance information. Most common identity threshold is 95%. Alignment files are usually in .bam or .sam format.

Common mappers:

* Bowtie2

```
bowtie2-build --threads 32 {sequences}.fa {db_name}
bowtie2 --threads 32 -x {db_name} -1 {sample1_trimmed_R1}.fq \
-2 {sample1_trimmed_R2}.fq -S {sample1}.sam > {sample1}.bowtie2.log
```

* Minimap2 (preferred; faster)

```
minimap2 -d catalog.mmi {contigs.}fna.gz # make index
minimap2 -t {num_threads} -N {50} -ax sr catalog.mmi \
{sample1_F_trimmed}.fq.gz {sample1_R_trimmed}.fq.gz | \
samtools view -F 3584 -b --threads {num_threads} > {sample1}.bam
```

* MagicBLAST (output is BLAST format rather than SAM or BAM files)

For a more thorough tutorial on how to use MagicBLAST for coverage calculation, see Roth Conrad's GitHub repo: https://github.com/rotheconrad/00_in-situ_GeneCoverage

```
makeblastdb -in sequences.fasta -dbtype nucl -parse_seqids -out db_name

magicblast -query {sample1}_concatenated_reads.fq \
-db {db_name} -infmt fastq -outfmt tabular -no_unaligned -splice F \
-out {out_table}.tsv -parse_deflines T -num_threads {num_threads}
```

### 3.5. Genome binning

Use k-mer composition and coverage information to group contigs into likely genomes (metagenome-assembled genomes, or MAGs). Different binning programs will produce different results and we recommend using multiple programs on your metagenomic assemblies, then using the highest-quality dereplicated MAGs from all programs.

Common binning programs:

* MaxBin2

```
run_MaxBin.pl -contig {contigs}.fa -min_contig_length 1000 \
-reads {sample1_R1_trimmed}.fq -reads2 {sample1_R2_trimmed}.fq \
-out {sample1}_MaxBin -thread {num_threads}
```

* MetaBat2

```
metabat2 -i {assembled_contigs}.fa -a {depth.txt} -o {output_dir} \
-t {num_threads} -m 1500
```

* GroopM
* Concoct
* BinSanity
* [Vamb](https://github.com/RasmussenLab/vamb) (machine-learning based)

### 3.6. Genome quality 

Completion and contamination/redundancy values are calculated based on the presence of single-copy core genes.

* CheckM

```
checkm lineage_wf -t 8 -x fa --tmp tmp_folder input_folder output_folder >output.log
```

* anvi'o – Anvi'o uses slightly different SCG criteria to assess genome completion and redundancy (what CheckM calls "contamination"). See the anvi'o tutorial for more information.

### 3.7. Dereplicating MAGs using Average Nucleotide Identity

* Average Nucleotide Identity (ANI) calculated for all orthologous genes shared between two genomes
* 95% genome-wide ANI (gANI) considered species-level delineation
* Requires both alignment fraction (AF) and alignment identity (AI) of aligned fraction
* If AF < 70%, ANI is not useful and Average Amino acid Identity (AAI) should be used

Dereplication programs:

* [dRep](https://drep.readthedocs.io/en/latest/)
dRep is very fast because it uses Mash to initially dereplicate at a coarse level, then FastANI to further dereplicate within the "primary clusters."
The below command assumes you have all your high-quality MAG fasta files (.fa extension) in one directory (e.g., "fastas"). The "map" csv file has genome ID and quality information (not used in this example). The '-sa' parameter provides the ANI threshold for dereplication.

```
dRep dereplicate dRep_out -g fastas/*.fa --genomeInfo dRep_MAGs_map.csv -sa 0.95 -nc 0.2 --ignoreGenomeQuality
```

* fastANI
* pyANI

### 3.8. Genome taxonomy

Taxonomy is assigned using sequences of SCGs

* GTDB-Tk (taxonomy)
* [anvi’o](https://merenlab.org/software/anvio/)
	- quality
	- taxonomy (GTDB but only based on ribosomal proteins, not tree-based like GTDB-Tk)

## 4. Existing wrappers or pipelines

These tools (especially anvi'o) can be incorporated into our tutorial:

* [anvi’o](https://merenlab.org/software/anvio/vignette/) (post-assembly only)
* ATLAS
* StaG-mwc
* DASTool (for binning only; dereplication and scoring of bins from multiple programs)

```
DAS_Tool -i bin_table1.tsv,bin_table2.tsv -l method1,method2 -c assembled_contigs.fasta -o output_folder -t 24
```

## 5. Features wish list

The following list includes features that we might want to include but don't have yet:

* Map reads to a set of reference genomes (e.g., the top 5 species in a sample)
* MAG functional analysis (DRAM, MEBS, MirobeAnnotator)
