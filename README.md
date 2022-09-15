# Osprey

Osprey is software for phasing or imputing copy number variants using SNP genotypes.

## History

Osprey is derived from methods developed for imputation of VNTRs in a collaboration between Bob Handsaker
in the lab of Steve McCarroll (http://mccarrolllab.org) and Ronen Mukamel and Po-Ru Loh in the Loh lab
(http://statgen.hms.harvard.edu/) as described in Protein-coding repeat polymorphisms strongly shape diverse human phenotypes
(https://www.science.org/doi/10.1126/science.abg8289, PMID:34554798 https://pubmed.ncbi.nlm.nih.gov/34554798).
Osprey was strongly influenced by the MVNCALL software (Menelaou and Marchini, 2013, https://doi.org/10.1093/bioinformatics/bts632).

Conceptually, at each locus of interest Osprey matches the local phased SNP haplotypes in a reference panel to local phased SNP haplotypes
in a target panel. Copy number variants present in the reference panel are then imputed into the target panel based on the best-matching
haplotypes. This is done for one copy number variant at a time (or one locus with multiple copy number variants at a time).

## Usage

Osprey consists of two main programs, **ospreyIBS** and **osprey**, that are separated for ease of use.
The ospreyIBS program does haplotype matching and produces a matrix of the best N matching reference and target haplotypes for each of the haplotypes in the target panel.
The osprey program them uses this matrix to impute a set of copy number variants using the IBS (identity-by-state) matrix produced by ospreyIBS.

**ospreyIBS arguments**

|Argument|Description|
|:-------|:----------|
|--inputFile [-vcf] file           |Input VCF file of phased SNP genotypes (required).
|--ibsMatrixFile [-ibs] file.gz    |Output IBS matrix file (required).
|--ibsInterval [-L] interval       |Target genomic interval (required).
|--geneticMapFile [-gmap] file     |Genetic map file (required). The shapeit file format is expected.
|--refSampleFile [-rs] file        |File listing the reference samples to use for matching (one per line).
|--threads [-t] N                  |Number of threads to use (default 1).
|--version                         |Print software version.
|--help [-h]                       |Print help synopsis.
|--verbose N                       |Verbosity level.
|--debug N                         |Debug output level.

The input file must be in VCF format (optionally compressed). The genotypes must be phased.

The output ibs matrix is gzip compressed if the file name ends in .gz (recommended).

The target genomic interval indicates the region where the copy number variation occurs.
Haplotype matching is done outside of the target region.

Currently only diploid genomic regions are supported. The organism's genome must be diploid and the target chromosome must be diploid.

**osprey arguments**

|Argument|Description|
|:-------|:----------|
|--inputFile [-vcf] file         |Input VCF file of copy number variants (required).
|--outputFile [-o] file          |Output VCF file (required).
|--ibsMatrixFile [-ibs] file     |IBS matrix produced by ospreyIBS (required).
|--threads [-t] N                |Number of threads to use (default 1).
|--iterations [-iter] N          |Number of phasing iterations (default 250).
|--site ID                       |Specify only one input site to process by variant ID or multiple sites comma-separated.
|--siteIndex N@B or N[-M]        |Specify sites to process for parallelization (see below).
|--version                       |Print software version.
|--help [-h]                     |Print help synopsis.
|--verbose N                     |Verbosity level.
|--debug N                       |Debug output level.

The input file must be in VCF format (optionally compresed).
The input VCF is assumed to contain genotypes for a set of reference samples with copy number likelihoods that will be phased
and optionally a set of target samples with no copy number likelihoods that will be imputed.
Each input record should be a copy number variant with copy number likelihoods, if present, encoded in the GT (FORMAT) tags CNL or CNP.
If both CNL and CNP tags are present, CNL is used.

The output file is a VCF file, which will be compressed if the file name ends in .gz, of each phased or imputed variant.
This will be all variants in the input file, unless a subset of variant sites is selected with `--site` or `--siteIndex`.

The phasing/imputation results are encoded in the following GT (FORMAT) tags:

|Tag|Description|
|:---|:---|
|PST|A character value of either P or I indicating whether the samples was phased or imputed.|
|PCN|This indicates the phased (haploid) copy numbers in the format "n\|n". The order corresponds to the phased SNP genotypes.|
|PCNF|Phased copy number dosage estimates in the format "n.nn\|n.nn".|
|PCNQ|A Phred-scaled estimate of the confidence in the phasing/imputation quality, one value for each haploptype, in the format "q1\|q2".|
|PCNL|Two vectors of copy number likelihoods, one for each haplotype, separated by vertical bar.|

The last four fields are analogous to the copy-number genotype fields CN, CNF, CNQ and CNL used by GenomeSTRiP for diploid copy number.

## Parallelization

The osprey program has built-in support for parallelizing workflows by processing disjoint subsets of the input VCF in different runs.
This is typically done using `--siteIndex N@B` where N is the batch number and B is the number of batches desired.
For example, to process an input VCF file using 10 parallel processes, use `--siteIndex k@10` and run 10 processes with k ranging from 1 to 10
and naming the output files appropriately (e.g. with k@10 in the name).
The first process will process the first batch of the input sites (approximately 10% of them), the second process will do the second batch, and so on.
The output VCFs can then be concatenated with `bcftools concat`.

It is also possible to select a single site to process or a specific set of sites using `--site` and specifying the ID of the variant
or a comma-separated list of variant names. Alternatively, `--siteIndex N` will process the Nth variant in the input file (and `--siteIndex
N-M` will process all sites with indexes from N to M numbered from one). This functionality is more often used for debugging than for parallelization.

The combination of ospreyIBS and osprey are intended to be used on one locus (one target interval) at a time.
Parallelization at the locus level can be done by scripting the set of target intervals to process, but there is no built-in support
for this level of parallel processing.

## Credits

Bob Handsaker is the main developer/maintainer of Osprey.

Po-Ru Loh wrote the first version of the core algorithm.

Ronen Mukamel did much of the analysis and methods development on the VNTR imputation algorithm that preceded Osprey.
