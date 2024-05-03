### Authors

* Ben Weisburd, Broad Institute
* Egor Dolzhenko, PacBio

### Introduction

Tandem repeats are regions of the genome that consist of consecutive copies of some motif sequence. For example, CAGCAGCAG is a tandem repeats with motif CAG. Many types of genomic studies require annotations of tandem repeats in the reference genome, called repeat catalogs. Repeat catalogs typically consist of the repeat reference coordinates and one or multiple motif sequences that the repeat is composed of.

The purpose of this repo is to provides methods and best practices for defining high-quality repeat catalogs. Our initial focus will be on the human genome, however the tutorial should also be applicable to other closely-related genomes. We'd also love to extend this work to plants and other species. Please consider creating a GitHub issue or reaching out by email if you are interested in this.

### Defining a genome-wide tandem repeat catalog

The steps below start with a reference genome fasta file and proceed to generate a comprehensive genome-wide TR catalog that can be used as input to the following TR genotyping tools:
- ExpansionHunter
- GangSTR
- HipSTR
- TRGT
- LongTR


### Step 1: Detect all perfect (ie. non-interrupted) tandem repeats in the reference genome

Run [colab-repeat-finder](https://github.com/broadinstitute/colab-repeat-finder) on the reference genome: 

```
git clone git@github.com:broadinstitute/colab-repeat-finder.git
cd colab-repeat-finder/python
python3 perfect_repeat_finder.py  --min-repeats 3  --min-span 9  --min-motif-size 1  --max-motif-size 50  --output-prefix perfect_repeats.hg38  --show-progress-bar   /path/to/hg38.fa   
```

The output of this for hg38 is already available @
https://storage.cloud.google.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp.bed.gz

*NOTE:* colab-repeat-finder is a new, simple tool that finds all perfect repeats in a given input sequence that satify user parameters. We could have, instead, run TandemRepeatFinder (TRF) with very large mismatch and indel penalties to make it return only perfect repeats, but in this mode TRF fails to detect ~3% of perfect repeats (2-50bp motifs) for unclear reasons. We therefore created colab-repeat-finder to maximise sensitivity. 


### Step 2: Merge loci from step 1 (and optionally add other repeat catalogs)

Here, we fix some representation issues by doing the following:
- collapsing adjacent loci that have the same motif but were reported as two separate loci by colab-repeat-finder due to a single base pair interruption between them.
- filtering out repeats that include non-A,C,G,T characters (such as homopolymer repeats of 'N')

Also, we can optionally augment our catalog with loci from additional sources. Specifically, while the catalog from step 1 includes all tandem repeat loci that have at least 3 repeats of some motif and span at least 9bp in the reference genome, it misses loci that have 2 or fewer repeats in the reference (while having 3 or more repeats in other genomes within the population). To capture these loci as well, we can merge the catalog from step 1 with any available catalogs of polymorphic tandem repeat loci that were generated via 
orthogonal methods. For the human genome, these include the following catalogs:

* Illumina catalog of 174k TR loci that are polymporphic in 2.4k diverse population samples from 1kGP.  
* Truth set from [weisburd 2023] of polymorphic TR loci in 51 samples from the HPRC. 


```
ALL_PERFECT_REPEATS_CATALOG_URL=https://storage.cloud.google.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_50bp.repeats_3x_and_spans_9bp.bed.gz
ILLUMINA_CATALOG_URL=https://storage.cloud.google.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz
TRUTH_SET_CATALOG_URL=https://storage.cloud.google.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/combined.51_samples.variants.bed.gz


wget $ALL_PERFECT_REPEATS_CATALOG_URL
wget $ILLUMINA_CATALOG_URL
wget $TRUTH_SET_CATALOG_URL

# TODO rename to merge_loci
python3 -u -m str_analysis.merge_loci --verbose \
  --merge-adjacent-loci-with-same-motif \
  $(ALL_PERFECT_REPEATS_CATALOG_URL) \
  $(basename $ILLUMINA_CATALOG_URL) \
  $(basename TRUTH_SET_CATALOG_URL)
```



### Step 3: Annotate loci with gene names, etc.

Use the `python3 -u -m str_analysis.annotate_and_filter_str_catalog` script.

### Step 4: Convert the combined catalog into tool-specific catalog formats, while adding tool-specific optimizations to the locus definitions 

Use the scrips in the str_analysis repo.

