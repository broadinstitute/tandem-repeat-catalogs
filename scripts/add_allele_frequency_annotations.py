import argparse
import collections
import gzip
import ijson
from intervaltree import Interval, IntervalTree
import json
import os
import pandas as pd
import re
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.file_utils import download_local_copy
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--skip-illumina174k-frequencies", action="store_true",
					help="Skip annotating with allele frequencies from Illumina 174k catalog")
parser.add_argument("--skip-t2t-assembly-frequencies", action="store_true",
					help="Skip annotating with allele frequencies from T2T assemblies")
parser.add_argument("--annotate-overlapping-loci", action="store_true",
					help="By default, this script will only annotate loci that exactly match the boundaries of loci in "
						 "the T2T assemblies catalog. This option enables annotating loci that overlap with loci in the "
						 "T2T assemblies catalog and have the same canonical motif.")
parser.add_argument("-o", "--output-path", help="Output JSON path for annotated catalog")
parser.add_argument("input_variant_catalog", help="Variant catalog in JSON or BED format")
args, _ = parser.parse_known_args()

if not args.output_path:
	args.output_path = re.sub("(.bed|.json)(.gz)?$", "", os.path.expanduser(args.input_variant_catalog))
	args.output_path += ".with_allele_frequences.json"

def convert_allele_histogram_dict_to_string(allele_histogram_dict):
	data = sorted(allele_histogram_dict.items())
	return ",".join(f"{repeat_number}x:{allele_count}" for repeat_number, allele_count in data)

def get_percentile_from_allele_histogram_dict(allele_histogram_dict, percentile):
	total = sum(allele_histogram_dict.values())
	cutoff = total * percentile
	for repeat_number, count in sorted(allele_histogram_dict.items(), reverse=True):
		total -= count
		if total <= cutoff:
			return repeat_number

def get_stdev_of_allele_histogram_dict(allele_histogram_dict):
	total = sum(allele_histogram_dict.values())
	mean = sum(repeat_number * count for repeat_number, count in allele_histogram_dict.items()) / total
	return (sum((repeat_number - mean) ** 2 * count for repeat_number, count in allele_histogram_dict.items()) / total) ** 0.5

#%%

histograms_from_illumina_174k = {}
stdev_from_illumina_174k = {}
if not args.skip_illumina174k_frequencies:
	# download illumina table
	url = "https://github.com/Illumina/RepeatCatalogs/raw/master/hg38/genotype/1000genomes/1kg.gt.hist.tsv.gz"
	print(f"Loading allele frequencies for the Illumina 174k catalog from {url}")
	df1 = pd.read_table(download_local_copy(url))
	print(f"Parsed {len(df1):,d} rows")
	print("Computing histograms")
	for _, record in df1.iterrows():
		chrom, start_0based, end = record.VariantId.split("_")
		chrom = chrom.replace("chr", "")
		repeat_numbers = [int(x) for x in record.RepeatNumbers.split(",")]
		allele_counts = [int(x) for x in record.AlleleCounts.split(",")]
		if len(repeat_numbers) != len(allele_counts):
			raise ValueError(f"RepeatNumbers and AlleleCounts have different lengths: {record.to_dict()}")

		key = (chrom, int(start_0based), int(end))
		histogram_dict = dict(zip(repeat_numbers, allele_counts))
		histograms_from_illumina_174k[key] = convert_allele_histogram_dict_to_string(histogram_dict)
		stdev_from_illumina_174k[key] = get_stdev_of_allele_histogram_dict(histogram_dict)

	print(f"Processed allele frequency histograms for {len(df1):,d} rows and computed {len(histograms_from_illumina_174k):,d} records")
df1 = None

histograms_from_t2t_assemblies = {}
stdev_from_t2t_assemblies = {}
interval_trees_for_t2t_assemblies = collections.defaultdict(IntervalTree)
if not args.skip_t2t_assembly_frequencies:
	# download table of genotypes from T2T assemblies
	url2 = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/joined.78_samples.variants.tsv.gz"
	print(f"Loading allele frequencies for the catalog of polymorphic loci in T2T assemblies from {url2}")
	df2 = pd.read_table(download_local_copy(url2))
	print(f"Parsed {len(df2):,d} rows")
	print("Computing histograms")
	allele_columns = [c for c in df2.columns if c.startswith("NumRepeats") and c != "NumRepeatsInReference"]
	df2.rename(columns={c: c.replace(":", "_") for c in allele_columns}, inplace=True)
	allele_columns = [c.replace(":", "_") for c in allele_columns]

	for record in df2.itertuples():
		chrom, start_1based, end = parse_interval(record.Locus)
		start_0based = start_1based - 1
		chrom = chrom.replace("chr", "")
		key = (chrom, start_0based, end)
		histogram_dict = collections.Counter()
		for c in allele_columns:
			allele_size = getattr(record, c) if not pd.isna(getattr(record, c)) else record.NumRepeatsInReference
			histogram_dict[int(float(allele_size))] += 1

		histograms_from_t2t_assemblies[key] = convert_allele_histogram_dict_to_string(histogram_dict)
		stdev_from_t2t_assemblies[key] = get_stdev_of_allele_histogram_dict(histogram_dict)

		if args.annotate_overlapping_loci and end > start_0based:
			interval_trees_for_t2t_assemblies[chrom].add(Interval(start_0based, end, data = {
				"HistogramDict": histogram_dict,
				"CanonicalMotif": record.CanonicalMotif,
			}))

	print(f"Processed allele frequency histograms from {len(df2):,d} rows and computed {len(histograms_from_t2t_assemblies):,d} records")
df2 = None

#%%
print(f"Parsing and annotating {args.input_variant_catalog}")
counters = collections.Counter()
fopen = gzip.open if args.output_path.endswith("gz") else open
with fopen(os.path.expanduser(args.output_path), "wt") as out:
	out.write("[")
	for i, record in enumerate(get_variant_catalog_iterator(os.path.expanduser(args.input_variant_catalog))):
		record = dict(record)
		counters["total"] += 1
		if isinstance(record["ReferenceRegion"], list):
			raise ValueError(f"ReferenceRegion is a list in {record.to_dict()}")
		chrom, start_0based, end = parse_interval(record["ReferenceRegion"])
		chrom = chrom.replace("chr", "")
		key = (chrom, start_0based, end)
		if key in histograms_from_illumina_174k:
			counters["found_illumina174_histogram"] += 1
			record["AlleleFrequenciesFromIllumina174k"] = histograms_from_illumina_174k[key]
			record["VariationStdevFromIllumina174k"] = stdev_from_illumina_174k[key]

		if key in histograms_from_t2t_assemblies:
			counters["found_t2t_assemblies_histogram"] += 1
			record["AlleleFrequenciesFromT2TAssemblies"] = histograms_from_t2t_assemblies[key]
			record["VariationStdevFromT2TAssemblies"] = stdev_from_t2t_assemblies[key]
		elif args.annotate_overlapping_loci:
			# check for overlap with nearby interval
			matching_interval = None
			for interval in interval_trees_for_t2t_assemblies[chrom].overlap(start_0based, end):
				if not record.get("CanonicalMotif"):
					record["CanonicalMotif"] = compute_canonical_motif(record["LocusStructure"].strip("()*+").upper())

				if interval.data["CanonicalMotif"] != record["CanonicalMotif"]:
					continue
				if interval.overlap_size(start_0based, end) < 2*len(record["CanonicalMotif"]):
					continue
				if (interval.begin > start_0based and interval.end > end) or (interval.begin < start_0based and interval.end < end):
					continue

				matching_interval = interval
				break

			if matching_interval:
				counters["found_t2t_assemblies_histogram_via_overlap"] += 1
				histogram_dict = matching_interval.data["HistogramDict"]
				locus_boundary_diff = ((matching_interval.end - matching_interval.begin) - (end - start_0based))//len(record["CanonicalMotif"])

				# Adjust the genotype repeat count by the difference in locus boundaries since changes in locus
				# bondaries affect the overall repeat count in each allele.
				histogram_dict_adjusted = {
					max(0, repeat_number - locus_boundary_diff): count for repeat_number, count in histogram_dict.items()
				}

				record["AlleleFrequenciesFromT2TAssemblies"] = convert_allele_histogram_dict_to_string(histogram_dict_adjusted)
				record["VariationStdevFromT2TAssemblies"] = get_stdev_of_allele_histogram_dict(histogram_dict_adjusted)

		if i > 0:
			out.write(", ")
		json.dump(record, out, indent=1)
	out.write("]\n")

print(f"Annotated {counters['found_illumina174_histogram']:,d} out of {counters['total']:,d} loci in the Illumina 174k allele frequency catalog")
print(f"Annotated {counters['found_t2t_assemblies_histogram']:,d} out of {counters['total']:,d} loci in the T2T assemblies allele frequency catalog")
if args.annotate_overlapping_loci:
	print(f"Annotated {counters['found_t2t_assemblies_histogram_via_overlap']:,d} out of {counters['total']:,d} loci in the T2T assemblies allele frequency catalog based on overlap")
print(f"Wrote {counters['total']:,d} records to {args.output_path}")

#%%



# gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/combined/joined_only_pure_repeats.78_samples.variants.tsv.gz