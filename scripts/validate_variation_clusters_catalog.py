"""Checks a variation_clusters_and_all_repeats TRGT catalog against the repeat catalog it was built
from.

The catalog holds one row per repeat plus one extra row per variation cluster, so the repeat rows
should reproduce the repeat catalog's own TRGT file exactly, apart from the VC_FILTER key that marks
loci the variation clusters table could not process.
"""

import argparse
import collections
import gzip
import re


def open_maybe_gzipped(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def strip_filter_key(info_field):
    """Returns the info field without the VC_FILTER key that this catalog adds to flagged loci."""
    return re.sub(r";VC_FILTER=[^;]*", "", info_field)


def locus_id_of(info_field):
    return re.search(r"ID=([^;]+)", info_field).group(1)


def is_variation_cluster(info_field):
    return locus_id_of(info_field).startswith("VC:")


def read_catalog(path):
    """Reads a TRGT catalog BED file.

    Args:
        path (str): Path of the BED file, optionally gzipped.

    Returns:
        tuple: (row counter keyed by the repeat rows with VC_FILTER stripped, cluster row count,
            counter of VC_FILTER values, counter of locus ids, VC_FILTER value by locus id for the
            flagged rows only)
    """
    repeat_rows = collections.Counter()
    filter_values = collections.Counter()
    locus_ids = collections.Counter()
    filter_by_locus_id = {}
    cluster_row_count = 0
    with open_maybe_gzipped(path) as f:
        for line in f:
            chrom, start, end, info_field = line.rstrip("\n").split("\t")[:4]
            locus_ids[locus_id_of(info_field)] += 1
            if is_variation_cluster(info_field):
                cluster_row_count += 1
                continue
            match = re.search(r";VC_FILTER=([^;]*)", info_field)
            if match:
                filter_values[match.group(1)] += 1
                filter_by_locus_id[locus_id_of(info_field)] = match.group(1)
            repeat_rows[f"{chrom}\t{start}\t{end}\t{strip_filter_key(info_field)}"] += 1
    return repeat_rows, cluster_row_count, filter_values, locus_ids, filter_by_locus_id


def check_repeat_rows_match_repeat_catalog(repeat_rows, repeat_catalog_path):
    """Returns a list of failure messages, empty if the repeat rows reproduce the repeat catalog."""
    remaining = collections.Counter(repeat_rows)
    missing = 0
    with open_maybe_gzipped(repeat_catalog_path) as f:
        for line in f:
            row = line.rstrip("\n")
            if remaining[row] > 0:
                remaining[row] -= 1
            else:
                missing += 1
    extra = sum(count for count in remaining.values() if count > 0)
    failures = []
    if missing:
        failures.append(f"{missing:,d} rows in {repeat_catalog_path} have no matching repeat row")
    if extra:
        example = next(row for row, count in remaining.items() if count > 0)
        failures.append(f"{extra:,d} repeat rows are not in {repeat_catalog_path}, e.g. {example}")
    return failures


def check_locus_ids_are_unique(locus_ids):
    duplicates = [locus_id for locus_id, count in locus_ids.items() if count > 1]
    if not duplicates:
        return []
    return [f"{len(duplicates):,d} locus ids appear on more than one row, e.g. {duplicates[:3]}"]


def check_catalog_adds_clusters_and_flags(cluster_row_count, filter_values):
    """Guards against a file that is really just a copy of the repeat catalog.

    The cluster rows and the VC_FILTER keys are the only things this catalog adds over the repeat
    catalog's own TRGT file, and check_repeat_rows_match_repeat_catalog passes either way, so
    without this a catalog missing both would be reported as correct.
    """
    failures = []
    if cluster_row_count == 0:
        failures.append("the catalog has no variation cluster rows")
    if not filter_values:
        failures.append("no row carries a VC_FILTER key, so flagged loci are not being marked")
    return failures


def check_expected_filtered_locus(filter_by_locus_id, expected_locus_id, expected_filter):
    """Guards the regression where loci flagged by the clusters table vanished from the catalog."""
    actual = filter_by_locus_id.get(expected_locus_id)
    if actual == expected_filter:
        return []
    if actual is None:
        return [f"{expected_locus_id} carries no VC_FILTER key, so {expected_filter}-flagged loci "
                f"are being dropped or left unflagged rather than marked"]
    return [f"{expected_locus_id} carries VC_FILTER={actual}, expected {expected_filter}"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expect-filtered-locus", default="1-146228800-146228821-GCC",
                        help="Locus id that the variation clusters table flags and that must still "
                             "have a row of its own. Defaults to the NOTCH2NLA repeat.")
    parser.add_argument("variation_clusters_and_all_repeats_bed_path",
                        help="Path of the combined TRGT catalog to check")
    parser.add_argument("repeat_catalog_trgt_bed_path",
                        help="Path of the repeat catalog's own TRGT BED file")
    args = parser.parse_args()

    repeat_rows, cluster_row_count, filter_values, locus_ids, filter_by_locus_id = read_catalog(
        args.variation_clusters_and_all_repeats_bed_path)

    print(f"{sum(repeat_rows.values()):,d} repeat rows, {cluster_row_count:,d} variation cluster rows")
    for value, count in sorted(filter_values.items()):
        print(f"  {count:>12,d}  rows flagged VC_FILTER={value}")

    failures = (
        check_repeat_rows_match_repeat_catalog(repeat_rows, args.repeat_catalog_trgt_bed_path)
        + check_locus_ids_are_unique(locus_ids)
        + check_catalog_adds_clusters_and_flags(cluster_row_count, filter_values)
        + check_expected_filtered_locus(filter_by_locus_id, args.expect_filtered_locus, "DEPTH")
    )

    for failure in failures:
        print(f"FAILED: {failure}")
    if failures:
        raise SystemExit(1)
    print("All checks passed")


if __name__ == "__main__":
    main()
