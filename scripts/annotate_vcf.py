#!/usr/bin/env python3
"""
VCF Annotation Helper

Utility script to prepare VCF files for the pathway subtyping framework.
Adds required INFO fields (GENE, CONSEQUENCE, CADD) from various annotation sources.

Usage:
    python scripts/annotate_vcf.py input.vcf output.vcf --vep-tsv annotations.tsv
    python scripts/annotate_vcf.py input.vcf output.vcf --annovar-file input.hg38_multianno.txt
"""

import argparse
import gzip
import logging
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_vep_tsv(tsv_path: str) -> Dict[str, Tuple[str, str, float]]:
    """
    Parse VEP TSV output to extract gene, consequence, and CADD.

    Expected columns: #Uploaded_variation, Gene, Consequence, CADD_PHRED
    """
    annotations = {}

    with open(tsv_path) as f:
        header = None
        for line in f:
            if line.startswith("#Uploaded"):
                header = line.strip().split("\t")
                continue
            if header is None or line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < len(header):
                continue

            row = dict(zip(header, parts))

            var_id = row.get("#Uploaded_variation", "")
            gene = row.get("Gene", row.get("SYMBOL", ""))
            consequence = row.get("Consequence", "").split(",")[0]
            cadd = row.get("CADD_PHRED", "0")

            try:
                cadd_score = float(cadd) if cadd and cadd != "-" else 0.0
            except ValueError:
                cadd_score = 0.0

            # Use variant ID as key
            if var_id and gene:
                annotations[var_id] = (gene, consequence, cadd_score)

    logger.info(f"Loaded {len(annotations)} annotations from VEP TSV")
    return annotations


def parse_annovar_multianno(annovar_path: str) -> Dict[str, Tuple[str, str, float]]:
    """
    Parse ANNOVAR multianno output.

    Expected format: table_annovar.pl output with Gene.refGene, Func.refGene, CADD_phred
    """
    annotations = {}

    with open(annovar_path) as f:
        header = None
        for line in f:
            if header is None:
                header = line.strip().split("\t")
                continue

            parts = line.strip().split("\t")
            if len(parts) < len(header):
                continue

            row = dict(zip(header, parts))

            # Build variant key (chr:pos:ref:alt)
            chrom = row.get("Chr", "")
            pos = row.get("Start", "")
            ref = row.get("Ref", "")
            alt = row.get("Alt", "")
            var_key = f"{chrom}:{pos}:{ref}:{alt}"

            gene = row.get("Gene.refGene", row.get("Gene", "")).split(";")[0]
            func = row.get("ExonicFunc.refGene", row.get("Func.refGene", ""))
            cadd = row.get("CADD_phred", row.get("CADD13_PHRED", "0"))

            # Map ANNOVAR function to VEP-style consequence
            consequence = _map_annovar_func(func)

            try:
                cadd_score = float(cadd) if cadd and cadd != "." else 0.0
            except ValueError:
                cadd_score = 0.0

            if gene and gene != ".":
                annotations[var_key] = (gene, consequence, cadd_score)

    logger.info(f"Loaded {len(annotations)} annotations from ANNOVAR")
    return annotations


def _map_annovar_func(func: str) -> str:
    """Map ANNOVAR ExonicFunc to VEP-style consequence."""
    mapping = {
        "frameshift insertion": "frameshift_variant",
        "frameshift deletion": "frameshift_variant",
        "frameshift block substitution": "frameshift_variant",
        "stopgain": "stop_gained",
        "stoploss": "stop_lost",
        "nonsynonymous SNV": "missense_variant",
        "synonymous SNV": "synonymous_variant",
        "splicing": "splice_acceptor_variant",
        "nonframeshift insertion": "inframe_insertion",
        "nonframeshift deletion": "inframe_deletion",
    }
    return mapping.get(func, "unknown")


def annotate_vcf(
    input_vcf: str,
    output_vcf: str,
    annotations: Dict[str, Tuple[str, str, float]],
) -> int:
    """
    Add annotations to VCF file.

    Returns number of annotated variants.
    """
    open_func = gzip.open if input_vcf.endswith(".gz") else open
    write_func = gzip.open if output_vcf.endswith(".gz") else open

    annotated_count = 0
    total_count = 0

    with open_func(input_vcf, "rt") as fin, write_func(output_vcf, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                # Add INFO header lines before #CHROM
                if line.startswith("#CHROM"):
                    fout.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n')
                    fout.write('##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">\n')
                    fout.write('##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD phred score">\n')
                fout.write(line)
                continue

            total_count += 1
            parts = line.strip().split("\t")
            if len(parts) < 8:
                fout.write(line)
                continue

            chrom, pos, vid, ref, alt = parts[:5]

            # Try different key formats
            var_keys = [
                vid,  # ID field
                f"{chrom}:{pos}:{ref}:{alt}",  # chr:pos:ref:alt
                f"{chrom}_{pos}_{ref}_{alt}",  # chr_pos_ref_alt
            ]

            gene, consequence, cadd = None, None, None
            for key in var_keys:
                if key in annotations:
                    gene, consequence, cadd = annotations[key]
                    break

            if gene:
                # Add annotations to INFO field
                info = parts[7]
                new_info = f"{info};GENE={gene};CONSEQUENCE={consequence};CADD={cadd:.1f}"
                parts[7] = new_info
                annotated_count += 1

            fout.write("\t".join(parts) + "\n")

    logger.info(f"Annotated {annotated_count}/{total_count} variants")
    return annotated_count


def validate_vcf(vcf_path: str) -> Dict[str, int]:
    """
    Validate VCF has required annotations.

    Returns dict with counts of each annotation field.
    """
    counts = {"total": 0, "GENE": 0, "CONSEQUENCE": 0, "CADD": 0}

    open_func = gzip.open if vcf_path.endswith(".gz") else open

    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            counts["total"] += 1
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            info = parts[7]
            if "GENE=" in info:
                counts["GENE"] += 1
            if "CONSEQUENCE=" in info:
                counts["CONSEQUENCE"] += 1
            if "CADD=" in info:
                counts["CADD"] += 1

    return counts


def main():
    parser = argparse.ArgumentParser(
        description="Annotate VCF with GENE, CONSEQUENCE, CADD for pathway subtyping"
    )
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output annotated VCF file")
    parser.add_argument("--vep-tsv", help="VEP TSV annotation file")
    parser.add_argument("--annovar-file", help="ANNOVAR multianno file")
    parser.add_argument("--validate-only", action="store_true", help="Only validate, don't annotate")

    args = parser.parse_args()

    if args.validate_only:
        counts = validate_vcf(args.input_vcf)
        print(f"VCF Validation Report for: {args.input_vcf}")
        print(f"  Total variants: {counts['total']}")
        print(f"  With GENE: {counts['GENE']} ({100*counts['GENE']/max(1,counts['total']):.1f}%)")
        print(f"  With CONSEQUENCE: {counts['CONSEQUENCE']} ({100*counts['CONSEQUENCE']/max(1,counts['total']):.1f}%)")
        print(f"  With CADD: {counts['CADD']} ({100*counts['CADD']/max(1,counts['total']):.1f}%)")

        if counts["GENE"] < counts["total"] * 0.5:
            print("\nWARNING: Less than 50% of variants have GENE annotation")
            sys.exit(1)
        return

    # Load annotations
    annotations = {}
    if args.vep_tsv:
        annotations = parse_vep_tsv(args.vep_tsv)
    elif args.annovar_file:
        annotations = parse_annovar_multianno(args.annovar_file)
    else:
        parser.error("Must provide either --vep-tsv or --annovar-file")

    # Annotate VCF
    annotate_vcf(args.input_vcf, args.output_vcf, annotations)
    print(f"Annotated VCF written to: {args.output_vcf}")

    # Validate output
    counts = validate_vcf(args.output_vcf)
    print(f"Annotation summary:")
    print(f"  With GENE: {counts['GENE']}/{counts['total']}")


if __name__ == "__main__":
    main()
