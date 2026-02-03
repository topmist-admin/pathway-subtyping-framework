#!/usr/bin/env python3
"""
VCF Annotation Helper (v0.2)

Utility script to prepare VCF files for the pathway subtyping framework.
Adds required INFO fields (GENE, CONSEQUENCE, CADD) from various annotation sources.

Features:
- Support for VEP TSV and ANNOVAR multianno formats
- Multi-allelic variant handling
- Validation mode to check existing annotations
- Detailed error messages with fix suggestions

Usage:
    # Annotate using VEP output
    python scripts/annotate_vcf.py input.vcf output.vcf --vep-tsv annotations.tsv

    # Annotate using ANNOVAR output
    python scripts/annotate_vcf.py input.vcf output.vcf --annovar-file input.hg38_multianno.txt

    # Validate existing VCF annotations
    python scripts/annotate_vcf.py input.vcf --validate-only

    # Show detailed statistics
    python scripts/annotate_vcf.py input.vcf output.vcf --vep-tsv annotations.tsv --verbose

For more help: python scripts/annotate_vcf.py --help
"""

import argparse
import gzip
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


class AnnotationError(Exception):
    """Exception raised for annotation-related errors."""

    def __init__(self, message: str, suggestions: List[str]):
        self.message = message
        self.suggestions = suggestions
        super().__init__(self._format())

    def _format(self) -> str:
        lines = [self.message, "", "Suggestions:"]
        for i, s in enumerate(self.suggestions, 1):
            lines.append(f"  {i}. {s}")
        return "\n".join(lines)


def parse_vep_tsv(tsv_path: str, verbose: bool = False) -> Dict[str, Tuple[str, str, float]]:
    """
    Parse VEP TSV output to extract gene, consequence, and CADD.

    Supported VEP output formats:
    - Standard TSV (--tab)
    - With CADD plugin
    - With LoFtool scores

    Expected columns: #Uploaded_variation, Gene, Consequence, CADD_PHRED (optional)
    """
    annotations = {}
    stats = Counter()

    path = Path(tsv_path)
    if not path.exists():
        raise AnnotationError(
            f"VEP TSV file not found: {tsv_path}",
            [
                "Check the file path is correct",
                "Run VEP first: vep -i input.vcf -o output.tsv --tab --symbol --fields "
                "'Uploaded_variation,Gene,SYMBOL,Consequence,CADD_PHRED'",
            ],
        )

    with open(tsv_path) as f:
        header = None
        for line_num, line in enumerate(f, 1):
            # Handle VEP header comments
            if line.startswith("##"):
                continue

            if line.startswith("#Uploaded") or line.startswith("#Location"):
                header = line.strip().lstrip("#").split("\t")
                continue

            if header is None:
                if line_num == 1:
                    # Try to use first line as header
                    header = line.strip().split("\t")
                    continue
                raise AnnotationError(
                    f"Could not find header line in VEP TSV: {tsv_path}",
                    [
                        "VEP TSV should have a header line starting with '#Uploaded_variation'",
                        "Re-run VEP with: vep -i input.vcf -o output.tsv --tab",
                    ],
                )

            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < len(header):
                stats["skipped_short"] += 1
                continue

            row = dict(zip(header, parts))
            stats["total"] += 1

            # Extract variant ID (try multiple column names)
            var_id = (
                row.get("Uploaded_variation")
                or row.get("Location")
                or row.get("#Uploaded_variation")
                or ""
            )

            # Extract gene symbol (try multiple column names)
            gene = row.get("SYMBOL") or row.get("Gene") or row.get("SYMBOL_SOURCE") or ""

            # Extract consequence
            consequence = row.get("Consequence", "").split(",")[0]

            # Extract CADD score (try multiple column names)
            cadd_str = (
                row.get("CADD_PHRED") or row.get("CADD_phred") or row.get("CADD") or "0"
            )

            try:
                cadd_score = float(cadd_str) if cadd_str and cadd_str not in ("-", ".") else 0.0
            except ValueError:
                cadd_score = 0.0

            # Track statistics
            if gene and gene not in ("-", "."):
                stats["has_gene"] += 1
            if consequence and consequence not in ("-", "."):
                stats["has_consequence"] += 1
            if cadd_score > 0:
                stats["has_cadd"] += 1

            # Use variant ID as key (handle VEP-style IDs)
            if var_id and gene and gene not in ("-", "."):
                # VEP may produce multiple entries per variant (multiple transcripts)
                # Keep the first one or one with highest CADD
                if var_id not in annotations or cadd_score > annotations[var_id][2]:
                    annotations[var_id] = (gene, consequence, cadd_score)

    if verbose:
        total = max(stats.get("total", 1), 1)
        logger.info("VEP TSV Statistics:")
        logger.info(f"  Total entries: {stats.get('total', 0)}")
        logger.info(f"  With gene: {stats.get('has_gene', 0)} ({100*stats.get('has_gene', 0)/total:.1f}%)")
        logger.info(
            f"  With consequence: {stats.get('has_consequence', 0)} "
            f"({100*stats.get('has_consequence', 0)/total:.1f}%)"
        )
        logger.info(f"  With CADD: {stats.get('has_cadd', 0)} ({100*stats.get('has_cadd', 0)/total:.1f}%)")
        logger.info(f"  Unique variants: {len(annotations)}")

    logger.info(f"Loaded {len(annotations)} annotations from VEP TSV")
    return annotations


def parse_annovar_multianno(annovar_path: str, verbose: bool = False) -> Dict[str, Tuple[str, str, float]]:
    """
    Parse ANNOVAR multianno output.

    Supported formats:
    - table_annovar.pl output
    - With CADD scores (CADD_phred column)

    Expected format: Gene.refGene, ExonicFunc.refGene, CADD_phred
    """
    annotations = {}
    stats = Counter()

    path = Path(annovar_path)
    if not path.exists():
        raise AnnotationError(
            f"ANNOVAR file not found: {annovar_path}",
            [
                "Check the file path is correct",
                "Run ANNOVAR first: table_annovar.pl input.avinput humandb/ "
                "-buildver hg38 -out output -protocol refGene,cadd -operation g,f",
            ],
        )

    with open(annovar_path) as f:
        header = None
        for line_num, line in enumerate(f, 1):
            if header is None:
                header = line.strip().split("\t")
                # Validate required columns
                required = {"Chr", "Start", "Ref", "Alt"}
                found = set(header) & required
                if len(found) < 3:
                    raise AnnotationError(
                        f"ANNOVAR file missing required columns: {required - found}",
                        [
                            "Ensure ANNOVAR output has Chr, Start, Ref, Alt columns",
                            "Check that you used the correct output format",
                        ],
                    )
                continue

            parts = line.strip().split("\t")
            if len(parts) < len(header):
                stats["skipped_short"] += 1
                continue

            row = dict(zip(header, parts))
            stats["total"] += 1

            # Build variant key (chr:pos:ref:alt)
            chrom = row.get("Chr", "")
            pos = row.get("Start", "")
            ref = row.get("Ref", "")
            alt = row.get("Alt", "")
            var_key = f"{chrom}:{pos}:{ref}:{alt}"

            # Extract gene (handle multiple genes)
            gene = row.get("Gene.refGene", row.get("Gene", "")).split(";")[0]
            if gene in (".", "NONE", ""):
                gene = row.get("Gene.ensGene", "").split(";")[0]

            # Extract functional annotation
            func = row.get("ExonicFunc.refGene", row.get("Func.refGene", ""))

            # Map ANNOVAR function to VEP-style consequence
            consequence = _map_annovar_func(func)

            # Extract CADD score (try multiple column names)
            cadd_str = row.get("CADD_phred") or row.get("CADD13_PHRED") or row.get("CADD") or "0"

            try:
                cadd_score = float(cadd_str) if cadd_str and cadd_str not in (".", "") else 0.0
            except ValueError:
                cadd_score = 0.0

            # Track statistics
            if gene and gene not in (".", "NONE"):
                stats["has_gene"] += 1
            if consequence and consequence != "unknown":
                stats["has_consequence"] += 1
            if cadd_score > 0:
                stats["has_cadd"] += 1

            if gene and gene not in (".", "NONE"):
                annotations[var_key] = (gene, consequence, cadd_score)

    if verbose:
        total = max(stats.get("total", 1), 1)
        logger.info("ANNOVAR Statistics:")
        logger.info(f"  Total entries: {stats.get('total', 0)}")
        logger.info(f"  With gene: {stats.get('has_gene', 0)} ({100*stats.get('has_gene', 0)/total:.1f}%)")
        logger.info(
            f"  With consequence: {stats.get('has_consequence', 0)} "
            f"({100*stats.get('has_consequence', 0)/total:.1f}%)"
        )
        logger.info(f"  With CADD: {stats.get('has_cadd', 0)} ({100*stats.get('has_cadd', 0)/total:.1f}%)")
        logger.info(f"  Unique variants: {len(annotations)}")

    logger.info(f"Loaded {len(annotations)} annotations from ANNOVAR")
    return annotations


def _map_annovar_func(func: str) -> str:
    """Map ANNOVAR ExonicFunc to VEP-style consequence."""
    if not func:
        return "unknown"

    func_lower = func.lower()
    mapping = {
        "frameshift insertion": "frameshift_variant",
        "frameshift deletion": "frameshift_variant",
        "frameshift block substitution": "frameshift_variant",
        "frameshift_variant": "frameshift_variant",
        "stopgain": "stop_gained",
        "stop_gained": "stop_gained",
        "stoploss": "stop_lost",
        "stop_lost": "stop_lost",
        "nonsynonymous snv": "missense_variant",
        "missense_variant": "missense_variant",
        "synonymous snv": "synonymous_variant",
        "synonymous_variant": "synonymous_variant",
        "splicing": "splice_acceptor_variant",
        "splice_acceptor_variant": "splice_acceptor_variant",
        "splice_donor_variant": "splice_donor_variant",
        "nonframeshift insertion": "inframe_insertion",
        "inframe_insertion": "inframe_insertion",
        "nonframeshift deletion": "inframe_deletion",
        "inframe_deletion": "inframe_deletion",
        "startloss": "start_lost",
        "start_lost": "start_lost",
    }

    for key, value in mapping.items():
        if key in func_lower:
            return value

    return "unknown"


def annotate_vcf(
    input_vcf: str,
    output_vcf: str,
    annotations: Dict[str, Tuple[str, str, float]],
    verbose: bool = False,
) -> Dict[str, int]:
    """
    Add annotations to VCF file.

    Handles:
    - Gzipped and plain VCF files
    - Multi-allelic variants
    - Multiple key formats for matching

    Returns dict with annotation statistics.
    """
    stats = {
        "total": 0,
        "annotated": 0,
        "multi_allelic": 0,
        "already_annotated": 0,
    }

    open_func = gzip.open if input_vcf.endswith(".gz") else open
    write_func = gzip.open if output_vcf.endswith(".gz") else open

    with open_func(input_vcf, "rt") as fin, write_func(output_vcf, "wt") as fout:
        header_written = False

        for line in fin:
            if line.startswith("#"):
                # Add INFO header lines before #CHROM
                if line.startswith("#CHROM") and not header_written:
                    fout.write(
                        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol from annotation">\n'
                    )
                    fout.write(
                        '##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">\n'
                    )
                    fout.write('##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD phred score">\n')
                    header_written = True
                fout.write(line)
                continue

            stats["total"] += 1
            parts = line.strip().split("\t")
            if len(parts) < 8:
                fout.write(line)
                continue

            chrom, pos, vid, ref, alts_str = parts[:5]
            info = parts[7]

            # Check if already annotated
            if "GENE=" in info:
                stats["already_annotated"] += 1
                fout.write(line)
                continue

            # Handle multi-allelic
            alts = alts_str.split(",")
            if len(alts) > 1:
                stats["multi_allelic"] += 1

            # Try different key formats to find annotation
            gene, consequence, cadd = None, None, None

            for alt in alts:
                var_keys = [
                    vid,  # ID field
                    f"{chrom}:{pos}:{ref}:{alt}",  # chr:pos:ref:alt
                    f"{chrom}_{pos}_{ref}_{alt}",  # chr_pos_ref_alt (VEP style)
                    f"{chrom}_{pos}_{ref}/{alt}",  # chr_pos_ref/alt (alternative)
                    f"{chrom}:{pos}:{ref}:{alts_str}",  # with all alts
                ]

                for key in var_keys:
                    if key in annotations:
                        gene, consequence, cadd = annotations[key]
                        break
                if gene:
                    break

            if gene:
                # Add annotations to INFO field
                new_info = f"{info};GENE={gene};CONSEQUENCE={consequence};CADD={cadd:.1f}"
                parts[7] = new_info
                stats["annotated"] += 1

            fout.write("\t".join(parts) + "\n")

    if verbose:
        logger.info("Annotation Statistics:")
        logger.info(f"  Total variants: {stats['total']}")
        logger.info(
            f"  Annotated: {stats['annotated']} ({100*stats['annotated']/max(1, stats['total']):.1f}%)"
        )
        logger.info(f"  Already annotated: {stats['already_annotated']}")
        logger.info(f"  Multi-allelic: {stats['multi_allelic']}")

    logger.info(f"Annotated {stats['annotated']}/{stats['total']} variants")
    return stats


def validate_vcf(vcf_path: str, detailed: bool = False) -> Dict[str, Any]:
    """
    Validate VCF has required annotations.

    Returns comprehensive validation report.
    """
    report: Dict[str, Any] = {
        "total": 0,
        "GENE": 0,
        "CONSEQUENCE": 0,
        "CADD": 0,
        "multi_allelic": 0,
        "samples": 0,
        "genes": set(),
        "consequences": Counter(),
        "cadd_distribution": {"low": 0, "medium": 0, "high": 0},
    }

    if not Path(vcf_path).exists():
        raise AnnotationError(
            f"VCF file not found: {vcf_path}",
            ["Check the file path is correct"],
        )

    open_func = gzip.open if vcf_path.endswith(".gz") else open

    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                report["samples"] = len(parts) - 9 if len(parts) > 9 else 0
                continue

            report["total"] += 1
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            alt = parts[4]
            if "," in alt:
                report["multi_allelic"] += 1

            info = parts[7]

            # Check for GENE
            if "GENE=" in info:
                report["GENE"] += 1
                # Extract gene name
                for item in info.split(";"):
                    if item.startswith("GENE="):
                        gene = item.split("=")[1]
                        if detailed:
                            report["genes"].add(gene)

            # Check for CONSEQUENCE
            if "CONSEQUENCE=" in info:
                report["CONSEQUENCE"] += 1
                for item in info.split(";"):
                    if item.startswith("CONSEQUENCE="):
                        consequence = item.split("=")[1]
                        if detailed:
                            report["consequences"][consequence] += 1

            # Check for CADD
            if "CADD=" in info:
                report["CADD"] += 1
                for item in info.split(";"):
                    if item.startswith("CADD="):
                        try:
                            cadd = float(item.split("=")[1])
                            if cadd < 15:
                                report["cadd_distribution"]["low"] += 1
                            elif cadd < 25:
                                report["cadd_distribution"]["medium"] += 1
                            else:
                                report["cadd_distribution"]["high"] += 1
                        except ValueError:
                            pass

    # Convert sets to lists for JSON serialization
    report["genes"] = list(report["genes"]) if detailed else len(report["genes"])
    report["consequences"] = dict(report["consequences"]) if detailed else len(report["consequences"])

    return report


def print_validation_report(report: Dict[str, Any], vcf_path: str) -> bool:
    """Print formatted validation report and return pass/fail status."""
    total = max(1, report["total"])

    print(f"\n{'=' * 60}")
    print("VCF Validation Report")
    print(f"{'=' * 60}")
    print(f"File: {vcf_path}")
    print(f"{'=' * 60}")

    print("\nBasic Statistics:")
    print(f"  Total variants: {report['total']:,}")
    print(f"  Samples: {report['samples']}")
    print(f"  Multi-allelic: {report['multi_allelic']:,}")

    print("\nAnnotation Coverage:")
    gene_pct = 100 * report["GENE"] / total
    cons_pct = 100 * report["CONSEQUENCE"] / total
    cadd_pct = 100 * report["CADD"] / total

    gene_status = "PASS" if gene_pct >= 50 else "FAIL"
    cons_status = "PASS" if cons_pct >= 50 else "WARNING"
    cadd_status = "PASS" if cadd_pct >= 30 else "WARNING"

    print(f"  GENE: {report['GENE']:,}/{report['total']:,} ({gene_pct:.1f}%) [{gene_status}]")
    print(f"  CONSEQUENCE: {report['CONSEQUENCE']:,}/{report['total']:,} ({cons_pct:.1f}%) [{cons_status}]")
    print(f"  CADD: {report['CADD']:,}/{report['total']:,} ({cadd_pct:.1f}%) [{cadd_status}]")

    if report["CADD"] > 0:
        print("\nCADD Score Distribution:")
        print(f"  Low (<15): {report['cadd_distribution']['low']:,}")
        print(f"  Medium (15-25): {report['cadd_distribution']['medium']:,}")
        print(f"  High (>25): {report['cadd_distribution']['high']:,}")

    # Overall status
    is_valid = gene_pct >= 50

    print(f"\n{'=' * 60}")
    if is_valid:
        print("Overall Status: PASS")
        print("VCF is ready for pathway subtyping analysis")
    else:
        print("Overall Status: FAIL")
        print("\nRequired actions:")
        if gene_pct < 50:
            print("  1. Add gene annotations using VEP or ANNOVAR")
            print("     VEP: vep -i input.vcf -o output.tsv --tab --symbol")
            print("     ANNOVAR: table_annovar.pl input.avinput humandb/ -buildver hg38 -protocol refGene")
            print("  2. Re-run this script with --vep-tsv or --annovar-file to annotate")
    print(f"{'=' * 60}\n")

    return is_valid


def main():
    parser = argparse.ArgumentParser(
        description="Annotate VCF with GENE, CONSEQUENCE, CADD for pathway subtyping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate existing VCF
  python scripts/annotate_vcf.py input.vcf --validate-only

  # Annotate using VEP output
  python scripts/annotate_vcf.py input.vcf output.vcf --vep-tsv vep_annotations.tsv

  # Annotate using ANNOVAR output
  python scripts/annotate_vcf.py input.vcf output.vcf --annovar-file sample.hg38_multianno.txt

  # Verbose mode for debugging
  python scripts/annotate_vcf.py input.vcf output.vcf --vep-tsv annotations.tsv --verbose

For more information, see: docs/troubleshooting.md#data-issues
        """,
    )
    parser.add_argument("input_vcf", help="Input VCF file (can be .vcf or .vcf.gz)")
    parser.add_argument("output_vcf", nargs="?", help="Output annotated VCF file")
    parser.add_argument("--vep-tsv", help="VEP TSV annotation file (from vep --tab)")
    parser.add_argument("--annovar-file", help="ANNOVAR multianno file (from table_annovar.pl)")
    parser.add_argument("--validate-only", action="store_true", help="Only validate annotations, don't annotate")
    parser.add_argument("--detailed", action="store_true", help="Show detailed validation statistics")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")

    args = parser.parse_args()

    try:
        if args.validate_only:
            report = validate_vcf(args.input_vcf, detailed=args.detailed)
            is_valid = print_validation_report(report, args.input_vcf)
            sys.exit(0 if is_valid else 1)

        if not args.output_vcf:
            parser.error("output_vcf is required when not using --validate-only")

        # Load annotations
        annotations = {}
        if args.vep_tsv:
            annotations = parse_vep_tsv(args.vep_tsv, verbose=args.verbose)
        elif args.annovar_file:
            annotations = parse_annovar_multianno(args.annovar_file, verbose=args.verbose)
        else:
            parser.error("Must provide either --vep-tsv or --annovar-file")

        if len(annotations) == 0:
            raise AnnotationError(
                "No annotations were loaded from the annotation file",
                [
                    "Check that the annotation file is not empty",
                    "Verify the file format matches expected VEP/ANNOVAR output",
                    "Run with --verbose to see detailed parsing statistics",
                ],
            )

        # Annotate VCF
        annotate_vcf(args.input_vcf, args.output_vcf, annotations, verbose=args.verbose)

        print(f"\nAnnotated VCF written to: {args.output_vcf}")

        # Validate output
        print("\nValidating output VCF...")
        report = validate_vcf(args.output_vcf)
        is_valid = print_validation_report(report, args.output_vcf)

        if not is_valid:
            print("WARNING: Output VCF still has low annotation coverage")
            print("Consider using a different annotation source or checking variant ID matching")
            sys.exit(1)

    except AnnotationError as e:
        logger.error(str(e))
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
