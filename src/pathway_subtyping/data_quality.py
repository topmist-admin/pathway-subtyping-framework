"""
Real-World Data Quality Module

Handles messy real-world VCF data with:
- Multi-allelic variant support
- Graceful handling of missing annotations
- Data quality reporting and warnings
- User-friendly error messages with fix guidance

Week 4 deliverable for v0.2 roadmap.
"""

import gzip
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class DataQualityReport:
    """Report on data quality issues found during VCF parsing."""

    total_variants: int = 0
    parsed_variants: int = 0
    skipped_variants: int = 0

    # Annotation coverage
    variants_with_gene: int = 0
    variants_with_consequence: int = 0
    variants_with_cadd: int = 0

    # Multi-allelic handling
    multi_allelic_variants: int = 0
    multi_allelic_expanded: int = 0

    # Genotype issues
    missing_genotypes: int = 0
    malformed_genotypes: int = 0

    # Warnings and errors
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)

    @property
    def gene_coverage(self) -> float:
        """Percentage of variants with gene annotation."""
        if self.parsed_variants == 0:
            return 0.0
        return 100.0 * self.variants_with_gene / self.parsed_variants

    @property
    def consequence_coverage(self) -> float:
        """Percentage of variants with consequence annotation."""
        if self.parsed_variants == 0:
            return 0.0
        return 100.0 * self.variants_with_consequence / self.parsed_variants

    @property
    def cadd_coverage(self) -> float:
        """Percentage of variants with CADD score."""
        if self.parsed_variants == 0:
            return 0.0
        return 100.0 * self.variants_with_cadd / self.parsed_variants

    @property
    def is_usable(self) -> bool:
        """Check if data quality is sufficient for analysis."""
        # Require at least 50% gene annotation coverage
        return self.gene_coverage >= 50.0 and self.parsed_variants > 0

    def add_warning(self, message: str) -> None:
        """Add a warning message."""
        if message not in self.warnings:
            self.warnings.append(message)
            logger.warning(message)

    def add_error(self, message: str) -> None:
        """Add an error message."""
        if message not in self.errors:
            self.errors.append(message)
            logger.error(message)

    def to_dict(self) -> Dict[str, Any]:
        """Convert report to dictionary."""
        return {
            "total_variants": self.total_variants,
            "parsed_variants": self.parsed_variants,
            "skipped_variants": self.skipped_variants,
            "annotation_coverage": {
                "gene": f"{self.gene_coverage:.1f}%",
                "consequence": f"{self.consequence_coverage:.1f}%",
                "cadd": f"{self.cadd_coverage:.1f}%",
            },
            "multi_allelic": {
                "original": self.multi_allelic_variants,
                "expanded": self.multi_allelic_expanded,
            },
            "genotype_issues": {
                "missing": self.missing_genotypes,
                "malformed": self.malformed_genotypes,
            },
            "warnings": self.warnings,
            "errors": self.errors,
            "is_usable": self.is_usable,
        }

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "Data Quality Report",
            "=" * 40,
            f"Total variants: {self.total_variants}",
            f"Parsed variants: {self.parsed_variants}",
            f"Skipped variants: {self.skipped_variants}",
            "",
            "Annotation Coverage:",
            f"  GENE: {self.gene_coverage:.1f}%"
            f" ({self.variants_with_gene}/{self.parsed_variants})",
            f"  CONSEQUENCE: {self.consequence_coverage:.1f}%"
            f" ({self.variants_with_consequence}/{self.parsed_variants})",
            f"  CADD: {self.cadd_coverage:.1f}%"
            f" ({self.variants_with_cadd}/{self.parsed_variants})",
            "",
            f"Multi-allelic variants: {self.multi_allelic_variants}"
            f" (expanded to {self.multi_allelic_expanded})",
            f"Missing genotypes: {self.missing_genotypes}",
            f"Malformed genotypes: {self.malformed_genotypes}",
        ]

        if self.warnings:
            lines.extend(["", "Warnings:"])
            for w in self.warnings:
                lines.append(f"  - {w}")

        if self.errors:
            lines.extend(["", "Errors:"])
            for e in self.errors:
                lines.append(f"  - {e}")

        status = "PASS" if self.is_usable else "FAIL"
        lines.extend(["", f"Data Quality Status: {status}"])

        return "\n".join(lines)


class VCFDataQualityError(Exception):
    """Exception raised when VCF data quality is insufficient."""

    def __init__(self, message: str, report: DataQualityReport, fix_suggestions: List[str]):
        self.message = message
        self.report = report
        self.fix_suggestions = fix_suggestions
        super().__init__(self._format_message())

    def _format_message(self) -> str:
        lines = [
            self.message,
            "",
            "Data Quality Issues:",
        ]

        if self.report.gene_coverage < 50:
            lines.append(
                f"  - Only {self.report.gene_coverage:.1f}% of variants have GENE annotation"
            )

        if self.report.warnings:
            for w in self.report.warnings[:3]:  # Show first 3 warnings
                lines.append(f"  - {w}")

        lines.extend(["", "Suggested Fixes:"])
        for i, suggestion in enumerate(self.fix_suggestions, 1):
            lines.append(f"  {i}. {suggestion}")

        return "\n".join(lines)


def parse_genotype(gt_string: str, target_allele: int = 1) -> Tuple[int, bool]:
    """
    Parse a VCF genotype string and return allele count for a specific allele.

    Supports:
    - Bi-allelic: 0/0, 0/1, 1/1, 0|0, 0|1, 1|1
    - Multi-allelic: 0/2, 1/2, 2/2, etc.
    - Missing: ./., .|.

    Args:
        gt_string: Genotype string (e.g., "0/1", "1|2")
        target_allele: Which alternate allele to count (1 for first ALT, 2 for second, etc.)
                       For bi-allelic sites, this should always be 1.

    Returns:
        Tuple of (allele_count, is_valid)
        allele_count: Count of the target allele (0, 1, or 2 for diploid)
        is_valid: Whether the genotype was successfully parsed

    Note:
        This function counts occurrences of a SPECIFIC allele number, not any non-reference.
        For bi-allelic sites with target_allele=1:
          - 0/0 -> 0 (no ALT alleles)
          - 0/1 -> 1 (one ALT allele)
          - 1/1 -> 2 (two ALT alleles)
        For multi-allelic with target_allele=2:
          - 0/2 -> 1 (one copy of second ALT)
          - 1/2 -> 1 (one copy of second ALT)
          - 2/2 -> 2 (two copies of second ALT)
    """
    if not gt_string:
        return 0, False

    # Extract GT field (first colon-separated value)
    gt = gt_string.split(":")[0]

    # Handle missing genotypes
    if gt in (".", "./.", ".|."):
        return 0, False

    # Parse phased (|) or unphased (/) genotypes
    separator = "|" if "|" in gt else "/"
    try:
        alleles = gt.split(separator)
        if len(alleles) != 2:
            return 0, False

        # Convert to integers, handling "." for missing
        allele_calls = []
        for a in alleles:
            if a == ".":
                return 0, False
            allele_calls.append(int(a))

        # Count occurrences of the specific target allele
        # This is consistent with multi-allelic handling in expand_multi_allelic()
        count = sum(1 for a in allele_calls if a == target_allele)
        return count, True

    except (ValueError, IndexError):
        return 0, False


def parse_info_field(info_string: str) -> Dict[str, Any]:
    """
    Parse VCF INFO field into a dictionary.

    Handles:
    - Key=Value pairs
    - Flag fields (no value)
    - Multiple values (comma-separated)
    - Missing values (.)

    Args:
        info_string: INFO field string

    Returns:
        Dictionary of INFO field values
    """
    info_dict = {}

    if not info_string or info_string == ".":
        return info_dict

    for item in info_string.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            # Handle multiple values
            if "," in val:
                info_dict[key] = val.split(",")
            elif val == ".":
                info_dict[key] = None
            else:
                info_dict[key] = val
        else:
            # Flag field
            info_dict[item] = True

    return info_dict


def expand_multi_allelic(
    chrom: str,
    pos: int,
    vid: str,
    ref: str,
    alts: List[str],
    info_dict: Dict[str, Any],
    sample_genotypes: List[str],
    samples: List[str],
) -> List[Dict[str, Any]]:
    """
    Expand a multi-allelic variant into multiple bi-allelic records.

    Args:
        chrom: Chromosome
        pos: Position
        vid: Variant ID
        ref: Reference allele
        alts: List of alternate alleles
        info_dict: Parsed INFO field
        sample_genotypes: Raw genotype strings for each sample
        samples: Sample IDs

    Returns:
        List of expanded variant records
    """
    expanded = []

    for alt_idx, alt in enumerate(alts):
        alt_num = alt_idx + 1  # ALT alleles are 1-indexed

        # Create variant record for this allele
        variant = {
            "chrom": chrom,
            "pos": pos,
            "id": f"{vid}_{alt_num}" if len(alts) > 1 else vid,
            "ref": ref,
            "alt": alt,
            "original_alts": ",".join(alts) if len(alts) > 1 else None,
        }

        # Extract allele-specific INFO fields
        for key, val in info_dict.items():
            if isinstance(val, list) and len(val) == len(alts):
                # Allele-specific field (e.g., AC, AF)
                variant[key.lower()] = val[alt_idx]
            elif isinstance(val, list):
                # Keep as is
                variant[key.lower()] = val[0] if val else None
            else:
                variant[key.lower()] = val

        # Parse genotypes for this specific allele
        gt_row = {}
        for sample, gt_data in zip(samples, sample_genotypes):
            gt = gt_data.split(":")[0]

            # For multi-allelic, check if this specific allele is present
            if gt in (".", "./.", ".|."):
                gt_row[sample] = 0
                continue

            separator = "|" if "|" in gt else "/"
            try:
                allele_calls = [int(a) if a != "." else 0 for a in gt.split(separator)]
                # Count occurrences of this specific alternate allele
                count = sum(1 for a in allele_calls if a == alt_num)
                gt_row[sample] = count
            except (ValueError, IndexError):
                gt_row[sample] = 0

        variant["genotypes"] = gt_row
        expanded.append(variant)

    return expanded


def load_vcf_with_quality_check(
    vcf_path: str,
    strict: bool = False,
    min_gene_coverage: float = 50.0,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], DataQualityReport]:
    """
    Load VCF file with comprehensive data quality checking.

    Features:
    - Multi-allelic variant expansion
    - Missing annotation handling
    - Data quality reporting
    - User-friendly error messages

    Args:
        vcf_path: Path to VCF file
        strict: If True, raise exception on quality issues
        min_gene_coverage: Minimum percentage of variants with gene annotation

    Returns:
        Tuple of (variants_df, genotypes_df, samples, quality_report)

    Raises:
        VCFDataQualityError: If data quality is insufficient and strict=True
        FileNotFoundError: If VCF file doesn't exist
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(
            f"VCF file not found: {vcf_path}\n\n"
            f"Please check:\n"
            f"  1. The file path is correct\n"
            f"  2. You have read permissions\n"
            f"  3. The file hasn't been moved or deleted"
        )

    report = DataQualityReport()
    variants = []
    samples = []
    all_genotypes = []

    # Determine if gzipped
    open_func = gzip.open if str(vcf_path).endswith(".gz") else open

    try:
        with open_func(vcf_path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # Skip meta-information lines
                if line.startswith("##"):
                    continue

                # Parse header line
                if line.startswith("#CHROM"):
                    parts = line.split("\t")
                    if len(parts) < 10:
                        report.add_error(
                            "VCF header has fewer than 10 columns. "
                            "Expected: #CHROM POS ID REF ALT QUAL FILTER"
                            " INFO FORMAT <samples>"
                        )
                        continue
                    samples = parts[9:]
                    if not samples:
                        report.add_error("No sample columns found in VCF header")
                    continue

                # Parse variant line
                report.total_variants += 1
                parts = line.split("\t")

                if len(parts) < 9:
                    report.skipped_variants += 1
                    report.add_warning(f"Line {line_num}: Skipped - fewer than 9 columns")
                    continue

                chrom, pos, vid, ref, alt, qual, filt, info, fmt = parts[:9]

                # Parse INFO field
                info_dict = parse_info_field(info)

                # Check for multi-allelic
                alts = alt.split(",")
                is_multi_allelic = len(alts) > 1

                if is_multi_allelic:
                    report.multi_allelic_variants += 1

                # Get sample genotypes
                sample_gts = parts[9:] if len(parts) > 9 else []

                # Expand multi-allelic variants
                if is_multi_allelic:
                    expanded_variants = expand_multi_allelic(
                        chrom, int(pos), vid, ref, alts, info_dict, sample_gts, samples
                    )
                    report.multi_allelic_expanded += len(expanded_variants)
                else:
                    # Single allele - create one record
                    gt_row = {}
                    for sample, gt_data in zip(samples, sample_gts):
                        allele_count, is_valid = parse_genotype(gt_data)
                        gt_row[sample] = allele_count
                        if not is_valid and gt_data not in (".", "./.", ".|."):
                            report.malformed_genotypes += 1
                        elif not is_valid:
                            report.missing_genotypes += 1

                    expanded_variants = [
                        {
                            "chrom": chrom,
                            "pos": int(pos),
                            "id": vid,
                            "ref": ref,
                            "alt": alt,
                            "original_alts": None,
                            "genotypes": gt_row,
                            **{k.lower(): v for k, v in info_dict.items()},
                        }
                    ]

                # Process expanded variants
                for var in expanded_variants:
                    report.parsed_variants += 1

                    # Check annotation coverage
                    gene = var.get("gene", "")
                    if gene and gene != "." and gene is not True:
                        report.variants_with_gene += 1

                    consequence = var.get("consequence", "")
                    if consequence and consequence != "." and consequence is not True:
                        report.variants_with_consequence += 1

                    cadd = var.get("cadd")
                    if cadd is not None and cadd != ".":
                        try:
                            float(cadd)
                            report.variants_with_cadd += 1
                        except (ValueError, TypeError):
                            pass

                    # Store genotypes separately
                    gt_row = var.pop("genotypes", {})
                    all_genotypes.append(gt_row)

                    # Clean up variant record
                    variants.append(
                        {
                            "chrom": var.get("chrom", ""),
                            "pos": var.get("pos", 0),
                            "id": var.get("id", "."),
                            "ref": var.get("ref", ""),
                            "alt": var.get("alt", ""),
                            "qual": _safe_float(qual, 0.0),
                            "filter": filt,
                            "gene": _clean_annotation(var.get("gene")),
                            "consequence": _clean_annotation(var.get("consequence")),
                            "cadd": _safe_float(var.get("cadd"), 0.0),
                        }
                    )

    except gzip.BadGzipFile:
        raise VCFDataQualityError(
            "File appears to be corrupted or not a valid gzip file",
            report,
            [
                "Verify the file is not truncated: `gzip -t your_file.vcf.gz`",
                "Re-download or re-compress the file",
                "Use uncompressed VCF if gzip is causing issues",
            ],
        )
    except UnicodeDecodeError:
        raise VCFDataQualityError(
            "File contains invalid characters (encoding issue)",
            report,
            [
                "Ensure the file is UTF-8 encoded",
                "Check if the file is binary (not text VCF)",
                "Try converting with: `iconv -f ISO-8859-1 -t UTF-8 input.vcf > output.vcf`",
            ],
        )

    # Generate warnings based on coverage
    if report.gene_coverage < 50:
        report.add_warning(
            f"Low GENE annotation coverage ({report.gene_coverage:.1f}%). "
            f"Run VEP or ANNOVAR to annotate variants."
        )

    if report.consequence_coverage < 50:
        report.add_warning(
            f"Low CONSEQUENCE annotation coverage ({report.consequence_coverage:.1f}%). "
            f"Variant effect predictions may be incomplete."
        )

    if report.cadd_coverage < 30:
        report.add_warning(
            f"Low CADD score coverage ({report.cadd_coverage:.1f}%). "
            f"Consider adding CADD annotations for better variant prioritization."
        )

    if report.multi_allelic_variants > 0:
        logger.info(
            f"Expanded {report.multi_allelic_variants} multi-allelic variants "
            f"to {report.multi_allelic_expanded} bi-allelic records"
        )

    # Check if data is usable
    if strict and not report.is_usable:
        fix_suggestions = [
            "Annotate your VCF with VEP: " "`vep -i input.vcf -o annotated.vcf --symbol --pick`",
            "Or use ANNOVAR: "
            "`table_annovar.pl input.vcf humandb/ -buildver hg38"
            " -out output -protocol refGene -operation g`",
            "Use the provided annotation helper: "
            "`python scripts/annotate_vcf.py input.vcf output.vcf"
            " --vep-tsv vep_output.tsv`",
            "See troubleshooting guide: docs/troubleshooting.md#data-issues",
        ]
        raise VCFDataQualityError(
            "VCF data quality is insufficient for analysis",
            report,
            fix_suggestions,
        )

    # Create DataFrames
    variants_df = pd.DataFrame(variants)
    genotypes_df = pd.DataFrame(all_genotypes)

    return variants_df, genotypes_df, samples, report


def _safe_float(value: Any, default: float = 0.0) -> float:
    """Safely convert a value to float."""
    if value is None or value == "." or value is True:
        return default
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def _clean_annotation(value: Any) -> str:
    """Clean an annotation value to a string."""
    if value is None or value == "." or value is True or value is False:
        return ""
    if isinstance(value, list):
        return str(value[0]) if value else ""
    return str(value)


def validate_vcf_for_pipeline(
    vcf_path: str,
    verbose: bool = True,
) -> Tuple[bool, DataQualityReport, List[str]]:
    """
    Validate a VCF file before running the pipeline.

    Performs comprehensive checks and provides actionable feedback.

    Args:
        vcf_path: Path to VCF file
        verbose: Print detailed output

    Returns:
        Tuple of (is_valid, report, fix_suggestions)
    """
    fix_suggestions = []

    try:
        _, _, samples, report = load_vcf_with_quality_check(vcf_path, strict=False)
    except FileNotFoundError as e:
        report = DataQualityReport()
        report.add_error(str(e))
        return False, report, ["Verify the file path and permissions"]
    except VCFDataQualityError as e:
        return False, e.report, e.fix_suggestions

    # Check sample count
    if len(samples) < 10:
        report.add_warning(
            f"Only {len(samples)} samples found. "
            f"Clustering may not be meaningful with <10 samples."
        )

    # Check variant count
    if report.parsed_variants < 100:
        report.add_warning(
            f"Only {report.parsed_variants} variants parsed. "
            f"Analysis may be underpowered with <100 variants."
        )

    # Generate fix suggestions based on issues
    if report.gene_coverage < 50:
        fix_suggestions.append("Annotate variants with gene symbols using VEP or ANNOVAR")

    if report.consequence_coverage < 50:
        fix_suggestions.append("Add consequence annotations (missense, frameshift, etc.)")

    if report.cadd_coverage < 30:
        fix_suggestions.append(
            "Add CADD scores for better variant weighting (optional but recommended)"
        )

    is_valid = report.is_usable and not report.errors

    if verbose:
        print(report.summary())
        if fix_suggestions:
            print("\nSuggested improvements:")
            for i, suggestion in enumerate(fix_suggestions, 1):
                print(f"  {i}. {suggestion}")

    return is_valid, report, fix_suggestions
