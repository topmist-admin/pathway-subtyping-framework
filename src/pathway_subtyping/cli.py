"""
Command-line interface for the Pathway Subtyping Framework.

Usage:
    python -m pathway_subtyping --config configs/example_autism.yaml
    psf --config configs/example_schizophrenia.yaml
"""

import sys
from pathlib import Path

import click

from . import __version__
from .pipeline import DemoPipeline, PipelineConfig


@click.command()
@click.option(
    "--config",
    "-c",
    type=click.Path(exists=True),
    required=True,
    help="Path to YAML configuration file",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    default=None,
    help="Override output directory from config",
)
@click.option(
    "--seed",
    "-s",
    type=int,
    default=None,
    help="Override random seed from config",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    default=True,
    help="Enable/disable verbose output",
)
@click.version_option(version=__version__, prog_name="pathway-subtyping-framework")
def main(config: str, output: str, seed: int, verbose: bool) -> None:
    """
    Pathway Subtyping Framework - Genetic Subtype Analysis Pipeline

    Run the pathway-based analysis pipeline for molecular subtype discovery.

    Example:
        python -m pathway_subtyping --config configs/example_autism.yaml
    """
    click.echo(f"Pathway Subtyping Framework v{__version__}")
    click.echo("=" * 50)

    # Load configuration
    config_path = Path(config)
    click.echo(f"Loading config: {config_path}")

    try:
        pipeline_config = PipelineConfig.from_yaml(str(config_path))

        # Apply overrides
        if output:
            pipeline_config.output_dir = output
        if seed is not None:
            pipeline_config.seed = seed
        pipeline_config.verbose = verbose

        # Run pipeline
        pipeline = DemoPipeline(pipeline_config)
        pipeline.run()

        click.echo("")
        click.echo("Pipeline completed successfully!")
        click.echo(f"Results: {pipeline_config.output_dir}")

    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Pipeline failed: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
