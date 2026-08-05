"""Shared provenance helper for the revision analyses.

Absolute paths are machine-specific and must not reach the deposit: a re-run once wrote a
local scratch path into a deposited JSONL, and hand-scrubbing it did not prevent the next
re-run from putting it back. Scripts call `safe_argv(args)` instead of `vars(args)`.
"""
import os


def safe_argv(a):
    """vars(a) with any path-like value reduced to its basename."""
    return {k: (os.path.basename(v) if isinstance(v, str) and os.sep in v else v)
            for k, v in vars(a).items()}
