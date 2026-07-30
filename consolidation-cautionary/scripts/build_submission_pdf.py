#!/usr/bin/env python3
"""Build a submission-ready PDF from the manuscript markdown.

Pipeline: markdown --(pandoc)--> standalone HTML with print CSS --(headless
Chrome)--> PDF. Chosen because no LaTeX engine is installed on this host; Chrome
gives reliable pagination and table rendering.

Usage:
    python build_submission_pdf.py --md MANUSCRIPT.md --out OUTPUT.pdf
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

CSS = """
@page { size: A4; margin: 22mm 20mm 20mm 20mm; }
html { -webkit-print-color-adjust: exact; }
body {
  font-family: "Times New Roman", Times, serif;
  font-size: 11.5pt; line-height: 1.55; color: #000; max-width: none;
}
h1 { font-size: 17pt; line-height: 1.25; margin: 0 0 4pt 0; text-align: left; }
h2 {
  font-size: 13pt; margin: 20pt 0 6pt 0; padding-top: 4pt;
  border-top: 0.6pt solid #bbb; page-break-after: avoid;
}
h3 { font-size: 11.8pt; margin: 14pt 0 4pt 0; page-break-after: avoid; }
h2 + p, h3 + p { margin-top: 0; }
p { margin: 0 0 7pt 0; text-align: justify; hyphens: auto; }
em { font-style: italic; }
hr { border: none; border-top: 0.6pt solid #ccc; margin: 14pt 0; }
table {
  border-collapse: collapse; width: 100%; margin: 9pt 0 11pt 0;
  font-size: 9.8pt; page-break-inside: avoid;
}
th, td {
  border: 0.5pt solid #999; padding: 4pt 6pt; text-align: left;
  vertical-align: top;
}
th { background: #eee; font-weight: bold; }
code {
  font-family: "SF Mono", Menlo, Consolas, monospace; font-size: 9.2pt;
  background: #f2f2f2; padding: 0.5pt 2.5pt; border-radius: 2px;
}
blockquote {
  margin: 8pt 0 8pt 12pt; padding-left: 10pt; border-left: 2pt solid #ccc;
  color: #333;
}
ul, ol { margin: 0 0 8pt 0; padding-left: 18pt; }
li { margin-bottom: 3pt; }
/* keep a heading with its first table/paragraph */
h2, h3 { break-after: avoid-page; }
"""


def strip_draft_banner(md_text: str) -> tuple[str, int]:
    """Remove the leading internal draft/review-status blockquote.

    The working manuscript opens with a blockquote recording draft version, hostile-review
    verdicts, PI-gated blockers, "Not yet submitted", and internal section references. That
    block is for the authors; it must never reach a journal or a third-party reviewer.

    This used to be stripped by an ad-hoc `awk` invocation OUTSIDE this script, which meant
    the protection was only as reliable as whoever remembered to pipe through it. It did not
    survive: a PDF built by calling this script directly shipped the whole banner, including
    the private Bitbucket remote and the review history. The strip now lives here, so the
    only way to emit the banner is to ask for it explicitly with --keep-banner.

    Only a blockquote that appears BEFORE the first heading is removed, so a legitimate
    blockquote inside the body is never touched.
    """
    lines = md_text.splitlines(keepends=True)
    i = 0
    while i < len(lines) and not lines[i].lstrip().startswith("#"):
        i += 1
    head = lines[:i]
    if not any(ln.lstrip().startswith(">") for ln in head):
        return md_text, 0
    return "".join(lines[i:]), len(head)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--keep-banner", action="store_true",
                    help="keep the internal draft banner (author preview only; never for "
                         "a journal or third-party reviewer)")
    args = ap.parse_args()

    if not shutil.which("pandoc"):
        sys.exit("pandoc not found")
    if not os.path.exists(CHROME):
        sys.exit(f"Chrome not found at {CHROME}")

    tmp = tempfile.mkdtemp(prefix="pdfbuild_")
    css_path = os.path.join(tmp, "print.css")
    html_path = os.path.join(tmp, "manuscript.html")
    with open(css_path, "w") as fh:
        fh.write(CSS)

    md_src = args.md
    if not args.keep_banner:
        text = open(args.md, encoding="utf-8").read()
        stripped, n = strip_draft_banner(text)
        if n:
            md_src = os.path.join(tmp, os.path.basename(args.md))
            with open(md_src, "w", encoding="utf-8") as fh:
                fh.write(stripped)
            print(f"Stripped internal draft banner ({n} lines) before the first heading.")
        else:
            print("No leading draft banner found; nothing stripped.")

    # --resource-path anchors relative image paths to the MARKDOWN's directory.
    # Without it pandoc resolves them against the working directory, which is why
    # the manuscript previously carried absolute /Users/... paths: they "worked"
    # for the author and rendered as broken images for everyone else (and leaked
    # a local filesystem layout into a submitted PDF).
    cmd = ["pandoc", md_src, "-f", "gfm", "-t", "html5", "-s",
           "--resource-path", os.path.dirname(os.path.abspath(args.md)),
           "--css", css_path, "--embed-resources", "-o", html_path]
    if args.title:
        cmd += ["--metadata", f"title={args.title}"]
    subprocess.run(cmd, check=True)

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    subprocess.run([
        CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer",
        f"--print-to-pdf={out}", f"file://{html_path}",
    ], check=True, capture_output=True, timeout=180)

    if not os.path.exists(out):
        sys.exit("Chrome did not produce a PDF")
    print(f"Wrote {out} ({os.path.getsize(out)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
