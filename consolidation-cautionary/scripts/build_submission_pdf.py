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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
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

    cmd = ["pandoc", args.md, "-f", "gfm", "-t", "html5", "-s",
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
