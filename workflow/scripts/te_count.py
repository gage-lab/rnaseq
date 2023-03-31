#!/usr/bin/env python
# Created on: Mar 30, 2023 at 9:41:07 PM
__author__ = "Michael Cuoco"

from pathlib import Path
import os, sys, json
from snakemake.shell import shell

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.input.meta_info, "r") as f:
    meta_info = json.load(f)

if "F" in meta_info["library_types"][0]:
    strandedness = "forward"
elif "R" in meta_info["library_types"][0]:
    strandedness = "reverse"
elif "U" in meta_info["library_types"][0]:
    strandedness = "no"
else:
    ValueError("No strandedness information found in meta_info.json")

outdir = Path(snakemake.output[0]).parent.__str__()
if not os.path.exists(outdir):
    os.makedirs(outdir)

# set variables for shell command
stem = outdir + "/" + Path(snakemake.output[0]).stem
cmd = "TElocal" if "telocal" in snakemake.rule else "TEcount"

shell(
    "{cmd} "
    "--BAM {snakemake.input.bam} "
    "--GTF {snakemake.input.txome_gtf} "
    "--TE {snakemake.input.rmsk} "
    "--stranded {strandedness} "
    "--mode {snakemake.params.mode} "
    "--project {stem} "
    "--sortByPos "
    " --verbose 3 2> {snakemake.log}"
)

sys.stderr.close()
