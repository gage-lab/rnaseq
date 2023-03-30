#!/usr/bin/env python
# Created on: Mar 29, 2023 at 4:42:23 PM
__author__ = "Michael Cuoco"

import gtfparse
import pandas as pd
import sys


def gene_len(gtf):
    gtf["length"] = abs(gtf["end"] - gtf["start"]) + 1
    return gtf.groupby("gene_id").mean("length").reset_index()[["gene_id", "length"]]


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    sys.stderr.write("Reading in files...")

    # Read in gtf files
    genes = gtfparse.read_gtf(snakemake.input.txome_gtf[0])
    rmsk = gtfparse.read_gtf(snakemake.input.rmsk_gtf[0])
    rmsk["gene_id"] = (
        rmsk["transcript_id"]
        + ":"
        + rmsk["gene_id"]
        + ":"
        + rmsk["family_id"]
        + ":"
        + rmsk["class_id"]
    )
    joint = pd.concat(
        [gene_len(genes[genes["feature"] == "transcript"]), gene_len(rmsk)]
    )

    # merge with counts
    counts = pd.read_csv(
        snakemake.input.counts[0],
        sep="\t",
        header=0,
        names=["gene_id", "count"],
        dtype={"gene_id": str, "count": float},
    ).merge(joint, on="gene_id", how="left")

    # compute TPM
    counts["rpk"] = counts["count"] / counts["length"]
    scale_factor = counts["rpk"].sum() / 1e6
    counts["TPM"] = counts["rpk"] / scale_factor
    counts.drop(["rpk"], axis=1, inplace=True)

    # write out
    counts.to_csv(snakemake.output[0], sep="\t", index=False)

    sys.stderr.close()
