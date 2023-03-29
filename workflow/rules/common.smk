import pandas as pd
import numpy as np

samples = (
    pd.read_csv(config["samples"], dtype={"sample_name": str}, sep="\t")
    .set_index("sample_name", drop=False)
    .sort_index()
)

outdir = config["outdir"]


def is_paired_end(sample):
    sample_units = samples.loc[sample]
    null_units = sample_units.isnull()
    paired = ~null_units["fq2"]
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), f"invalid units for sample {sample}, must be all paired end or all single end"
    return all_paired
