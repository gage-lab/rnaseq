import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def is_paired_end(sample):
    sample_units = samples.loc[sample]
    null_units = sample_units.isnull()
    paired = ~null_units["fq2"] 
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired

def get_fq(wildcards):
	s = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
	if not is_paired_end(wildcards.sample):
		return {"fq1": f"{s.fq1}"}
	else:
		return {"fq1": f"{s.fq1}", "fq2": f"{s.fq2}"}