#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = 'Michael Cuoco'

import os
from snakemake.shell import shell

revision = snakemake.params.get("revision")
profile = snakemake.params.get("profile", [])
params = snakemake.params.get("params-file", )
extra = snakemake.params.get("extra", "")
if isinstance(profile, str):
    profile = [profile]

args = []

if revision:
    args += ["-revision", revision]
if profile:
    args += ["-profile", ",".join(profile)]
if params:
    args += ["-params-file", params]
print(args)

# TODO pass threads in case of single job
# TODO limit parallelism in case of pipeline
# TODO handle other resources

add_parameter = lambda name, value: args.append("--{} {}".format(name, value))

for name, files in snakemake.input.items():
    if isinstance(files, list):
        # TODO check how multiple input files under a single arg are usually passed to nextflow
        files = ",".join(files)
    add_parameter(name, files)
for name, value in snakemake.params.items():
    if (
        name != "pipeline"
        and name != "revision"
        and name != "profile"
        and name != "extra"
        and name != "params-file"
    ):
        add_parameter(name, value)

add_parameter("max-cpus", snakemake.threads)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
args = " ".join(args)
pipeline = snakemake.params.pipeline

shell("nextflow run {pipeline} {args} {extra} {log}")