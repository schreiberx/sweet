#! /usr/bin/env python3

import glob
from path import Path
import os

benchref = Path('job_benchref_RT_bgalewsky_fsph0_dt00090.00_W-00001_pf_nlev1_pf_nit8_pf_nnod5_SDC_GAUSS_LOBATTO_pf_rk0_M0256/output_prog_vrt_t00000000144.00000000.sweet')
cmd = "../../mule_local/python/mule_local/postprocessing/SphereDataSpectralDiff.py"

for path in glob.iglob('job_bench_*/output_prog_vrt_t00000000144.00000000.sweet', recursive=True):
    full_cmd = f"{cmd} {benchref} {path} > diffs/diff_{Path(path).parent.name}.txt"
    os.system(full_cmd)
