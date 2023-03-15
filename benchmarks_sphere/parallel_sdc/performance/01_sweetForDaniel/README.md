# First string scaling test with space-time parallelization

## Description

_First strong scaling tests of Parallel SDC implementation_

- fixed time step, fixed space mesh resolution
- measurements :
    1. space parallelization only, increase number of OMP threads $[1,2,4,8,16,32]$
    2. time-parallelization with 4 nodes, increase space OMP threads from 1 to 8.

## Running benchmark

```bash
mule.benchmark.cleanup_job_dirs; ./create_jobs.py; ./compile_platform*.sh; mule.benchmark.jobs_run_directly
```