# First string scaling test with space-time parallelization

## Description

_First strong scaling tests of Parallel SDC implementation_

- fixed time step, fixed space mesh resolution ($1024 \times 1024$ grid)
- measurements :
    1. space parallelization only, increase number of OMP threads $[1,2,4,8,16, ..., N_{max}]$ where $ N_{max}$ is the maximum number of CPU on the system, given by the python function `multiprocessing.cpu_count()`.
    2. time-parallelization with 4 nodes, increase space OMP threads from 1 to $N_{max}/2$.

## Running benchmark

After activating `default_llvm` environment in sweet root :

```bash
$ source activate default_llvm
```

Benchmark steps by step :

```bash
$ ./step0_clean.sh      # clean slate
$ ./step1_setup.py      # setup compilation and jobs
$ ./step2_compile.sh    # compile associated sweet program
$ ./step3_run.sh        # run all simulations
$ ./step4_process.py    # post process, generate figures, etc ...
```

Full run :

```bash
$ for s in step*; do ./"$s"; done
```