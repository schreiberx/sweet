#! /usr/bin/env bash

TAGS="prog_h prog_vort prog_div"

for ANALYTICAL in 0 1; do

    JOB_DIFF_DIRS="job_bench_sweet_analytical_$ANALYTICAL job_bench_sweet_analytical_$ANALYTICAL"

    # Directory with reference solution
    JOB_REF_DIR="job_benchref_solution_analytical_$ANALYTICAL"

    for JOB_DIFF_DIR in $JOB_DIFF_DIRS; do
            echo
            echo "Job $JOB_DIFF_DIR"
            echo

            for TAG in $TAGS; do
                    echo "***************************************"
                    echo "Processing $TAG"
                    FILES=$(ls -1 ${JOB_REF_DIR}/*_${TAG}_*.csv)
                    for F in $FILES; do
                            F=$(basename $F)
                            echo "$JOB_REF_DIR/$F"
                            echo "$JOB_DIFF_DIR/$F"
                            echo -n "$F: "
                            if [[ ! -e "$JOB_REF_DIR/$F" ]]; then
                                    echo "ERROR"
                                    exit
                            fi
                            ./cmp_solutions.py "$JOB_REF_DIR/$F" "$JOB_DIFF_DIR/$F" || exit 1
                    done
            done
    done

done
