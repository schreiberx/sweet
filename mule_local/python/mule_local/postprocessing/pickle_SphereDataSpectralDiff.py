#! /usr/bin/env python3

import sys
import math
import glob

from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from mule_local.postprocessing.SphereDataSpectralDiff import *
from mule.parhelper import *


class pickle_SphereDataSpectralDiff:

    def __init__(
            self,
            job_directories = None,
            ref_file_ending = None,
            jobdir_pattern = None,
        ):
        """
        Generate the .pickle files in each job directory based on simulation output files given by 'ref_file_ending'

        Parameters
        ----------
        ref_file_ending: str
            string with ending of reference files
        jobdir_pattern: str
            pattern to detect job directories
            Default: './job_bench*'
        """
        self._setup(job_directories, ref_file_ending, jobdir_pattern)
        pass





    def _setup(
            self,
            job_directories = None,
            ref_file_ending = None,
            jobdir_pattern = None,
        ):

        if job_directories != None:
            j = JobsData(job_dirs = job_directories, verbosity=0)

        else:
            if jobdir_pattern == None:
                jobdir_pattern = './job_bench*'

            j = JobsData(jobdir_pattern, verbosity=0)


        jobs = j.get_flattened_data()

        no_reference_job_unique_id_found = True

        if len(jobs) == 0:
            raise Exception("No jobs found!")


        def proc_job(key, job):
            print("Processing "+key)

            # Sort out jobs which don't have a reference job id
            # These jobs are likely the reference jobs themselves
            if 'jobgeneration.reference_job_unique_id' not in job:
                return

            no_reference_job_unique_id_found = False

            reference_job_unique_id = job['jobgeneration.reference_job_unique_id']
            print(" + ref job id: "+reference_job_unique_id)

            ref_key = None
            for skey, sjob in jobs.items():
                if sjob['jobgeneration.job_unique_id'] == reference_job_unique_id:
                    ref_key = skey

            if ref_key == None:
                print("Fatal: missing reference job with id "+reference_job_unique_id)
                print("Fatal: reference job was intended for job with dirpath: "+job['jobgeneration.job_dirpath'])
                raise Exception("Reference job not found!")

            # Load reference job
            ref_job = jobs[ref_key]

            if ref_file_ending != None:
                use_ref_file_ending = ref_file_ending
            else:
                # "output_%s_t%020.8f.csv"
                use_ref_file_ending = "_t{:020.8f}.sweet".format(float(ref_job['runtime.max_simulation_time'])/(60*60))

            if use_ref_file_ending == "":
                raise Exception("No reference file ending provided / found")

            # Load reference files
            ref_files = []
            files = os.listdir(ref_job['jobgeneration.job_dirpath'])
            for f in files:
                if use_ref_file_ending in f:
                    ref_files.append(f)

            if len(ref_files) == 0:
                print("No reference files found!")
                print("*"*80)
                print("Reference directory: "+ref_job['jobgeneration.job_dirpath'])
                print("Reference file endings: "+use_ref_file_ending)
                print("*"*80)
                raise Exception("Reference files not found!")

            for ref_file in ref_files:
                s = None
                try:
                    s = SphereDataSpectralDiff(
                            ref_job['jobgeneration.job_dirpath']+'/'+ref_file,
                            job['jobgeneration.job_dirpath']+'/'+ref_file
                    )
                except Exception as e:
                    raise e
                    #print(str(e))
                    #raise Exception("Error occured which is ignored")
                    # Ignore missing files
                    return

                s.print()

                pickle_filename = 'sphere_data_diff_'+ref_file.replace('output_', '').replace(use_ref_file_ending, '')+'.pickle'
                print("Writing file "+pickle_filename)
                s.write_file(job['jobgeneration.job_dirpath']+'/'+pickle_filename)

            print(ref_key)
            print("")


        k = list(jobs.keys())
        v = list(jobs.values())

        n = len(jobs)
        def fun(i):
            try:
                proc_job(k[i], v[i])
            except:
                raise

        if 0:
            # Doesn't really provide any speedup :-(
            parhelper(n, fun, None, use_tqdm=True, verbose=1)

        else:
            for i in range(len(k)):
                fun(i)


        if no_reference_job_unique_id_found:
            print("*"*80)
            print("Warning: No data generated")
            print("No job with a reference_job_unique_id found!")
            print("*"*80)



if __name__ == '__main__':

    if len(sys.argv) > 1:
        p = pickle_SphereDataSpectralDiff(sys.argv[1:])
    else:
        p = pickle_SphereDataSpectralDiff()

