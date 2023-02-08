#! /usr/bin/env python3

import sys
import math
import glob

from mule.InfoError import *
from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.PlaneDataPhysicalDiff import *



class pickle_PlaneDataPhysicalDiff(InfoError):

    def __init__(
            self,
            ref_file_ending = None,
            jobdir_pattern = None,
            job_dirs = None,
            params = [],
        ):
        """
        Generate the .pickle files in each job directory based on simulation output files given by 'ref_file_ending'

        Parameters
        ----------
        ref_file_ending: str
            string with ending of reference files

        jobdir_pattern: str
            string with matching pattern for job directories

        params: list of strings
            list with optional parameters

            'ignore_missing_file':
                Don't throw an error if a file is missing

            ohter parameters:
                see PlaneDataPhysicalDiff.py
        """

        InfoError.__init__(self, "pickle_PlaneDataPhysicalDiff")
        self._setup(ref_file_ending=ref_file_ending, jobdir_pattern=jobdir_pattern, job_dirs=job_dirs, params=params)



    def _setup(
            self,
            ref_file_ending = None,
            jobdir_pattern = None,
            job_dirs = None,
            params = [],
        ):

        if job_dirs != None:
            j = JobsData(job_dirs=job_dirs, verbosity=0)

        else:
            if jobdir_pattern == None:
                jobdir_pattern = './job_bench*'

            j = JobsData(jobdir_pattern, verbosity=0)

        jobs = j.get_flattened_data()


        no_reference_job_unique_id_found = True

        if len(jobs) == 0:
            raise Exception("No jobs found!")

        for key, job in jobs.items():
            print("Processing "+key)

            # Sort out jobs which don't have a reference job id
            # These jobs are likely the reference jobs themselves
            if 'jobgeneration.reference_job_unique_id' not in job:
                continue

            no_reference_job_unique_id_found = False

            reference_job_unique_id = job['jobgeneration.reference_job_unique_id']
            print(" + ref job id: "+reference_job_unique_id)

            ref_key = None
            for skey, sjob in jobs.items():
                if sjob['jobgeneration.job_unique_id'] == reference_job_unique_id:
                    if ref_key != None:
                        raise Exception("FATAL: Reference job already found and now there's another one? This is probably not what you wanted, there might be 2 reference jobs")

                    ref_key = skey

            if ref_key == None:
                print("")
                print("FATAL: missing reference job with id "+reference_job_unique_id)
                print("")
                print("FATAL: reference job was intended for job with dirpath: "+job['jobgeneration.job_dirpath'])
                print("")
                print("FATAL: Hint: If specifying job directories manually, reference job *MUST* be included in the provided job directories!")
                print("")
                raise Exception("Reference job not found!")

            # Load reference job
            ref_job = jobs[ref_key]

            #
            # Load
            #     ref_files:        list of reference files
            #    use_ref_file_ending:    file ending for pickle output file
            #
            # Were the reference filenames provided?
            # ...
            if 'output.reference_filenames' in job:
                # ... then we use the reference files

                # They are available in 'output.reference_filenames' and separated by ';'
                ref_files = job['output.reference_filenames'].split(";")

                #
                # Now we have to find the file ending without the time stamp
                # We do this to generate a unique pickle file which is independent of the time
                #
                # We guess that this starts at the last '_' character in the filename
                # E.g. '_t00000864000.00000000.csv'
                pos = ref_files[0].rfind('_')
                if pos < 0:
                    raise Exception("File ending not found for reference file '"+ref_files[0]+"'")

                use_ref_file_ending = ref_files[0][pos:]
                print("use_ref_file_ending: "+use_ref_file_ending)

                if len(ref_files) == 0:
                    print("Meta tags with name of reference files not found!")
                    print("*"*80)
                    print("Reference directory: "+ref_job['jobgeneration.job_dirpath'])
                    print("Job directory: "+job['jobgeneration.job_dirpath'])
                    print("Reference file endings: "+use_ref_file_ending)
                    print("*"*80)
                    print("* Skipping this job data")
                    print("*"*80)
                    #raise Exception("Meta tags with names of reference files not found!")
                    continue

            else:

                if ref_file_ending != None:
                    use_ref_file_ending = ref_file_ending
                else:
                    print("*"*80)

                    # "output_%s_t%020.8f.csv"
                    use_ref_file_ending = "_t{:020.8f}.csv".format(float(ref_job['runtime.max_simulation_time'])/(60*60))

                if use_ref_file_ending == "":
                    raise Exception("No reference file ending provided / found")

                # Load reference files
                ref_files = []
                files = os.listdir(ref_job['jobgeneration.job_dirpath'])
                for f in files:
                    if use_ref_file_ending in f:
                        ref_files.append(f)


                if len(ref_files) == 0:
                    print("No reference files found or meta tag with reference file names not detected!")
                    print("*"*80)
                    print("Reference directory: "+ref_job['jobgeneration.job_dirpath'])
                    print("Job directory: "+job['jobgeneration.job_dirpath'])
                    print("Reference file endings: "+use_ref_file_ending)
                    print("*"*80)

                    # Search for tag which indicates that simulation was successfully finished
                    if not 'output.simulation_successfully_finished' in job:
                        print("(ignoring error, since simulation was not successfully finished)")
                    else:
                        raise Exception("No reference files not found or meta tag with reference file names not detected!")


            for ref_file in ref_files:
                print("")
                print("Reference file: "+ref_file)

                if '_spec_' in ref_file:
                    self.info("WARNING: Skipping '"+ref_file+"', since this is spectral data")

                else:
                    s = None

                    try:
                        s = PlaneDataPhysicalDiff(
                                ref_job['jobgeneration.job_dirpath']+'/'+ref_file,
                                job['jobgeneration.job_dirpath']+'/'+ref_file,
                                params
                        )

                    except FileNotFoundError as e:
                        # Ignoring missing files should be configured via "ignore_missing_file" parameter, see above
                        if "ignore_missing_file" in params:
                            self.info("Ignoring Error:")
                            self.info(str(e))
                            continue

                        raise Exception(e)

                    except IOError as e:
                        # Ignoring missing files should be configured via "ignore_missing_file" parameter, see above
                        if "ignore_missing_file" in params:
                            self.info("Ignoring Error:")
                            self.info(str(e))
                            continue

                        raise Exception(e)


                    s.print()

                    pickle_filename = 'plane_data_diff_'+ref_file.replace('output_', '').replace(use_ref_file_ending, '')+'.pickle'

                    print("Writing file "+pickle_filename)
                    s.write_file(job['jobgeneration.job_dirpath']+'/'+pickle_filename)

            print(ref_key)
            print("")



        if no_reference_job_unique_id_found:
            print("*"*80)
            print("Warning: No data generated")
            print("No job with a reference_job_unique_id found!")
            print("*"*80)



if __name__ == '__main__':

    p = pickle_PlaneDataPhysicalDiff()

