#! /usr/bin/env python3

import sys
import math
import glob

from mule.InfoError import *
from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from mule_local.postprocessing.PlaneData_Spectrum import *



class pickle_PlaneData_Spectrum(InfoError):

    def __init__(
            self,
            ref_file_ending = None,
            jobdir_pattern = None,
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

            'only_last_file':
                Process only last file

            other parameters:
                see PlaneData_Spectrum.py
        """

        InfoError.__init__(self, "pickle_PlaneData_Spectrum")
        self._setup(ref_file_ending=ref_file_ending, jobdir_pattern=jobdir_pattern, params=params)



    def _setup(
            self,
            ref_file_ending = None,
            jobdir_pattern = None,
            params = [],
        ):

        if jobdir_pattern == None:
            jobdir_pattern = './job_bench*'

        j = JobsData(jobdir_pattern, verbosity=0)
        jobs = j.get_flattened_data()

        no_reference_job_unique_id_found = True

        if len(jobs) == 0:
            raise Exception("No jobs found!")

        for key, job in jobs.items():
            print("Processing "+key)

            # job directory path
            job_dirpath = job['jobgeneration.job_dirpath']

            # files 
            data_files = glob.glob(job_dirpath+'/output*prog*.csv')

            if len(data_files) == 0:
                self.info("WARNING")
                self.info("WARNING: No files found")
                self.info("WARNING: However, there should be at least one file (the one at t=0)")
                self.info("WARNING")


            # Iterate over all files
            for data_file in data_files:

                # Skip files stored in spectral space
                if '_spec' in data_file:
                    continue

                try:
                    s = PlaneData_Spectrum(
                            data_file,
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

                # determine (time-depending) ending of reference file
                pickle_filename = data_file.replace('.csv', '_spectrum.pickle')

                print("Writing file "+pickle_filename)
                s.write_file(pickle_filename)



if __name__ == '__main__':

    #p = pickle_PlaneData_Spectrum(params=['only_last_file'])
    p = pickle_PlaneData_Spectrum()

