#! /usr/bin/env python3

import sys
import math
import glob

from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
from mule_local.postprocessing.SphereDataPhysicalDiff import *
from mule_local.postprocessing.PlaneDataPhysicalDiff import *
from mule.parhelper import *


class PInT_Errors:

    def __init__(
            self,
            job_directories = None,
            pint_type = None,
            geometry = None,
            precomputed_errors = None,
            file_type = None,
            ref_type = None,
            jobdir_pattern = None
        ):
        """
        Compute or read PInT errors (parareal or MGRIT) along iterations

        Parameters
        ----------
        pint_type: str
            'parareal' or 'mgrit'
        geometry: str
            'plane' or 'sphere'
        precomputed_errors: bool
            True: read errors already computed during simulation
            False: read reference and PInT solutions and compute error
        file_type: str
            'csv' or 'bin'
        ref_type: reference solution (given ref or fine solution)
            'ref' or 'fine'
        jobdir_pattern: str
            pattern to detect job directories
            Default: './job_bench*'
        """
        self._setup(job_directories, pint_type, geometry, precomputed_errors, file_type, ref_type, jobdir_pattern)

    def store_err_in_dict(self, err_dict, errL1, errL2, errLinf, var, t, niter):

        ## store in dict
        ## dict to store errors: err[var][t][err_order] = [[0, err_0], [1, err_1], ... , [niter, err_niter]]
        ## dict to store errors: err[var][t][err_order] = val
        if var not in err_dict.keys():
            err_dict[var] = {}
        if t not in err_dict[var].keys():
            err_dict[var][t] = {}
            err_dict[var][t]["errL1"] = errL1
            err_dict[var][t]["errL2"] = errL2
            err_dict[var][t]["errLinf"] = errLinf

        ###err_dict[var][t]["errL1"] = np.vstack((err_dict[var][t]["errL1"], np.array([niter, errL1])))
        ###err_dict[var][t]["errL2"] = np.vstack((err_dict[var][t]["errL2"], np.array([niter, errL2])))
        ###err_dict[var][t]["errLinf"] = np.vstack((err_dict[var][t]["errLinf"], np.array([niter, errLinf])))

        return


    def read_pint_errors_job(self, job):

        ## dict to store errors: err[var][t][err_order] = [[0, err_0], [1, err_1], ... , [niter, err_niter]]
        err = {};

        path = self.jobs[job]['jobgeneration.job_dirpath'];

        ## get list of error files
        list_files = glob.glob(path + "/" + self.pint_type + "_error*.csv")

        ## get iteration corresponding to this simulation
        jj = job.split("_iter")
        niter_job = int(jj[1])

        print ('Reading ' + str(len(list_files)) + ' error files of job ' + path)
        for f in list_files:
            ## f has the form path/[sim_type]_error_[ref_type]_[var_name]_tTTTTTTTTTTT.TTTTTTT_iterIII.csv

            ## skip files containing residual
            if "residual" in f:
                continue;

            ## identify variable
            ff = os.path.basename(f).split("_t0")
            var = ff[0].split("_" + self.ref_type + "_")[1]

            ## identify time
            ff = ff[1].split("_iter")
            t = float(ff[0])

            ## identify iteration
            niter = int(ff[1].split("." + self.file_type)[0])

            ## read only iteration corresponding to this job
            if niter != niter_job:
                continue;

            lines = [line.rstrip() for line in open(f)]
            assert len(lines) == 8

            ## check info

            spl = lines[0].split();
            assert spl[0] == "#BASESOLUTION", spl
            assert spl[1] == self.ref_type, spl

            spl = lines[1].split();
            assert spl[0] == "#VAR", spl
            assert spl[1] == var, (spl, var)

            spl = lines[2].split();
            assert spl[0] == "#ITERATION", spl
            assert int(spl[1]) == niter, (spl, niter)

            spl = lines[3].split();
            assert spl[0] == "#TIMESLICE", spl

            spl = lines[4].split();
            assert spl[0] == "#TIMEFRAMEEND", spl
            assert float(spl[1]) == t, (spl, t)

            ## get errors

            errL1 = -1;
            errL2 = -1;
            errLinf = -1;

            spl = lines[5].split();
            assert spl[0] == "errL1"
            errL1 = float(spl[1]);

            spl = lines[6].split();
            assert spl[0] == "errL2"
            errL2 = float(spl[1]);

            spl = lines[7].split();
            assert spl[0] == "errLinf"
            errLinf = float(spl[1]);

            if errL1 < 0:
                print("ERROR: err_L1 not found in " + path);
                sys.exit();
            if errL2 < 0:
                print("ERROR: err_L2 not found in " + path);
                sys.exit();
            if errLinf < 0:
                print("ERROR: err_Linf not found in " + path);
                sys.exit();


            self.store_err_in_dict(err, errL1, errL2, errLinf, var, t, niter)

        if len(list_files) > 0:
            print(job, path)
            print("err", err)

        return err


    def compute_pint_errors_job(self, job):

        ## dict to store errors: err[var][t][err_order] = [[0, err_0], [1, err_1], ... , [niter, err_niter]]
        err = {};

        path = self.jobs[job]['jobgeneration.job_dirpath'];

        ## get list of solution files
        list_files = glob.glob(path + "/" + self.pint_type + "output_*.csv")

        ## get iteration corresponding to this simulation
        jj = job.split("_iter")
        niter_job = int(jj[0])

        print ('Reading ' + str(len(list_files)) + ' error files of job ' + path)
        for f in list_files:
            ## f has the form path/[sim_type]_error_[ref_type]_[var_name]_tTTTTTTTTTTT.TTTTTTT_iterIII.csv

            ## skip files containing residual
            if "residual" in f:
                continue;

            ## skip files containing computed errors (should not exist!)
            if var[:14] == self.pint_type + "_error":
                continue;

            ## identify variable
            ff = os.path.basename(f).split("_t0")
            var = ff[0]

            ## identify time
            ff = ff[1].split("_iter")
            t = float(fff[0])

            ## identify iteration
            niter = int(ff[1].split("." + self.file_type)[0])

            ## read only iteration corresponding to this job
            if niter != niter_job:
                continue;

            ## compute errors
            if file_type == "csv":
                if self.geometry == "plane":
                    pass
                else:
                    e = SphereDataPhysicalDiff(
                                                  self.ref_job['jobgeneration.job_dirpath']+'/' + ref_file,
                                                  job['jobgeneration.job_dirpath']+'/' + f
                    )
            elif file_type == "sweet":
                pass
                ###sol, m_max, n_max = read_bin_file(f);

            errL1 = e["errL1"]
            errL2 = e["errL2"]
            errLinf = e["errLinf"]

            self.store_err_in_dict(err, errL1, errL2, errLinf, var, t, niter)

        pass


    def get_pint_errors(self):

        self.ref_key = self.get_reference_job()

        err_all_jobs = {}

        if self.precomputed_errors:
            for job in self.jobs.keys():
                ## skip reference simulation
                if job == self.ref_key:
                    continue
                ## skip original parareal jobs (read only copies created per iteration)
                if not "_iter" in job:
                    continue;
                err_all_jobs[job] = self.read_pint_errors_job(job)

                ##pickle_filename = 'sphere_data_diff_'+ref_file.replace('output_', '').replace(use_ref_file_ending, '')+'.pickle'
                pickle_filename = self.pint_type + '_errors.pickle'
                print("Writing file " + self.jobs[job]['jobgeneration.job_dirpath'] + '/' + pickle_filename)
                self.write_file(self.jobs[job]['jobgeneration.job_dirpath']+'/'+pickle_filename, err_all_jobs[job])

        else:
            self.read_reference_files()
            for job in self.jobs.keys():
                if job == self.ref_key:
                    continue
                err_all_jobs[job] = self.compute_pint_errors_job(job)

        return err_all_jobs

    def write_file(
            self,
            picklefile,
            err,
            tagname = None
        ):

        #
        # If picklefile is specified, write norm data to pickle file.
        # This can be later on further postprocessed!
        #
        if picklefile != None:
            import pickle

            if tagname != None:
                tagname += '.'
            else:
                tagname = ''

            pickle_data = {};

            for var in err.keys():

                for t in err[var].keys():

                    pickle_data[tagname +  var + '.t' + str(t) + '.norm_l1'] = err[var][t]['errL1']
                    pickle_data[tagname +  var + '.t' + str(t) + '.norm_l2'] = err[var][t]['errL2']
                    pickle_data[tagname +  var + '.t' + str(t) + '.norm_linf'] = err[var][t]['errLinf']

            print(" + picklefile: "+str(picklefile))

            with open(picklefile, 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(pickle_data, f)

        print("")


    def get_reference_job(self):

        for job in self.jobs.values():

            # Sort out jobs which don't have a reference job id
            # These jobs are likely the reference jobs themselves
            if 'jobgeneration.reference_job_unique_id' not in job:
                continue

            no_reference_job_unique_id_found = False

            reference_job_unique_id = job['jobgeneration.reference_job_unique_id']
            print(" + ref job id: "+reference_job_unique_id)

            ref_key = None
            for skey, sjob in self.jobs.items():
                if sjob['jobgeneration.job_unique_id'] == reference_job_unique_id:
                    ref_key = skey

            if ref_key == None:
                print("Fatal: missing reference job with id "+reference_job_unique_id)
                print("Fatal: reference job was intended for job with dirpath: "+job['jobgeneration.job_dirpath'])
                raise Exception("Reference job not found!")

        return ref_key


    def read_reference_files(self):


        # Load reference job
        ref_job = self.jobs[self.ref_key]

        #####if ref_file_ending != None:
        #####    use_ref_file_ending = ref_file_ending
        #####else:
        #####    # "output_%s_t%020.8f.csv"
        #####    use_ref_file_ending = "_t{:020.8f}.csv".format(float(ref_job['runtime.max_simulation_time'])/(60*60))

        #####if use_ref_file_ending == "":
        #####    raise Exception("No reference file ending provided / found")

        # Load reference files
        ref_files = []
        files = os.listdir(ref_job['jobgeneration.job_dirpath'])
        for f in files:
            ####if use_ref_file_ending in f:
            ####    ref_files.append(f)
            ref_files.append(f)
        if len(ref_files) == 0:
            print("No reference files found!")
            print("*"*80)
            print("Reference directory: "+ref_job['jobgeneration.job_dirpath'])
            print("Reference file endings: "+use_ref_file_ending)
            print("*"*80)
            raise Exception("Reference files not found!")



    def _setup(
            self,
            job_directories = None,
            pint_type = None,
            geometry = None,
            precomputed_errors = None,
            file_type = None,
            ref_type = None,
            jobdir_pattern = None
        ):

        self.pint_type = pint_type
        self.geometry = geometry
        self.precomputed_errors = precomputed_errors
        if file_type == "csv":
            self.file_type = "csv"
        elif file_type == "bin":
            self.file_type = "sweet"
        self.ref_type = ref_type

        if job_directories != None:
            j = JobsData(job_dirs = job_directories, verbosity=0)

        else:
            if jobdir_pattern == None:
                jobdir_pattern = './job_bench*'

            j = JobsData(jobdir_pattern, verbosity=0)

        self.jobs = j.get_flattened_data()

        no_reference_job_unique_id_found = True

        if len(self.jobs) == 0:
            raise Exception("No jobs found!")


if __name__ == '__main__':

    p = PInT_Errors()
